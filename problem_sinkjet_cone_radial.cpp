#include <cmath>
#include <fstream>
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include <algorithm>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// point to  MeshBlock(unit that makes up physical grid), real time(double), how much real time each step, multi-dimensional arrays: primitive variables(affected gas), magnetic field, primitive scalars (passive), conservative variables(mass, momentum, energy),conservative scalars
// Since MeshBlock is only a unit, enrollUserExplicitSourceFunction will call SinkJetSource function for every single Meshblock
void SinkJetSource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, const AthenaArray<Real> &prim_scalar, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);

// Params settings for all MeshBlocks
namespace
{
    struct SinkParams
    {
        Real x0, y0, z0; // position of the center of sink particle -
        Real M0;         // initial total mass
        Real r_sink;     // radius
        Real rho_thr;    // density threshold inside the sink particle, when it goes above the threshold, the gas mass transfer from plain grid cells to part of Msink. It means that the inner gas has collapsed enough to be considered as accreted
        Real rho_floor;  // the density after accretion instead of making cells totally empty, so that new gas can flow in
        Real Msink;      // Current sink mass
        Real gsoft;      // gravity softening length: so that when r is very small, gravity force grow smoothly instead of infinitely
        Real G;          // will set to 1 for now

        int jet_switch; // 0 is OFF, 1 is ON
        Real r_noz;     // radius of nozzle
        Real L_noz;     // length of nozzle
        Real v_jet;     // jet speed along z axis
        Real rho_jet;   // density of jet
        Real A_noz;     // nozzle area
        int lock_pdot;  // 1 for fixed momentum
        Real Pdot;      // momentum flux
        Real v_cap;

        int exclude_nozzle_from_accretion; // for debug, default on
        // for data logging
        std::FILE *fp = nullptr; // a pointer that will fill a CSV file only when open it, std:: means it's a standard C++ library namespace

        // these are summed over all meshblocks per timestep
        Real dm_acc_step = 0.0; // total accretion mass
        Real dm_jet_step = 0.0; // total jet mass removed from sink
        Real mgas_step = 0.0;   // total gas mass

        // 6 faces, 12 directions to track mass inflow and outflow
        Real dM_in_x1min = 0, dM_out_x1min = 0;
        Real dM_in_x1max = 0, dM_out_x1max = 0;
        Real dM_in_x2min = 0, dM_out_x2min = 0;
        Real dM_in_x2max = 0, dM_out_x2max = 0;
        Real dM_in_x3min = 0, dM_out_x3min = 0;
        Real dM_in_x3max = 0, dM_out_x3max = 0;

    } S; // S is a single storage box for all the sink+jet settings, like an object

    inline bool in_cone_plus(const Real x, const Real y, const Real z)
    {
        const Real r_perp = std::sqrt((x - S.x0) * (x - S.x0) + (y - S.y0) * (y - S.y0));
        const Real zapex = S.z0 + S.r_sink;
        const Real h = z - zapex;
        if (h > S.L_noz || h < 0)
            return false;
        const Real r_max = h * (S.r_noz / S.L_noz);
        return (r_perp <= r_max);
    }
    inline bool in_cone_minus(const Real x, const Real y, const Real z)
    {
        const Real r_perp = std::sqrt((x - S.x0) * (x - S.x0) + (y - S.y0) * (y - S.y0));
        const Real zapex = S.z0 - S.r_sink;
        const Real h = zapex - z;
        const Real r_max = h * (S.r_noz / S.L_noz);
        if (h > S.L_noz || h < 0)
            return false;
        return (r_perp <= r_max);
    }

    inline void AccumulateBoundaryFlux(MeshBlock *pmb,
                                       const AthenaArray<Real> &prim,
                                       const Real dt)
    {
        auto *pm = pmb->pmy_mesh;
        auto &ms = pm->mesh_size;

        // x1 boundaries
        // inner x1 boundary:
        if (pmb->block_size.x1min == ms.x1min)
        {
            const int i_in = pmb->is;  // first interior
            const int i_gh = i_in - 1; // first ghost

            for (int k = pmb->ks; k <= pmb->ke; ++k)
            {
                for (int j = pmb->js; j <= pmb->je; ++j)
                {

                    // face area
                    const Real A = pmb->pcoord->GetFace1Area(k, j, i_in);

                    // vx
                    const Real vx_face = 0.5 * (prim(IVX, k, j, i_gh) + prim(IVX, k, j, i_in));

                    // upwind density at face
                    const Real rho_face = (vx_face > 0.0)
                                              ? prim(IDN, k, j, i_gh)  // flow in
                                              : prim(IDN, k, j, i_in); // flow out
                    const Real mdot = rho_face * vx_face * A;
                    if (mdot > 0.0)
                        S.dM_in_x1min += mdot * dt; // inflow
                    else
                        S.dM_out_x1min += (-mdot) * dt; // outflow
                }
            }
        }

        // outer x1 boundary:
        if (pmb->block_size.x1max == ms.x1max)
        {
            const int i_in = pmb->ie;  // last interior
            const int i_gh = i_in + 1; // ghost just outside

            for (int k = pmb->ks; k <= pmb->ke; ++k)
            {
                for (int j = pmb->js; j <= pmb->je; ++j)
                {
                    const Real A = pmb->pcoord->GetFace1Area(k, j, i_in + 1);

                    const Real vx_face = 0.5 * (prim(IVX, k, j, i_in) + prim(IVX, k, j, i_gh));

                    const Real rho_face = (vx_face > 0.0)
                                              ? prim(IDN, k, j, i_in)
                                              : prim(IDN, k, j, i_gh); // flipped from x1i

                    const Real mdot = rho_face * vx_face * A;
                    if (mdot < 0.0)
                        S.dM_in_x1max += (-mdot) * dt; // inflow
                    else
                        S.dM_out_x1max += (mdot)*dt; // outflow
                }
            }
        }
        // x2 boundaries
        if (pmb->block_size.x2min == ms.x2min)
        {
            const int j_in = pmb->js;  // first interior
            const int j_gh = j_in - 1; // first ghost
            for (int k = pmb->ks; k <= pmb->ke; ++k)
            {
                for (int i = pmb->is; i <= pmb->ie; ++i)
                {
                    const Real A = pmb->pcoord->GetFace2Area(k, j_in, i);
                    const Real vy_face = 0.5 * (prim(IVY, k, j_gh, i) + prim(IVY, k, j_in, i));
                    const Real rho_face = (vy_face > 0.0)
                                              ? prim(IDN, k, j_gh, i)  // flow in
                                              : prim(IDN, k, j_in, i); // flow out
                    const Real mdot = rho_face * vy_face * A;
                    if (mdot > 0.0)
                        S.dM_in_x2min += mdot * dt; // inflow
                    else
                        S.dM_out_x2min += (-mdot) * dt; // outflow
                }
            }
        }

        if (pmb->block_size.x2max == ms.x2max)
        {
            const int j_in = pmb->je;
            const int j_gh = j_in + 1;

            for (int k = pmb->ks; k <= pmb->ke; ++k)
            {
                for (int i = pmb->is; i <= pmb->ie; ++i)
                {

                    const Real A = pmb->pcoord->GetFace2Area(k, j_gh, i);
                    const Real vy_face = 0.5 * (prim(IVY, k, j_in, i) + prim(IVY, k, j_gh, i));
                    const Real rho_face = (vy_face > 0.0)
                                              ? prim(IDN, k, j_in, i)  // flow out
                                              : prim(IDN, k, j_gh, i); // flow in
                    const Real mdot = rho_face * vy_face * A;
                    if (mdot < 0.0)
                        S.dM_in_x2max += (-mdot) * dt; // inflow
                    else
                        S.dM_out_x2max += (mdot)*dt; // outflow
                }
            }
        }

        // x3 boundaries
        if (pmb->block_size.x3min == ms.x3min)
        {
            const int k_in = pmb->ks;
            const int k_gh = k_in - 1;
            for (int j = pmb->js; j <= pmb->je; ++j)
            {
                for (int i = pmb->is; i <= pmb->ie; ++i)
                {
                    const Real A = pmb->pcoord->GetFace3Area(k_in, j, i);
                    const Real vz_face = 0.5 * (prim(IVZ, k_gh, j, i) + prim(IVZ, k_in, j, i));
                    const Real rho_face = (vz_face > 0.0)
                                              ? prim(IDN, k_gh, j, i)  // flow in
                                              : prim(IDN, k_in, j, i); // flow out
                    const Real mdot = rho_face * vz_face * A;
                    if (mdot > 0.0)
                        S.dM_in_x3min += mdot * dt; // inflow
                    else
                        S.dM_out_x3min += (-mdot) * dt; // outflow
                }
            }
        }

        if (pmb->block_size.x3max == ms.x3max)
        {
            const int k_in = pmb->ke;
            const int k_gh = k_in + 1;
            for (int j = pmb->js; j <= pmb->je; ++j)
            {
                for (int i = pmb->is; i <= pmb->ie; ++i)
                {
                    const Real A = pmb->pcoord->GetFace3Area(k_gh, j, i);
                    const Real vz_face = 0.5 * (prim(IVZ, k_in, j, i) + prim(IVZ, k_gh, j, i));
                    const Real rho_face = (vz_face > 0.0)
                                              ? prim(IDN, k_in, j, i)  // flow in
                                              : prim(IDN, k_gh, j, i); // flow out
                    const Real mdot = rho_face * vz_face * A;
                    if (mdot < 0.0)
                        S.dM_in_x3max += (-mdot) * dt; // inflow
                    else
                        S.dM_out_x3max += mdot * dt; // outflow
                }
            }
        }
    }

} // namespace ends

void DiodeInnerX1(MeshBlock *, Coordinates *, AthenaArray<Real> &, FaceField &,
                  Real, Real, int, int, int, int, int, int, int);
void DiodeOuterX1(MeshBlock *, Coordinates *, AthenaArray<Real> &, FaceField &,
                  Real, Real, int, int, int, int, int, int, int);
void DiodeInnerX2(MeshBlock *, Coordinates *, AthenaArray<Real> &, FaceField &,
                  Real, Real, int, int, int, int, int, int, int);
void DiodeOuterX2(MeshBlock *, Coordinates *, AthenaArray<Real> &, FaceField &,
                  Real, Real, int, int, int, int, int, int, int);
void DiodeInnerX3(MeshBlock *, Coordinates *, AthenaArray<Real> &, FaceField &,
                  Real, Real, int, int, int, int, int, int, int);
void DiodeOuterX3(MeshBlock *, Coordinates *, AthenaArray<Real> &, FaceField &,
                  Real, Real, int, int, int, int, int, int, int);

// Function: read num values for params from athinput.sinkjet, store each value of S object, carry it around for the whole simulation
void Mesh::InitUserMeshData(ParameterInput *pin)
{ // a pointer parameter that allow us to read the data from athinput.sinkjet
    // all sink params
    S.x0 = pin->GetReal("problem", "sink_x0");
    S.y0 = pin->GetReal("problem", "sink_y0");
    S.z0 = pin->GetReal("problem", "sink_z0");
    S.r_sink = pin->GetReal("problem", "r_sink");
    S.rho_thr = pin->GetReal("problem", "rho_thr");
    S.rho_floor = pin->GetReal("problem", "rho_floor");
    S.Msink = pin->GetReal("problem", "Msink0");
    S.gsoft = pin->GetReal("problem", "g_soft");
    S.G = pin->GetReal("problem", "Gconst");
    // for jets related specifically
    S.jet_switch = pin->GetInteger("problem", "jet_switch");

    std::cout << "DEBUG: jet_switch = " << S.jet_switch << std::endl;

    if (S.jet_switch)
    {
        S.r_noz = pin->GetReal("problem", "r_noz");
        S.L_noz = pin->GetReal("problem", "L_noz");
        S.v_jet = pin->GetReal("problem", "v_jet");
        S.lock_pdot = pin->GetOrAddInteger("problem", "lock_pdot", 1);
        S.rho_jet = pin->GetReal("problem", "rho_jet");
        S.v_cap = pin->GetOrAddReal("problem", "v_cap", 3.0);
        S.A_noz = M_PI * S.r_noz * S.r_noz;       // area per lobe
        S.Pdot = pin->GetReal("problem", "Pdot"); // S.Pdot = 2* S.rho_jet * S.v_jet^2 * S.A_noz;
    }
    // fixed momentum

    S.exclude_nozzle_from_accretion = 1; // avoid the sink accidentally accrete the gas we just injected as a jet
    S.M0 = -1.0;                         // set it later, bc it's impossible to find M0 = Msink +Mgas, without knowing Mgas(which is summed per MeshBlock inside SinkJetSource, requires integrating rho*V over every cell in every MeshBlock (and every MPI rank))

    // run my SinkJetSource function each timestep
    EnrollUserExplicitSourceFunction(SinkJetSource);
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, DiodeInnerX1);
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DiodeOuterX1);
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiodeInnerX2);
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, DiodeOuterX2);
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiodeInnerX3);
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiodeOuterX3);
}

// How to initialize gas on the grid
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    // // // // // // // // // // // // // // // // // // // //
    const Real rho0 = pin->GetReal("problem", "rho0");
    for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
            for (int i = is; i <= ie; ++i)
            {
                phydro->u(IDN, k, j, i) = rho0; // density
                phydro->u(IM1, k, j, i) = 0.0;  // x-p :set momentum to 0 means that velocity is also set to 0, bc v = p/density
                phydro->u(IM2, k, j, i) = 0.0;  // y-p
                phydro->u(IM3, k, j, i) = 0.0;  // z-p
            }
    // This loop goes over every cell in the grid to initialize rho to rho0 in the input file
    //  carve out the sink particle, set every cells to rho_floor
    for (int k = ks; k <= ke; ++k)
    {
        const Real zc = pcoord->x3v(k);
        for (int j = js; j <= je; ++j)
        {
            const Real yc = pcoord->x2v(j);
            for (int i = is; i <= ie; ++i)
            {
                const Real xc = pcoord->x1v(i);
                const Real dx = xc - S.x0, dy = yc - S.y0, dz = zc - S.z0; // compute cell's distance to the sink center
                const Real r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r < S.r_sink)
                { // check if inside the sink sphere
                    phydro->u(IDN, k, j, i) = S.rho_floor;
                    phydro->u(IM1, k, j, i) = 0.0;
                    phydro->u(IM2, k, j, i) = 0.0;
                    phydro->u(IM3, k, j, i) = 0.0;
                }
            }
        }
    }
}

void Mesh::UserWorkInLoop()
{
    std::cout << "DEBUG UserWorkInLoop time=" << time
              << " rank=" << Globals::my_rank << std::endl;
    Real dm_acc = S.dm_acc_step;
    Real dm_jet = S.dm_jet_step;
    Real mgas = S.mgas_step;

// This retrieve values and then combine all
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &dm_acc, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &dm_jet, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &mgas, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

    Real in_x1min = S.dM_in_x1min, out_x1min = S.dM_out_x1min;
    Real in_x1max = S.dM_in_x1max, out_x1max = S.dM_out_x1max;
    Real in_x2min = S.dM_in_x2min, out_x2min = S.dM_out_x2min;
    Real in_x2max = S.dM_in_x2max, out_x2max = S.dM_out_x2max;
    Real in_x3min = S.dM_in_x3min, out_x3min = S.dM_out_x3min;
    Real in_x3max = S.dM_in_x3max, out_x3max = S.dM_out_x3max;

#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &in_x1min, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &out_x1min, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &in_x1max, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &out_x1max, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &in_x2min, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &out_x2min, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &in_x2max, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &out_x2max, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &in_x3min, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &out_x3min, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &in_x3max, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &out_x3max, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

    if (Globals::my_rank == 0)
    {
        static bool wrote_header = false;
        std::ofstream fb("boundary_flux.csv", std::ios::app);
        fb.setf(std::ios::scientific);
        fb.precision(16);

        if (!wrote_header)
        {
            fb << "# time,"
               << "in_x1min,out_x1min,in_x1max,out_x1max,"
               << "in_x2min,out_x2min,in_x2max,out_x2max,"
               << "in_x3min,out_x3min,in_x3max,out_x3max\n";
            wrote_header = true;
        }

        fb << time << ","
           << in_x1min << "," << out_x1min << ","
           << in_x1max << "," << out_x1max << ","
           << in_x2min << "," << out_x2min << ","
           << in_x2max << "," << out_x2max << ","
           << in_x3min << "," << out_x3min << ","
           << in_x3max << "," << out_x3max << "\n";
    }
    S.dM_in_x1min = S.dM_out_x1min = 0;
    S.dM_in_x1max = S.dM_out_x1max = 0;
    S.dM_in_x2min = S.dM_out_x2min = 0;
    S.dM_in_x2max = S.dM_out_x2max = 0;
    S.dM_in_x3min = S.dM_out_x3min = 0;
    S.dM_in_x3max = S.dM_out_x3max = 0;

    S.Msink += (dm_acc - dm_jet);

    // set M0
    if (S.M0 < 0.0)
        S.M0 = S.Msink + mgas;

    // reset
    S.dm_acc_step = 0.0;
    S.dm_jet_step = 0.0;
    S.mgas_step = 0.0;

    // write CSV once
    if (Globals::my_rank == 0)
    {
        const Real sfe_M0 = S.Msink / std::max(S.M0, (Real)1e-60);
        const Real sfe_inst = S.Msink / std::max(S.Msink + mgas, (Real)1e-60);

        std::ofstream f("sinkjet_history.csv", std::ios::app);
        f.setf(std::ios::scientific);
        f.precision(16);
        f << time << "," << S.Msink << "," << mgas << "," << sfe_M0 << "," << sfe_inst << "\n";
    }
}

static inline void DiodeBCImpl(
    BoundaryFace face,
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt,
    int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    (void)pmb;
    (void)pco;
    (void)b;
    (void)time;
    (void)dt;

    const int vel =
        (face == BoundaryFace::inner_x1 || face == BoundaryFace::outer_x1) ? IVX : (face == BoundaryFace::inner_x2 || face == BoundaryFace::outer_x2) ? IVY
                                                                                                                                                      : IVZ;

    const Real inflow_sign =
        (face == BoundaryFace::inner_x1 || face == BoundaryFace::inner_x2 || face == BoundaryFace::inner_x3)
            ? +1.0
            : -1.0;

    auto copy_and_clamp = [&](int k, int j, int i, int kk, int jj, int ii)
    {
        prim(IDN, k, j, i) = prim(IDN, kk, jj, ii);
        prim(IVX, k, j, i) = prim(IVX, kk, jj, ii);
        prim(IVY, k, j, i) = prim(IVY, kk, jj, ii);
        prim(IVZ, k, j, i) = prim(IVZ, kk, jj, ii);
        if (inflow_sign * prim(vel, k, j, i) > 0.0)
            prim(vel, k, j, i) = 0.0;
    };

    if (face == BoundaryFace::inner_x1)
    {
        for (int k = kl; k <= ku; ++k)
            for (int j = jl; j <= ju; ++j)
                for (int g = 1; g <= ngh; ++g)
                    copy_and_clamp(k, j, il - g, k, j, il);
    }
    else if (face == BoundaryFace::outer_x1)
    {
        for (int k = kl; k <= ku; ++k)
            for (int j = jl; j <= ju; ++j)
                for (int g = 1; g <= ngh; ++g)
                    copy_and_clamp(k, j, iu + g, k, j, iu);
    }
    else if (face == BoundaryFace::inner_x2)
    {
        for (int k = kl; k <= ku; ++k)
            for (int i = il; i <= iu; ++i)
                for (int g = 1; g <= ngh; ++g)
                    copy_and_clamp(k, jl - g, i, k, jl, i);
    }
    else if (face == BoundaryFace::outer_x2)
    {
        for (int k = kl; k <= ku; ++k)
            for (int i = il; i <= iu; ++i)
                for (int g = 1; g <= ngh; ++g)
                    copy_and_clamp(k, ju + g, i, k, ju, i);
    }
    else if (face == BoundaryFace::inner_x3)
    {
        for (int j = jl; j <= ju; ++j)
            for (int i = il; i <= iu; ++i)
                for (int g = 1; g <= ngh; ++g)
                    copy_and_clamp(kl - g, j, i, kl, j, i);
    }
    else if (face == BoundaryFace::outer_x3)
    {
        for (int j = jl; j <= ju; ++j)
            for (int i = il; i <= iu; ++i)
                for (int g = 1; g <= ngh; ++g)
                    copy_and_clamp(ku + g, j, i, ku, j, i);
    }
}

void DiodeInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    DiodeBCImpl(BoundaryFace::inner_x1, pmb, pco, prim, b, time, dt, il, iu, jl, ju, kl, ku, ngh);
}

void DiodeOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    DiodeBCImpl(BoundaryFace::outer_x1, pmb, pco, prim, b, time, dt, il, iu, jl, ju, kl, ku, ngh);
}
void DiodeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    DiodeBCImpl(BoundaryFace::inner_x2, pmb, pco, prim, b, time, dt, il, iu, jl, ju, kl, ku, ngh);
}
void DiodeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    DiodeBCImpl(BoundaryFace::outer_x2, pmb, pco, prim, b, time, dt, il, iu, jl, ju, kl, ku, ngh);
}
void DiodeInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    DiodeBCImpl(BoundaryFace::inner_x3, pmb, pco, prim, b, time, dt, il, iu, jl, ju, kl, ku, ngh);
}
void DiodeOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
    DiodeBCImpl(BoundaryFace::outer_x3, pmb, pco, prim, b, time, dt, il, iu, jl, ju, kl, ku, ngh);
}

void SinkJetSource(MeshBlock *pmb, const Real time, const Real dt,
                   const AthenaArray<Real> &prim,
                   const AthenaArray<Real> &bcc,
                   const AthenaArray<Real> &prim_scalar,
                   AthenaArray<Real> &cons,
                   AthenaArray<Real> &cons_scalar)
{
    (void)bcc;                    // not used (no MHD) Three useless parameters bro
    (void)prim_scalar;            // not used (NSCALARS=0)
    (void)cons_scalar;            // not used (NSCALARS=0)
    auto &x1v = pmb->pcoord->x1v; // make it easier to extract x,y,z coordinates of each cell
    auto &x2v = pmb->pcoord->x2v;
    auto &x3v = pmb->pcoord->x3v;
    Real dMsink_local = 0.0;
    Real dMjet_local = 0.0;
    static Real last_time_sum = -1.0;

    if (time != last_time_sum)
    {
        last_time_sum = time;
    }
    // Accretion: remove gas inside sink region above threshold-------------------------------------------------------
    // loop over every cells in the grid excluding the jet nozzle cells so that we don't immediately reaccrete the gas we just injected
    for (int k = pmb->ks; k <= pmb->ke; ++k)
    {
        const Real zc = x3v(k);
        for (int j = pmb->js; j <= pmb->je; ++j)
        {
            const Real yc = x2v(j);
            for (int i = pmb->is; i <= pmb->ie; ++i)
            {
                const Real xc = x1v(i);

                // nozzle mask (bipolar cylinders)
                bool in_noz = false; // we have to check if the cell is inside the sink particle
                if (S.exclude_nozzle_from_accretion && S.jet_switch)
                {
                    in_noz = in_cone_plus(xc, yc, zc) || in_cone_minus(xc, yc, zc); // both nozzle
                }
                if (in_noz)
                    continue; // don't accrete nozzle cells

                const Real dx = xc - S.x0, dy = yc - S.y0, dz = zc - S.z0;
                const Real r = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (r >= S.r_sink)
                    continue; // if outside sink particle sphere do nothing

                const Real rho = cons(IDN, k, j, i);
                if (rho <= S.rho_thr)
                    continue; // below density threshold do nothing

                Real d_rho = rho - S.rho_floor; // Compute how much to remove
                d_rho = std::min(d_rho, 0.9 * rho);
                if (d_rho <= 0.0)
                    continue; // if density less or equal to 0, do nothing, if bigger, then good, proceed to actually remove

                // Msink += removed_gas_mass
                const Real Vcell = pmb->pcoord->GetCellVolume(k, j, i); // volume is read only, it's also k, j, i, not ijk, for Athena array reading, which i think is really wiered
                // DEBUG Statement
                const Real dM = d_rho * Vcell;
                dMsink_local += dM;
                // remove momentum of the removed gas
                cons(IDN, k, j, i) -= d_rho;
                cons(IM1, k, j, i) -= d_rho * prim(IVX, k, j, i);
                cons(IM2, k, j, i) -= d_rho * prim(IVY, k, j, i);
                cons(IM3, k, j, i) -= d_rho * prim(IVZ, k, j, i);
                // physics: velocity of remaining gas is conserved, since gas mass is removed, gas momentum is not conserved
            }
        }
    }

    // Gravity: the sink particle pulls on the surrounding gas, so give each gas cell a small momentum toward the sink ----------------------------------------------------
    // a = - GM r_hat / (r^2 + eps^2)^(3/2) ;;; update momentum = rho * a * dt

    for (int k = pmb->ks; k <= pmb->ke; ++k)
    { // loop over each cell
        const Real zc = x3v(k);
        for (int j = pmb->js; j <= pmb->je; ++j)
        {
            const Real yc = x2v(j);
            for (int i = pmb->is; i <= pmb->ie; ++i)
            {
                const Real xc = x1v(i);
                const Real dx = xc - S.x0, dy = yc - S.y0, dz = zc - S.z0; // displacement from cell to sink particle
                const Real r2_raw = dx * dx + dy * dy + dz * dz;
                if (r2_raw < S.r_sink * S.r_sink)
                    continue;
                const Real r2 = r2_raw + S.gsoft * S.gsoft;
                // gravitational softening to prevent gravity blow up near the center gsoft = r^2 + eps^2
                const Real r3 = r2 * sqrt(r2);
                // find acc  direction towards the sink
                const Real gm_r3 = (S.G * S.Msink) / r3; // compute scalar prefactor
                const Real ax = -gm_r3 * dx;             // ax
                const Real ay = -gm_r3 * dy;             // ay
                const Real az = -gm_r3 * dz;             // az
                // convert acc to momentum change
                const Real rho = cons(IDN, k, j, i);
                cons(IM1, k, j, i) += rho * ax * dt; // Delta(\rho*v) = \rho * \a * dt
                cons(IM2, k, j, i) += rho * ay * dt;
                cons(IM3, k, j, i) += rho * az * dt; // densor gas gets bigger momentum change towards the sink particle
            }
        }
    }

    // Jets: inject mass & momentum in cells from two short cylinders(the nozzle) ---------------------------------------------------------
    // Distribute the mass rate uniformly across all nozzle cells by volume
    /*
    definitions:
      dM_need = How much mass each cell needs
      M_need = How much total mass is needed to bring the nozzle region density to rho_jet
      M_injected = How much mass we ACTUALLY inject, to cap from removing more mass than the sink have
      dM = How much mass each cell receives
      */
    if (S.jet_switch)
    {
        // First: find the total nozzle volume of both lobes
        Real V_nozzle_total = 0.0;
        Real M_need = 0.0;
        Real M_nozzle_total = 0.0;
        Real rhoj = std::max(S.rho_jet, 1e-12);
        Real v_eff;
        // Real Px_inj = 0, Py_inj = 0, Pz_inj = 0;
        for (int k = pmb->ks; k <= pmb->ke; ++k)
        {
            const Real zc = x3v(k);
            for (int j = pmb->js; j <= pmb->je; ++j)
            {
                const Real yc = x2v(j);
                for (int i = pmb->is; i <= pmb->ie; ++i)
                {
                    const Real xc = x1v(i);
                    if (in_cone_plus(xc, yc, zc) || in_cone_minus(xc, yc, zc))
                    {
                        const Real dV = pmb->pcoord->GetCellVolume(k, j, i);
                        V_nozzle_total += dV;
                        const Real dM_need = std::max(rhoj - cons(IDN, k, j, i), 0.0) * dV;
                        M_need += dM_need;
                        M_nozzle_total += cons(IDN, k, j, i) * dV;
                    }
                }
            }
        }
        if (Globals::my_rank == 0)
        {
            std::cout << "DEBUG nozzle V=" << V_nozzle_total
                      << " M=" << M_nozzle_total
                      << "M_need =" << M_need << "\n";
        }
        const Real M_injected = std::min(M_need, std::max(S.Msink, 0.0));
        if (M_need <= 0.0)
            return;
        if (V_nozzle_total > 0.0)
        {
            // total mass injection rate = rho_jet * v_jet * A_face * 2 (two lobes)
            if (S.lock_pdot)
            {
                const Real P_total = S.Pdot * dt; // total, S.Pdot is for total not per lobe
                // calc v_eff so that momentum flux matches Pdot
                v_eff = P_total / M_injected;
                // v_eff = std::min(v_eff, S.v_cap); // cap the speed
            }
            else
            {
                v_eff = S.v_jet;
            }
            if (Globals::my_rank == 0)
            {
                std::cout << "DEBUG target P_total=" << (S.Pdot * dt)
                          << " v_eff=" << v_eff << "\n";
            }
            // Second: Distribute mass and momentum to each nozzle cell porportional to their volume; assign +z momentum in +lobe and -z in -lobe
            for (int k = pmb->ks; k <= pmb->ke; ++k)
            {
                const Real zc = x3v(k);
                for (int j = pmb->js; j <= pmb->je; ++j)
                {
                    const Real yc = x2v(j);
                    for (int i = pmb->is; i <= pmb->ie; ++i)
                    {
                        const Real xc = x1v(i);
                        const bool plus = in_cone_plus(xc, yc, zc);              // in the upper lobe
                        const bool minus = (!plus) && in_cone_minus(xc, yc, zc); // in the downer lobe
                        if (!(plus || minus))
                            continue;                                                 // do nothing if not in nozzle
                        const Real z0 = plus ? (S.z0 + S.r_sink) : (S.z0 - S.r_sink); // determine zapex depending on in which lobe to give it radial component velocity
                        const Real rx = xc - S.x0;
                        const Real ry = yc - S.y0;
                        const Real rz = zc - z0;
                        const Real vrmag = std::sqrt(rx * rx + ry * ry + rz * rz); // v radial magnitude
                        if (vrmag < 1e-12)
                            continue; // if too small, avoid dividing, continue = return
                        const Real ux = rx / vrmag;
                        const Real uy = ry / vrmag; // UNIT vectors
                        const Real uz = rz / vrmag;

                        const Real dV = pmb->pcoord->GetCellVolume(k, j, i);
                        const Real dM = std::max(rhoj - cons(IDN, k, j, i), 0.0) * dV / M_need * M_injected; // find the inject jet mass distributed per cell
                        const Real dRho = dM / dV;                                                           // compute density increment
                        // add to current density
                        cons(IDN, k, j, i) += dRho;

                        // calculate the injected momentum per cell
                        const Real P = dM * v_eff; // only magnitude
                        cons(IM1, k, j, i) += (P * ux) / dV;
                        cons(IM2, k, j, i) += (P * uy) / dV;
                        cons(IM3, k, j, i) += (P * uz) / dV;
                    }
                }
            }
            dMjet_local += M_injected;
        }
    }
    Real Mgas_local = 0.0;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
    {
        for (int j = pmb->js; j <= pmb->je; ++j)
        {
            for (int i = pmb->is; i <= pmb->ie; ++i)
            {
                const Real V = pmb->pcoord->GetCellVolume(k, j, i);
                Mgas_local += cons(IDN, k, j, i) * V;
            }
        }
    }
    S.mgas_step += Mgas_local;
    S.dm_acc_step += dMsink_local;
    S.dm_jet_step += dMjet_local;
    AccumulateBoundaryFlux(pmb, prim, dt);
}
