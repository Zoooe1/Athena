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

        int jet_switch;   // 0 is OFF, 1 is ON
        Real r_noz;       // radius of nozzle
        Real L_noz;       // length of nozzle
        Real v_jet;       // jet speed along z axis
        Real rho_jet;     // density of jet
        Real A_noz;       // nozzle area
        int lock_pdot;    // 1 for fixed momentum
        Real Pdot_target; // value of momentum
        Real v_cap;

        int exclude_nozzle_from_accretion; // for debug, default on
        // for data logging
        std::FILE *fp = nullptr; // a pointer that will fill a CSV file only when open it, std:: means it's a standard C++ library namespace

        // these are summed over all meshblocks per timestep
        Real dm_acc_step = 0.0; // total accretion mass
        Real dm_jet_step = 0.0; // total jet mass removed from sink
        Real mgas_step = 0.0;   // total gas mass
    } S;                        // S is a single storage box for all the sink+jet settings, like an object

    // Check if a given cell is within the nozzle on the top of the sink particle sphere
    inline bool in_nozzle_plus(const Real x, const Real y, const Real z)
    {
        // +z lobe: axial window [z0 + r_sink, z0 + r_sink + L_noz]
        const Real r_perp = std::sqrt((x - S.x0) * (x - S.x0) + (y - S.y0) * (y - S.y0)); // define the cell's distance to the center z axis
        if (r_perp > S.r_noz)
            return false; // if lie within Z-axis, continue to check if within the length of the cylinder
        const Real zmin = S.z0 + S.r_sink;
        const Real zmax = S.z0 + S.r_sink + S.L_noz;
        return (z >= zmin && z <= zmax);
    }
    // Check if a given cell is within the nozzle at the bottom of the sink particle sphere
    inline bool in_nozzle_minus(const Real x, const Real y, const Real z)
    {
        // -z lobe: axial window [z0 - r_sink - L_noz, z0 - r_sink]
        const Real r_perp = std::sqrt((x - S.x0) * (x - S.x0) + (y - S.y0) * (y - S.y0));
        if (r_perp > S.r_noz)
            return false;
        const Real zmin = S.z0 - S.r_sink - S.L_noz;
        const Real zmax = S.z0 - S.r_sink;
        return (z >= zmin && z <= zmax);
    }
} // namespace ends

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
        S.rho_jet = pin->GetReal("problem", "rho_jet");
        S.A_noz = M_PI * S.r_noz * S.r_noz;                      // area per lobe
        S.Pdot_target = S.rho_jet * S.v_jet * S.v_jet * S.A_noz; // per-lobe target
    }
    // fixed momentum
    S.lock_pdot = 1;
    S.v_cap = 3.0;
    S.exclude_nozzle_from_accretion = 1; // avoid the sink accidentally accrete the gas we just injected as a jet
    S.M0 = -1.0;                         // set it later, bc it's impossible to find M0 = Msink +Mgas, without knowing Mgas(which is summed per MeshBlock inside SinkJetSource, requires integrating rho*V over every cell in every MeshBlock (and every MPI rank))

    // run my SinkJetSource function each timestep
    EnrollUserExplicitSourceFunction(SinkJetSource);
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

// Start of the real physics part !!!!!!!!!!!!!!!!!!!!!!!!!! is called every timestep

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
                    in_noz = in_nozzle_plus(xc, yc, zc) || in_nozzle_minus(xc, yc, zc); // both nozzle
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
    if (S.jet_switch)
    {
        // First: find the total nozzle volume of both lobes
        Real V_nozzle_total = 0.0;
        for (int k = pmb->ks; k <= pmb->ke; ++k)
        {
            const Real zc = x3v(k);
            for (int j = pmb->js; j <= pmb->je; ++j)
            {
                const Real yc = x2v(j);
                for (int i = pmb->is; i <= pmb->ie; ++i)
                {
                    const Real xc = x1v(i);
                    if (in_nozzle_plus(xc, yc, zc) || in_nozzle_minus(xc, yc, zc))
                    {
                        V_nozzle_total += pmb->pcoord->GetCellVolume(k, j, i);
                    }
                }
            }
        }

        if (V_nozzle_total > 0.0)
        {
            // total mass injection rate = rho_jet * v_jet * A_face * 2 (two lobes)
            Real rhoj = std::max(S.rho_jet, 1e-12);
            Real v_eff;
            if (S.lock_pdot)
            {
                // choose v_eff so that momentum flux matches Pdot_target
                v_eff = std::sqrt(S.Pdot_target / (rhoj * S.A_noz));
                if (S.v_cap > 0.0)
                {
                    v_eff = std::min(v_eff, S.v_cap); // cap the speed
                }
            }
            else
            {
                v_eff = S.v_jet;
            }
            const Real Mdot_lobe = rhoj * v_eff * S.A_noz; // one lobe
            const Real dM_total = (2.0 * Mdot_lobe) * dt;  // both lobes
            // find the mass need to subtracted from the sink
            const Real dM_used = std::min(dM_total, std::max(S.Msink, 0.0));
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
                        const bool plus = in_nozzle_plus(xc, yc, zc);              // in the upper lobe
                        const bool minus = (!plus) && in_nozzle_minus(xc, yc, zc); // in the downer lobe
                        if (!(plus || minus))
                            continue; // do nothing if not in nozzle

                        // find volume of the grid
                        const Real Vcell = pmb->pcoord->GetCellVolume(k, j, i);
                        const Real dM = dM_used * (Vcell / V_nozzle_total); // find the inject jet mass distributed per cell
                        const Real dRho = dM / Vcell;                       // compute density increment

                        // add to conserved density
                        cons(IDN, k, j, i) += dRho;

                        // calculate the injected momentum
                        const Real dP = dM * v_eff * (plus ? +1.0 : -1.0); //+z lobe momentum is positive, -z lobe momentum is negative
                        cons(IM3, k, j, i) += dP / Vcell;
                    }
                }
            }
            dMjet_local += dM_used;
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
}
