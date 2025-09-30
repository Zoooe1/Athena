#include <cmath>
#include <fstream>
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declaration avoid compile error
void SinkJetSource(MeshBlock *pmb, const Real time, const Real dt,
                   const AthenaArray<Real> &prim,
                   const AthenaArray<Real> &bcc,
                   const AthenaArray<Real> &prim_scalar,
                   AthenaArray<Real> &cons,
                   AthenaArray<Real> &cons_scalar);
// Params
namespace{
struct SinkParams{
    Real x0, y0, z0; //position of the sink particle
    Real M0; //initial total mass
    Real r_sink; //radius
    Real rho_thr; //density threshold to accrete
    Real rho_floor; //the density leave behind after accretion instead of making cells totally empty
    Real Msink; //mass of sink
    Real gsoft; // gravity softening length
    Real G; //will set to 1 for now
    
    
    int jet_switch; //0 is OFF, 1 is ON
    Real r_noz; //radius of nozzle
    Real L_noz; //length of nozzle
    Real v_jet; //jet speed along z axis
    Real rho_jet; //density used to set mass injection rate
    int  exclude_nozzle_from_accretion; // for debug, default on
    std::FILE* fp = nullptr; // CSV log file handle
}S; //S is a single storage box for all the sink+jet settings, like an object
    


// Check if a given cell is within the nozzle on top of the sink particle sphere
inline bool in_nozzle_plus(const Real x, const Real y, const Real z) {
    // +z lobe: axial window [z0 + r_sink, z0 + r_sink + L_noz]
    const Real r_perp = std::sqrt((x - S.x0)*(x - S.x0) + (y - S.y0)*(y - S.y0));
    if (r_perp > S.r_noz) return false;
    const Real zmin = S.z0 + S.r_sink;
    const Real zmax = S.z0 + S.r_sink + S.L_noz;
    return (z >= zmin && z <= zmax);
  }
// Check if a given cell is within the nozzle at the bottom of the sink particle sphere
  inline bool in_nozzle_minus(const Real x, const Real y, const Real z) {
    // -z lobe: axial window [z0 - r_sink - L_noz, z0 - r_sink]
    const Real r_perp = std::sqrt((x - S.x0)*(x - S.x0) + (y - S.y0)*(y - S.y0));
    if (r_perp > S.r_noz) return false;
    const Real zmin = S.z0 - S.r_sink - S.L_noz;
    const Real zmax = S.z0 - S.r_sink;
    return (z >= zmin && z <= zmax);
  }
} // namespace ends



//Function: read num from input file athinput.sinkjet, store them inside global S, carry it around for the whole simulation
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // sink
  S.x0        = pin->GetReal("problem","sink_x0");
  S.y0        = pin->GetReal("problem","sink_y0");
  S.z0        = pin->GetReal("problem","sink_z0");
  S.r_sink    = pin->GetReal("problem","r_sink");
  S.rho_thr   = pin->GetReal("problem","rho_thr");
  S.rho_floor = pin->GetReal("problem","rho_floor");
  S.Msink     = pin->GetReal("problem","Msink0");
  S.gsoft     = pin->GetReal("problem","g_soft");
  S.G         = pin->GetReal("problem","Gconst");

    
  // jets
  S.jet_switch    = pin->GetInteger("problem","jet_switch");
  if (S.jet_switch) {
    S.r_noz   = pin->GetReal("problem","r_noz");
    S.L_noz   = pin->GetReal("problem","L_noz");
    S.v_jet   = pin->GetReal("problem","v_jet");
    S.rho_jet = pin->GetReal("problem","rho_jet");
  }
  S.exclude_nozzle_from_accretion = 1; // avoid the sink accidentally accrete the gas we just injected as a jet
  S.M0 = -1.0;
  S.fp = std::fopen("sinkjet_history.csv","w");
  if (S.fp) std::fprintf(S.fp, "# time,Msink,Mgas,SFE_M0,SFE_inst\n"); //store the mass in csv

  // run my SinkJetSource function each timestep
  EnrollUserExplicitSourceFunction(SinkJetSource);
}

// How to initialize gas on the grid
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // // // // // // // // // // // // // // // // // // // //
  const Real rho0 = pin->GetReal("problem","rho0");
  for (int k=ks; k<=ke; ++k)
  for (int j=js; j<=je; ++j)
  for (int i=is; i<=ie; ++i) {
    phydro->u(IDN,k,j,i) = rho0; // density
    phydro->u(IM1,k,j,i) = 0.0;  // x-p :set momentum to 0 means that velocity is also set to 0, bc v = p/density
    phydro->u(IM2,k,j,i) = 0.0;  // y-p
    phydro->u(IM3,k,j,i) = 0.0;  // z-p
  }
  //This loop goes over every cell in the grid to initialize rho to 0
  // carve out the sink particle, remove initial gas inside it and set every cells to rho_floor
  for (int k=ks; k<=ke; ++k) {
    const Real zc = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      const Real yc = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        const Real xc = pcoord->x1v(i);
        const Real dx = xc - S.x0, dy = yc - S.y0, dz = zc - S.z0; //compute cell's distance to the sink center
        const Real r = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (r < S.r_sink) { //check if inside the sink sphere
          phydro->u(IDN,k,j,i) = S.rho_floor;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
        }
      }
    }
  }
}




//Start of the real physics part !!!!!!!!!!!!!!!!!!!!!!!!!! is called every timestep

void SinkJetSource(MeshBlock *pmb, const Real time, const Real dt,
                   const AthenaArray<Real> &prim,
                   const AthenaArray<Real> &bcc,
                   const AthenaArray<Real> &prim_scalar,
                   AthenaArray<Real> &cons,
                   AthenaArray<Real> &cons_scalar) {
    (void)bcc;          // not used (no MHD) Three useless parameters
    (void)prim_scalar;  // not used (NSCALARS=0)
    (void)cons_scalar;  // not used (NSCALARS=0)
    auto &x1v = pmb->pcoord->x1v; //make it easier to extract x,y,z coordinates of each cell
    auto &x2v = pmb->pcoord->x2v;
    auto &x3v = pmb->pcoord->x3v;
    // Accretion: remove gas inside r_sink above threshold
    //-------------------------------------------------------
    // loop over every cells in the grid excluding the jet nozzle cells so that we don't immediately reaccrete the gas we just injected
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
        const Real zc = x3v(k);
        for (int j=pmb->js; j<=pmb->je; ++j) {
            const Real yc = x2v(j);
            for (int i=pmb->is; i<=pmb->ie; ++i) {
                const Real xc = x1v(i);
                
                // nozzle mask (bipolar cylinders)
                bool in_noz = false; //we have to check if the cell is inside the sink particle
                if (S.exclude_nozzle_from_accretion && S.jet_switch) {
                    in_noz = in_nozzle_plus(xc,yc,zc) || in_nozzle_minus(xc,yc,zc);
                }
                if (in_noz) continue; // don't accrete nozzle cells
                
                const Real dx = xc - S.x0, dy = yc - S.y0, dz = zc - S.z0;
                const Real r  = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (r >= S.r_sink) continue; //outside sink particle
                
                const Real rho = cons(IDN,k,j,i);
                if (rho <= S.rho_thr) continue; //below density threshold
                
                const Real d_rho = rho - S.rho_floor; // Compute how much to remove
                if (d_rho <= 0.0) continue;
                
                // Msink += removed_gas_mass
                const Real Vcell = pmb->pcoord->GetCellVolume(k,j,i); //volume
                S.Msink += d_rho * Vcell;
                
                // remove corresponding momentum of the removed gas
                cons(IDN,k,j,i) -= d_rho;
                cons(IM1,k,j,i) -= d_rho * prim(IVX,k,j,i);
                cons(IM2,k,j,i) -= d_rho * prim(IVY,k,j,i);
                cons(IM3,k,j,i) -= d_rho * prim(IVZ,k,j,i);
                // (isothermal: no energy variable here)
            }
        }
    }
    
    
    // Gravity: the sink particle pulls on the surrounding gas, so give each gas cell a small momentum toward the sink ----------------------------------------------------
    // a = - GM r_hat / (r^2 + eps^2)^(3/2) ;;; update momentum = rho * a * dt
    
    for (int k=pmb->ks; k<=pmb->ke; ++k) { //loop over each cell
        const Real zc = x3v(k);
        for (int j=pmb->js; j<=pmb->je; ++j) {
            const Real yc = x2v(j);
            for (int i=pmb->is; i<=pmb->ie; ++i) {
                const Real xc = x1v(i);
                const Real dx = xc - S.x0, dy = yc - S.y0, dz = zc - S.z0; //displacement from cell to sink particle
                const Real r2_raw = dx*dx + dy*dy + dz*dz;
                if (r2_raw < S.r_sink*S.r_sink) continue;
                const Real r2 = r2_raw+ S.gsoft*S.gsoft;
                //gravitational softening to prevent gravity blow up near the center gsoft = r^2 + eps^2
                const Real r3 = r2*sqrt(r2);
                // find acc  direction towards the sink
                const Real gm_r3 = (S.G * S.Msink) / r3;  //compute scalar prefactor
                const Real ax = -gm_r3 * dx; //ax
                const Real ay = -gm_r3 * dy; //ay
                const Real az = -gm_r3 * dz; //az
                //convert acc to momentum change
                const Real rho = cons(IDN,k,j,i);
                cons(IM1,k,j,i) += rho * ax * dt;//Delta(\rho*v) = \rho * \a * dt
                cons(IM2,k,j,i) += rho * ay * dt;
                cons(IM3,k,j,i) += rho * az * dt;//densor gas gets bigger momentum change towards the sink particle
            }
        }
    }
    
    // Jets: inject mass & momentum in cells from two short cylinders(the nozzle) ---------------------------------------------------------
    // We distribute the mass rate uniformly across all nozzle cells by volume
    if (S.jet_switch) {
        // First: find the total nozzle volume of both lobes
        Real V_nozzle_total = 0.0;
        for (int k=pmb->ks; k<=pmb->ke; ++k) {
            const Real zc = x3v(k);
            for (int j=pmb->js; j<=pmb->je; ++j) {
                const Real yc = x2v(j);
                for (int i=pmb->is; i<=pmb->ie; ++i) {
                    const Real xc = x1v(i);
                    if (in_nozzle_plus(xc,yc,zc) || in_nozzle_minus(xc,yc,zc)) {
                        V_nozzle_total += pmb->pcoord->GetCellVolume(k,j,i);
                    }
                }
            }
        }
        
        if (V_nozzle_total > 0.0) {
            // total mass injection rate = rho_jet * v_jet * A_face * 2 (two lobes)
            const Real A_face = M_PI * S.r_noz * S.r_noz;
            const Real Mdot_total = S.rho_jet * S.v_jet * A_face * 2.0; //compute mass flow
            const Real dM_total   = Mdot_total * dt; // mass to add
            
            // Second: Distribute mass and momentum to each nozzle cell porportional to their volume; assign +z momentum in +lobe and -z in -lobe
            for (int k=pmb->ks; k<=pmb->ke; ++k) {
                const Real zc = x3v(k);
                for (int j=pmb->js; j<=pmb->je; ++j) {
                    const Real yc = x2v(j);
                    for (int i=pmb->is; i<=pmb->ie; ++i) {
                        const Real xc = x1v(i);
                        const bool plus  = in_nozzle_plus(xc,yc,zc); //in the upper lobe?
                        const bool minus = (!plus) && in_nozzle_minus(xc,yc,zc); //in the downer lobe?
                        if (!(plus || minus)) continue; //not in nozzle?
                        //find volume of the grid
                        const Real Vcell = pmb->pcoord->GetCellVolume(k,j,i);
                        const Real dM = dM_total * (Vcell / V_nozzle_total);//inject this cell with some jet mass
                        const Real dRho = dM / Vcell; //compute density increment that is going to be added later
                        
                        // add to conserved density
                        cons(IDN,k,j,i) += dRho;
                        
                        // give the injected mass a push along the +/- z-axis
                        const Real dPz = dM * S.v_jet * (plus ? +1.0 : -1.0); //+z lobe momentum is positive, -z lobe momentum is negative
                        cons(IM3,k,j,i) += dPz;
                    }
                }
            }
        }
    }
    // SFE logging --------------------------------
    Real Mgas = 0.0;
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
            for (int i=pmb->is; i<=pmb->ie; ++i) {
                const Real V = pmb->pcoord->GetCellVolume(k,j,i);
                Mgas += cons(IDN,k,j,i) * V;
            }
        }
    }
    static Real last_t_written = -1.0;
    if (time > last_t_written) {
        if (S.M0 < 0.0) S.M0 = S.Msink + Mgas; //after this M0 stays constant
        
        // Two common SFE definitions:
        const Real sfe_M0   = S.Msink / (S.M0 > 0 ? S.M0 : 1.0);                   // What fraction of the original reservoir has become stars by       time t
        const Real sfe_inst = S.Msink / std::max(S.Msink + Mgas, 1e-60);           // SFE_intantaneous = Msink / (Msink + gas at the point now)
        
        std::ofstream f("sinkjet_history.csv", std::ios::app);
        f.setf(std::ios::scientific); f.precision(16);
        f << time << "," << S.Msink << "," << Mgas << ","
        << sfe_M0 << "," << sfe_inst << "\n";
        last_t_written = time;
    }
}
