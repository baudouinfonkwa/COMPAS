// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef IC_H_
#define IC_H_

#include <AMReX_ParmParse.H>
#include <Parm.H>
#include "PhysicsFunctions.H"
#include <Tools.H>

// =======================================
// Viscous and conductive fluxes (DIFFUSION = true): five problems in one case, chosen at run
// time by prob.problem (Parm.H). Each problem is in its own header, in a namespace of its
// name, with the code of each model under #if (PHYSICS == ...). The functions below are the
// ones COMPAS calls; they call those of the problem of the run.
// =======================================

#include "Couette2Layer.H"
#include "HeatConduction2Mat.H"
#include "ViscousShockTube.H"
#include "AdvectionViscousDrop.H"
#include "AdvectionViscousMaterials.H"


// =======================================
// Initial condition
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void initial_condition(amrex::Array<amrex::Real, NSTATE>& Conservative,
                       bool& sharp,
                       amrex::Real x,
                       amrex::Real y,
                       amrex::Real z,
                       AMREX_D_DECL(amrex::Real dx,
                                    amrex::Real dy,
                                    amrex::Real dz),
                       AMREX_D_DECL(amrex::Real dx_FinestLevel,
                                    amrex::Real dy_FinestLevel,
                                    amrex::Real dz_FinestLevel),
                       Parm const& parm)
{
    switch (parm.problem)
    {
#ifdef DIFFUSION_HAS_COUETTE2LAYER
    case Couette2Layer_ID:
        Couette2Layer::initial_condition(Conservative, sharp, x, y, z,
                                         AMREX_D_DECL(dx, dy, dz),
                                         AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                         parm);
        break;
#endif
#ifdef DIFFUSION_HAS_HEATCONDUCTION2MAT
    case HeatConduction2Mat_ID:
        HeatConduction2Mat::initial_condition(Conservative, sharp, x, y, z,
                                              AMREX_D_DECL(dx, dy, dz),
                                              AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                              parm);
        break;
#endif
#ifdef DIFFUSION_HAS_VISCOUSSHOCKTUBE
    case ViscousShockTube_ID:
        ViscousShockTube::initial_condition(Conservative, sharp, x, y, z,
                                            AMREX_D_DECL(dx, dy, dz),
                                            AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                            parm);
        break;
#endif
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSDROP
    case AdvectionViscousDrop_ID:
        AdvectionViscousDrop::initial_condition(Conservative, sharp, x, y, z,
                                                AMREX_D_DECL(dx, dy, dz),
                                                AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                                parm);
        break;
#endif
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSMATERIALS
    case AdvectionViscousMaterials_ID:
        AdvectionViscousMaterials::initial_condition(Conservative, sharp, x, y, z,
                                                     AMREX_D_DECL(dx, dy, dz),
                                                     AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel),
                                                     parm);
        break;
#endif
    default:
        amrex::ignore_unused(Conservative, sharp, x, y, z, AMREX_D_DECL(dx, dy, dz),
                             AMREX_D_DECL(dx_FinestLevel, dy_FinestLevel, dz_FinestLevel));
        amrex::Abort("Diffusion: the problem is not available in this build");
    }
}


// =======================================
// Boundary condition: only Couette2Layer has a user condition (its moving no-slip walls)
// =======================================

#ifdef DIFFUSION_HAS_COUETTE2LAYER
#define USER_BOUNDARY_FUNC_ROBIN
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void BoundaryFunction_Robin(int i, int j, int k, int n,
                            Real& BC_Dirichlet,
                            Real& BC_Neumann,
                            Real& BC_Function,
                            Real const& xBC,
                            Real const& yBC,
                            Real const& zBC,
                            Real const& time,
                            Array4<Real const> const& q,
               AMREX_D_DECL(int ilo, int jlo, int klo),
               AMREX_D_DECL(int ihi, int jhi, int khi),
               AMREX_D_DECL(Real const& dx,
                            Real const& dy,
                            Real const& dz),
                            Parm const* parm)
{
    if (parm->problem == Couette2Layer_ID) {
        Couette2Layer::BoundaryFunction_Robin(i, j, k, n, BC_Dirichlet, BC_Neumann, BC_Function,
                                              xBC, yBC, zBC, time, q,
                                              AMREX_D_DECL(ilo, jlo, klo),
                                              AMREX_D_DECL(ihi, jhi, khi),
                                              AMREX_D_DECL(dx, dy, dz),
                                              parm);
        return;
    }
    amrex::ignore_unused(i, j, k, n, BC_Dirichlet, BC_Neumann, BC_Function, xBC, yBC, zBC, time, q,
                         AMREX_D_DECL(ilo, jlo, klo), AMREX_D_DECL(ihi, jhi, khi), AMREX_D_DECL(dx, dy, dz));
    amrex::Abort("Diffusion: only Couette2Layer has a user boundary condition (run.lo_bc, run.hi_bc)");
}
#endif


// =======================================
// In-situ diagnostics: the CSV file of the problem (couette.csv, conduction.csv, box.csv,
// advection.csv or materials.csv)
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
#if (PHYSICS == SIXEQS_IE_NPHASE)
    // Diagnose the unchanged physical diffusion operators before projections:
    // total-energy/kinetic-energy flux identity versus the actual heating source.
    using namespace UserOutput;
    const auto audit=Measure(o,{Sum,Sum,Sum,Max,Sum,Sum,Min,Min},1,
    [=] AMREX_GPU_DEVICE (Cell const& c,State const& U,Parm const& p)->Values<8>{
        amrex::Array<Real,NSTATE> state;
        Real rho=0.,H=0.,alpha_min=1.,mass_min=1.e300;
        for(int n=0;n<NSTATE;++n)state[n]=U(c.i,c.j,c.k,n);
        for(int l=0;l<NPHASE;++l){rho+=state[INDEX_Mass1+l];H+=state[INDEX_Energy1+l];
            alpha_min=std::min(alpha_min,state[l]);mass_min=std::min(mass_min,state[INDEX_Mass1+l]);}
        amrex::Array<Real,AMREX_SPACEDIM> velocity;
        Conservative2Velocity(AMREX_D_DECL(velocity[0],velocity[1],velocity[2]),state);
        Real K=0.;for(int d=0;d<AMREX_SPACEDIM;++d)K+=.5*rho*velocity[d]*velocity[d];
        amrex::Array<Real,NPHASE> Phi;
        FDM_Conservative2ViscousDissipation_K(c.i,c.j,c.k,Phi,U,AMREX_D_DECL(c.dx[0],c.dx[1],c.dx[2]),p);
        Real source=0.,compatible=0.;
        for(int l=0;l<NPHASE;++l)source+=PhaseDiffusion_DissipationFraction_D(state[l],state[INDEX_Mass1+l])*Phi[l];
        for(int d=0;d<AMREX_SPACEDIM;++d){
            amrex::Array<Real,NSTATE> flux[2];
            SixEqIEDiffusionFace face_data[2];
            for(int side=0;side<2;++side){
                amrex::Array<amrex::Array<Real,NSTATE>,FaceStencil::N> stencil;
                for(int q=0;q<FaceStencil::N;++q){int di,dj,dk;FaceStencil::OffsetOnFace(q,d,di,dj,dk);
                    di+=side*(d==0);dj+=side*(d==1);dk+=side*(d==2);
                    for(int n=0;n<NSTATE;++n)stencil[q][n]=U(c.i+di,c.j+dj,c.k+dk,n);}
                Real bound;FDM_Conservative2FluxDiffusion_D(flux[side],bound,stencil,AMREX_D_DECL(c.dx[0],c.dx[1],c.dx[2]),p,d,&face_data[side]);
            }
            Real jump=flux[1][INDEX_Total_Energy]-flux[0][INDEX_Total_Energy];
            jump-=face_data[1].heat-face_data[0].heat;
            for(int b=0;b<AMREX_SPACEDIM;++b)jump-=velocity[b]*(flux[1][INDEX_MomentumX+b]-flux[0][INDEX_MomentumX+b]);
            compatible+=jump/c.dx[d];
        }
        if(p.Physics_Parm.compatible_diffusion)source=compatible;
        Real E=state[INDEX_Total_Energy];
        return {E*c.dV,H*c.dV,K*c.dV,std::abs(H+K-E)/std::max(Real(1.),std::abs(E)),
                source*c.dV,compatible*c.dV,alpha_min,mass_min};
    });
    WriteRow(o,"diffusion_state.csv",{{"energy",audit[0]},{"thermal",audit[1]},{"kinetic",audit[2]},
        {"closure_max",audit[3]},{"heating_source",audit[4]},{"heating_compatible",audit[5]},
        {"alpha_min",audit[6]},{"mass_min",audit[7]}});
#endif
    switch (UserOutput::HostParm(o).problem)
    {
#ifdef DIFFUSION_HAS_COUETTE2LAYER
    case Couette2Layer_ID:             Couette2Layer::UserOutputFunction(o);             break;
#endif
#ifdef DIFFUSION_HAS_HEATCONDUCTION2MAT
    case HeatConduction2Mat_ID:        HeatConduction2Mat::UserOutputFunction(o);        break;
#endif
#ifdef DIFFUSION_HAS_VISCOUSSHOCKTUBE
    case ViscousShockTube_ID:          ViscousShockTube::UserOutputFunction(o);          break;
#endif
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSDROP
    case AdvectionViscousDrop_ID:      AdvectionViscousDrop::UserOutputFunction(o);      break;
#endif
#ifdef DIFFUSION_HAS_ADVECTIONVISCOUSMATERIALS
    case AdvectionViscousMaterials_ID: AdvectionViscousMaterials::UserOutputFunction(o); break;
#endif
    default: break;
    }
}

#endif  // IC_H_
