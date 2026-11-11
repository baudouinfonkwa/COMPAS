// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PARM_H_
#define PARM_H_

#include <cmath>
#include <string>
#include <vector>
#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuMemory.H>

#include <AMReX_Array.H>
#include <AMReX_MLMG.H>

// One Parm.H for every build of this case (BUILD in the GNUmakefile). The values of
// FIVEEQS, SIXEQS, ... that PHYSICS is compared with are in Macros.H.
#include <Macros.H>
#if (PHYSICS == FIVEEQS)
#include <Physics_Parameter_Compressible5Eq.H>
#elif (PHYSICS == SIXEQS)
#include <Physics_Parameter_Compressible6Eq.H>
#elif (PHYSICS == FIVEEQS_NPHASE)
#include <Physics_Parameter_Compressible5Eq_NPhase.H>
#elif (PHYSICS == SIXEQS_IE_NPHASE)
#include <Physics_Parameter_Compressible6Eq_IE_NPhase.H>
#endif
#include <FiniteVolume_Parameter.H>
#include <PhaseField_Parameter.H>
#include <LinearSystem_Parameter.H>

struct Parm
{
    //
    // Physics control parameters
    //
    Physics_Parameter Physics_Parm;

    //
    // Finite Volume parameters
    //
    FiniteVolume_Parameter FiniteVolume_Parm;

    //
    // Phase-Field parameters
    //
#if (PHYSICS == SIXEQS)
    PhaseField_Parameter PhaseField_Parm{.ID_FluxPhase2Flux = 0};
#else
    PhaseField_Parameter PhaseField_Parm;
#endif

    //
    // Linear solver parameters (Grad_QPF)
    //
    LinearSystem_Parameter LinearSystem_Parm{.max_iter = 10'000, .max_fmg_iter = 10'000, .fixed_iter = 10'000, .verbose = 1};

    //
    // Equation of State (EOS)
    //
#if (PHYSICS == FIVEEQS)
    int EOS = 0;
#else
    int EOS = 1;
#endif

#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    // Ideal gas (EOS = 0)
    amrex::Real gamma_1;
    amrex::Real gamma_2;

    // Stiffened gas (EOS = 1)
    amrex::Real pinf_1, pinf_2;

    // Mie-Grüneisen (EOS = 2)
    amrex::Real gGamma_1 = 0.0, gGamma_2 = 0.0, pref_1 = 0.0, pref_2 = 0.0, eref_1 = 0.0, eref_2 = 0.0;

    // Thermal properties
    amrex::Real cv_1 = 1.0, cv_2 = 1.0;

    amrex::Real mu_1 = 0.0, mu_2 = 0.0, muB_1 = 0.0, muB_2 = 0.0, kappa_1 = 0.0, kappa_2 = 0.0;

    // Convex EOS fit parameters
    amrex::Real A1, A2, B1, B2, C1, C2, D1, D2;
#else
    // Convex EOS fit parameters
    amrex::Array<amrex::Real, NPHASE> A, B, C, D;
    amrex::Array<amrex::Real, NPHASE> gamma, pinf, gGamma, pref, eref, cv, mu, muB, kappa;
#endif

#if (PHYSICS == SIXEQS_IE_NPHASE)
    // Initial fractions: f_k + (1-sum f) alpha_raw,k. Not a transport clamp.
    amrex::Array<amrex::Real,NPHASE> alpha_floor = [] {
        amrex::Array<amrex::Real,NPHASE> values{};values.fill(1.e-6);return values;
    }();
    amrex::Real alpha_available=1.;

    void ReadInitialAlphaFloor()
    {
        amrex::ParmParse pp("prob");
        if(pp.contains("alpha_floor") && pp.countval("alpha_floor")!=NPHASE)
            amrex::Abort("prob.alpha_floor requires exactly NPHASE values");
        std::vector<amrex::Real> values(alpha_floor.begin(),alpha_floor.end());
        pp.queryarr("alpha_floor",values);
        amrex::Real sum=0.;
        for(int l=0;l<NPHASE;++l){
            if(!(values[l]>=0.) || !std::isfinite(values[l]))
                amrex::Abort("prob.alpha_floor requires finite nonnegative values for every phase");
            alpha_floor[l]=values[l];sum+=values[l];
        }
        if(!(sum<1.))amrex::Abort("prob.alpha_floor values must sum to less than one");
        alpha_available=1.-sum;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void ApplyInitialAlphaFloor(amrex::Array<amrex::Real,NPHASE>& alpha) const
    {
        for(int l=0;l<NPHASE;++l)alpha[l]=alpha_floor[l]+alpha_available*alpha[l];
    }
#endif

    //
    // Governing equation parameters
    //
    amrex::Real K = 0.0;
    amrex::Array<amrex::Real, NSTATE> UpperBound;
    amrex::Array<amrex::Real, NSTATE> LowerBound;
    amrex::Array<int,        NSTATE> ID_BoundType;
#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
    amrex::Array<amrex::Real, 2>     p_min;
#else
    amrex::Array<amrex::Real, NPHASE>     p_min;
#endif

    //
    // Problem parameters (read from prob.* in the inputs file)
    //
    int direction = 0; // Direction of the shock tube: 0 = x, 1 = y, 2 = z

    // Left (xi < x_diaphragm) and right states, with xi the coordinate along the tube.
    // The defaults are the Sod problem. The velocity is along the tube.
    amrex::Real rho_L = 1.0,   u_L = 0.0, P_L = 1.0;
    amrex::Real rho_R = 0.125, u_R = 0.0, P_R = 0.1;
    amrex::Real x_diaphragm = 0.5;

    // Material that fills each side (1 or 2)
    int mat_L = 1, mat_R = 1;

    //
    // Initialization routine
    //
    void Initialize();
    #define DYNAMIC_INIT
    void DynamicInit();
};

AMREX_FORCE_INLINE
void Parm::DynamicInit ()
{
#if (PHYSICS == SIXEQS_IE_NPHASE)
   ReadInitialAlphaFloor();
#endif
   amrex::ParmParse pp("prob.EOS");

#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
   pp.query("EOS",EOS);
   pp.query("gamma_1",gamma_1);
   pp.query("gamma_2",gamma_2);
   pp.query("pinf_1",pinf_1);
   pp.query("pinf_2",pinf_2);

   pp.query("gGamma_1",gGamma_1);
   pp.query("gGamma_2",gGamma_2);
   pp.query("pref_1",pref_1);
   pp.query("pref_2",pref_2);
   pp.query("eref_1",eref_1);
   pp.query("eref_2",eref_2);

   pp.query("cv_1",cv_1);
   pp.query("cv_2",cv_2);

   pp.query("mu_1",mu_1);
   pp.query("mu_2",mu_2);
   pp.query("muB_1",muB_1);
   pp.query("muB_2",muB_2);
   pp.query("kappa_1",kappa_1);
   pp.query("kappa_2",kappa_2);
#else
#if PHYSICS == SIXEQS_IE_NPHASE
   // Optional coefficients have physical defaults; every supplied array must
   // contain exactly one value per phase, including absent phases.
   std::vector<amrex::Real> gammaV(NPHASE,1.4), pinfV(NPHASE,0.), gGammaV(NPHASE,0.),
       prefV(NPHASE,0.), erefV(NPHASE,0.), cvV(NPHASE,1.),
       muV(NPHASE,0.), muBV(NPHASE,0.), kappaV(NPHASE,0.);
   auto phase_array = [&](char const* name, std::vector<amrex::Real>& values) {
       if(pp.contains(name) && pp.countval(name)!=NPHASE)
           amrex::Abort(std::string("prob.EOS.")+name+" requires exactly NPHASE values");
       pp.queryarr(name,values);
       if(values.size()!=NPHASE)
           amrex::Abort(std::string("prob.EOS.")+name+" requires exactly NPHASE values");
       for(auto value:values)if(!std::isfinite(value))
           amrex::Abort(std::string("prob.EOS.")+name+" contains a nonfinite value");
   };

   pp.query("EOS",EOS);
   phase_array("gamma",gammaV);
   phase_array("pinf",pinfV);

   phase_array("gGamma",gGammaV);
   phase_array("pref",prefV);
   phase_array("eref",erefV);

   phase_array("cv",cvV);

   phase_array("mu",muV);
   phase_array("muB",muBV);
   phase_array("kappa",kappaV);

#else
   std::vector<amrex::Real> gammaV(NPHASE), pinfV(NPHASE), gGammaV(NPHASE), prefV(NPHASE), erefV(NPHASE), cvV(NPHASE), muV(NPHASE), muBV(NPHASE), kappaV(NPHASE);

   pp.query("EOS",EOS);
   pp.queryarr("gamma",gammaV);
   pp.queryarr("pinf",pinfV);

   pp.queryarr("gGamma",gGammaV);
   pp.queryarr("pref",prefV);
   pp.queryarr("eref",erefV);

   pp.queryarr("cv",cvV);

   pp.queryarr("mu",muV);
   pp.queryarr("muB",muBV);
   pp.queryarr("kappa",kappaV);

#endif
   for (int ID_Phase = 0; ID_Phase < NPHASE; ID_Phase++){
      gamma[ID_Phase] = gammaV[ID_Phase];
      pinf[ID_Phase] = pinfV[ID_Phase];
      gGamma[ID_Phase] = gGammaV[ID_Phase];
      pref[ID_Phase] = prefV[ID_Phase];
      eref[ID_Phase] = erefV[ID_Phase];
      cv[ID_Phase] = cvV[ID_Phase];
      mu[ID_Phase] = muV[ID_Phase];
      muB[ID_Phase] = muBV[ID_Phase];
      kappa[ID_Phase] = kappaV[ID_Phase];
   }
#endif

   amrex::ParmParse pp_prob("prob");
   pp_prob.query("direction",direction);
   pp_prob.query("rho_L",rho_L);
   pp_prob.query("u_L",u_L);
   pp_prob.query("P_L",P_L);
   pp_prob.query("rho_R",rho_R);
   pp_prob.query("u_R",u_R);
   pp_prob.query("P_R",P_R);
   pp_prob.query("x_diaphragm",x_diaphragm);
   pp_prob.query("mat_L",mat_L);
   pp_prob.query("mat_R",mat_R);

#if PHYSICS == SIXEQS_IE_NPHASE
   if(mat_L<1 || mat_L>2 || mat_R<1 || mat_R>2)
       amrex::Abort("prob.mat_L and prob.mat_R must be 1 or 2");
   if(!(rho_L>0.) || !(rho_R>0.) || !std::isfinite(rho_L) || !std::isfinite(rho_R))
       amrex::Abort("Shock-tube densities must be finite and positive");
   if(direction<0 || direction>=AMREX_SPACEDIM)
       amrex::Abort("prob.direction must name a coordinate of this build");
#endif
   Initialize();

}

#endif
