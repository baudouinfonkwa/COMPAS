// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PARM_H_
#define PARM_H_

// EOS pressure/energy inputs are SI; entropy reference data use the same units.
#define PARM_EOS_SI_UNITS

#include <PhysicsFunctions.H>
#include <cmath>
#include <string>
#include <vector>
#include <AMReX_REAL.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuMemory.H>

#include <AMReX_Array.H>
#include <AMReX_MLMG.H>

#include <Physics_Parameter_Compressible6Eq_IE_NPhase.H>
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
    PhaseField_Parameter PhaseField_Parm;

    //
    // Linear solver parameters (Grad_QPF)
    //
    LinearSystem_Parameter LinearSystem_Parm{.max_iter = 10'000, .max_fmg_iter = 10'000, .fixed_iter = 10'000, .verbose = 1};

    //
    // Equation of State (EOS)
    //
    int EOS = 1;

    // Convex EOS fit parameters
    amrex::Array<amrex::Real, NPHASE> A, B, C, D;
    amrex::Array<amrex::Real, NPHASE> gamma, pinf, gGamma, pref, eref, cv, mu, muB, kappa;


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
    amrex::Array<amrex::Real, NPHASE>     p_min;

    //
    // Problem parameters (read from prob.* in the inputs file). The solver works with
    // densities scaled by rho0, pressures by P0 and velocities by u0 = sqrt(P0/rho0)
    //
    amrex::Real bubble_diameter = 500.e-6; // length reference, metres
    amrex::Real P0   = 101325.0;   // pre-shock pressure, Pa
    amrex::Real rho0 = 1000.0;     // pre-shock water density and reference density, kg/m^3
    amrex::Real initial_temperature = 0.; // optional K; 0 retains prescribed upstream densities
    amrex::Real bubble_vapor_fraction = 0.05; // volume fraction within the three-phase bubble
    amrex::Real upstream_rho=1.;
    amrex::Array<amrex::Real,NPHASE> upstream_phase_rho{};
    amrex::Real rho_air = 1.0;     // pre-shock air density, kg/m^3
    amrex::Real Mach = 1.1;        // shock Mach number in water
    amrex::Real u0   = 0.0;        // reference velocity, set in DynamicInit

    //
    // Initialization routine
    //
    // Dimensionless pure-water shock state, recomputed after every EOS setup.
    amrex::Real shock_rho=1., shock_u=0., shock_p=1., shock_speed=0.;
    #define PARM_POST_EOS_INITIALIZE
    void PostEOSInitialize()
    {
        // Initialize reloads b and q in SI units on every call.
        if(EOS==3)for(int l=0;l<NPHASE;++l){C[l]*=rho0;D[l]*=rho0/P0;}
        static_assert(NPHASE==2 || NPHASE==3,"ShockBubble supports two phases or liquid/vapor/air");
        if(NPHASE==3 && !(initial_temperature>0.))
            amrex::Abort("Three-phase ShockBubble requires prob.initial_temperature>0");
        for(int l=0;l<NPHASE;++l){
            const amrex::Real pi=B[l]/(A[l]+1.);
            const amrex::Real rho=initial_temperature>0.
                ? 1./(C[l]+cv[l]/A[l]*initial_temperature/(1.+pi))
                : (l==0 ? 1. : rho_air/rho0);
            upstream_phase_rho[l]=rho;
            if(!(A[l]>0.) || !std::isfinite(A[l]) ||
               !(1.-C[l]*rho>0.) || !std::isfinite(C[l]) ||
               !(1.+B[l]/(A[l]+1.)>0.) || !std::isfinite(B[l]) || !std::isfinite(D[l]))
                amrex::Abort("ShockBubble: inadmissible upstream EOS state (A, p+p_infinity or 1-b*rho)");
        }
        amrex::Real shock[5];
        upstream_rho=upstream_phase_rho[0];
        StiffenedShockConditions(Mach,0.,1.,upstream_rho,A[0],B[0],C[0],D[0],shock);
        shock_rho=shock[0];shock_u=shock[1];shock_p=shock[2];shock_speed=shock[4];
        if(!(shock_rho>0.) || !(1.-C[0]*shock_rho>0.) ||
           !std::isfinite(shock_rho) || !std::isfinite(shock_u) ||
           !std::isfinite(shock_p) || !std::isfinite(shock_speed))
            amrex::Abort("ShockBubble: invalid post-shock state");
    }
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

   // Problem parameters, and the stiffened-gas pressures in units of P0
   {
      amrex::ParmParse ppp("prob");
      ppp.query("bubble_diameter",bubble_diameter);
      ppp.query("rho0", rho0);
      ppp.query("rho_air", rho_air);
      ppp.query("initial_temperature",initial_temperature);
      ppp.query("bubble_vapor_fraction",bubble_vapor_fraction);
      ppp.query("P0", P0);
      ppp.query("Mach", Mach);
   }

   if(!(P0>0.) || !(rho0>0.) || !(rho_air>0.) || !std::isfinite(P0) || !std::isfinite(rho0) || !std::isfinite(rho_air))
       amrex::Abort("prob.P0, prob.rho0 and prob.rho_air must be finite and positive");
   if(!(bubble_diameter>0.) || !std::isfinite(bubble_diameter))
       amrex::Abort("prob.bubble_diameter must be finite and positive (metres)");
   if(initial_temperature<0. || !std::isfinite(initial_temperature) ||
      bubble_vapor_fraction<0. || bubble_vapor_fraction>1. || !std::isfinite(bubble_vapor_fraction))
       amrex::Abort("ShockBubble requires initial_temperature>=0 and 0<=bubble_vapor_fraction<=1, finite");
   u0 = std::sqrt(P0/rho0);

   for (int ID_Phase = 0; ID_Phase < NPHASE; ID_Phase++){
      gamma[ID_Phase] = gammaV[ID_Phase]; 
      pinf[ID_Phase] = pinfV[ID_Phase]/P0; 
      gGamma[ID_Phase] = gGammaV[ID_Phase]; 
      pref[ID_Phase] = prefV[ID_Phase]/P0;
      eref[ID_Phase] = erefV[ID_Phase]*rho0/P0;
      cv[ID_Phase] = cvV[ID_Phase]*rho0/P0;
      mu[ID_Phase] = muV[ID_Phase]/(rho0*u0*bubble_diameter);
      muB[ID_Phase] = muBV[ID_Phase]/(rho0*u0*bubble_diameter);
      kappa[ID_Phase] = kappaV[ID_Phase]/(P0*u0*bubble_diameter);
   }

   if(!(Mach>=1.) || !std::isfinite(Mach))amrex::Abort("prob.Mach must be finite and >= 1");
   Initialize();
   if(!(PhaseField_Parm.Eta_Multiplier>=0.) || !std::isfinite(PhaseField_Parm.Eta_Multiplier))
       amrex::Abort("Initial interface width must be finite and nonnegative");


}

#endif
