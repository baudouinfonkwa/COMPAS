// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <Macros.H>
#include <LoopMacros.H>
#include <Parm.H>

#include <string>
#include <cmath>
#include <vector>

namespace {

// Names the value that failed a check. i_Phase >= 0 is the entry of an N-phase list.
std::string EOSValueName (std::string const& key, int i_Phase)
{
    return (i_Phase < 0) ? key : "value " + std::to_string(i_Phase+1) + " of " + key;
}

// An EOS key that the inputs file does not set is left at 0, which gamma > 1 rejects.
void CheckGamma (amrex::Real gamma, std::string const& key, int EOS, int i_Phase = -1)
{
    if (!(gamma > 1.0)) {
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " requires gamma > 1, but " + EOSValueName(key, i_Phase)
                     + " is " + std::to_string(gamma) + ". Set " + key + " in the inputs file.");
    }
}

// The temperature is (rho e - rho D - B_T)/(rho cv), so it needs cv > 0.
[[maybe_unused]] void CheckCv (amrex::Real cv, std::string const& key, std::string const& user, int i_Phase = -1)
{
    if (!(cv > 0.0)) {
        amrex::Abort(user + " uses the temperature, which requires cv > 0, but " + EOSValueName(key, i_Phase)
                     + " is " + std::to_string(cv) + ". Set " + key + " in the inputs file.");
    }
}


#if PHYSICS == SIXEQS_IE_NPHASE && DIFFUSION == true
void CheckTransport(amrex::Real value,std::string const& key,int phase)
{
    if(!std::isfinite(value) || value<0.)
        amrex::Abort("SIXEQS_IE_NPHASE diffusion requires finite nonnegative "+EOSValueName(key,phase));
}
#endif

// The NASG covolume b (stored in C) must be >= 0; b = 0 is the stiffened gas.
void CheckCovolume (amrex::Real b, std::string const& key, int i_Phase = -1)
{
    if (!(b >= 0.0)) {
        amrex::Abort("prob.EOS.EOS = 3 (NASG) requires a covolume b >= 0, but " + EOSValueName(key, i_Phase)
                     + " is " + std::to_string(b) + ".");
    }
}

#if (PHYSICS == FIVEEQS || PHYSICS == SIXEQS)
// NASG covolumes and heats of formation of the two phases, prob.EOS.b_1, b_2, q_1, q_2 (default 0).
// They are read only for EOS = 3, so that AMReX reports them as unused for the other EOS.
void ReadNASG (amrex::Real& b_1, amrex::Real& b_2, amrex::Real& q_1, amrex::Real& q_2)
{
    amrex::ParmParse pp("prob.EOS");
    b_1 = 0.0; b_2 = 0.0; q_1 = 0.0; q_2 = 0.0;
    pp.query("b_1", b_1);
    pp.query("b_2", b_2);
    pp.query("q_1", q_1);
    pp.query("q_2", q_2);
    CheckCovolume(b_1, "prob.EOS.b_1");
    CheckCovolume(b_2, "prob.EOS.b_2");
}
#endif

#if (PHYSICS == FIVEEQS_NPHASE || PHYSICS == SIXEQS_IE_NPHASE)
// NASG covolumes and heats of formation, the lists prob.EOS.b and prob.EOS.q (default 0).
// They are read only for EOS = 3, so that AMReX reports them as unused for the other EOS.
void ReadNASG (amrex::Array<amrex::Real,NPHASE>& b, amrex::Array<amrex::Real,NPHASE>& q)
{
    amrex::ParmParse pp("prob.EOS");
    std::vector<amrex::Real> b_list(NPHASE, 0.0), q_list(NPHASE, 0.0);
    pp.queryarr("b", b_list);
    pp.queryarr("q", q_list);
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        b[i_Phase] = b_list[i_Phase];
        q[i_Phase] = q_list[i_Phase];
        CheckCovolume(b[i_Phase], "prob.EOS.b", i_Phase);
    }
}

// Each per-phase EOS list needs exactly NPHASE values. A shorter list leaves the last phases unset.
void CheckPhaseLists ()
{
    amrex::ParmParse pp("prob.EOS");
    for (std::string const key : {"gamma", "pinf", "gGamma", "pref", "eref", "cv", "mu", "muB", "kappa", "b", "q"}){
        if (pp.contains(key.c_str()) && pp.countval(key.c_str()) != NPHASE){
            amrex::Abort("prob.EOS." + key + " has " + std::to_string(pp.countval(key.c_str()))
                         + " values, but NPHASE = " + std::to_string(NPHASE) + " needs one value per phase.");
        }
    }
}
#endif

}



#if (PHYSICS == FIVEEQS)

void Parm::Initialize ()
{
    //Read Physics parameters from the inputs file
    Physics_Parm.Initialize();

    //Read Finite Volume parameters from the inputs file
    FiniteVolume_Parm.Initialize();

    //Read Phase-Field parameters from the inputs file
    PhaseField_Parm.Initialize();

    //Read linear solver parameters for Grad_QPF from the inputs file
    LinearSystem_Parm.Initialize();

    if (EOS == 0 || EOS == 1 || EOS == 3){
        CheckGamma(gamma_1, "prob.EOS.gamma_1", EOS);
        CheckGamma(gamma_2, "prob.EOS.gamma_2", EOS);
    }
#if (DIFFUSION == true)
    CheckCv(cv_1, "prob.EOS.cv_1", "DIFFUSION = true");
    CheckCv(cv_2, "prob.EOS.cv_2", "DIFFUSION = true");
#endif

    switch(EOS){
    case 0:  // Ideal gas
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = 0.; B2 = 0.;
        C1 = 0.; C2 = 0.; 
        D1 = 0.; D2 = 0.;
        break;
    case 1: // SG
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = gamma_1*pinf_1/(gamma_1 - 1.); B2 = gamma_2*pinf_2/(gamma_2 - 1.);
        C1 = 0.; C2 = 0.; 
        D1 = 0.; D2 = 0.;
        break;
    case 2: // Mie-Gruneisen
        if (gGamma_1 <= 0.0 || gGamma_2 <= 0.0){
            amrex::Abort("prob.EOS.EOS = 2 (Mie-Gruneisen) requires prob.EOS.gGamma_1 > 0 and prob.EOS.gGamma_2 > 0.");
        }
        A1 = 1./gGamma_1; A2 = 1./gGamma_2;
        B1 = -pref_1/gGamma_1; B2 = -pref_2/gGamma_2;
        C1 = 0.; C2 = 0.; 
        D1 = eref_1; D2 = eref_2;
        break;
    case 3: // NASG: the covolume b in C, the heat of formation q in D
        ReadNASG(C1, C2, D1, D2);
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = gamma_1*pinf_1/(gamma_1 - 1.); B2 = gamma_2*pinf_2/(gamma_2 - 1.);
        break;
    default:
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " is not recognized. "
                     "Valid values are 0 (ideal gas), 1 (stiffened gas), 2 (Mie-Gruneisen), "
                     "3 (Noble-Abel stiffened gas).");
        break;
    }

    // The characteristic Jacobians of the WENO5 quadrature reconstruction assume a stiffened gas
    if (EOS == 3 && FiniteVolume_Parm.ID_Quad == 1){
        amrex::Abort("FiniteVolume.Quad = 1 cannot be combined with prob.EOS.EOS = 3 (NASG): the characteristic "
                     "decomposition of the quadrature reconstruction assumes a stiffened gas.");
    }

    if (Physics_Parm.source_term == 0){
        p_min[0] = -B1/A1;
        p_min[1] = -B2/A2;
    }
    else if (Physics_Parm.source_term == 1){
        p_min[0] = -B1/(A1 + 1.0);
        p_min[1] = -B2/(A2 + 1.0);
    }
    else {
        amrex::Abort("Physics.source_term must be 0 or 1.");
    }

    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    UpperBound = {inf, inf, AMREX_D_DECL( inf,  inf,  inf),   inf, 1.0};
    LowerBound = {0.0, 0.0, AMREX_D_DECL(-inf, -inf, -inf),  -inf, 0.0};
    ID_BoundType   = { -1,  -1, AMREX_D_DECL(   2,    2,    2),    2,   0};

}

#endif

#if (PHYSICS == FIVEEQS_NPHASE)

void Parm::Initialize ()
{
    //Read Physics parameters from the inputs file
    Physics_Parm.Initialize();

    //Read Finite Volume parameters from the inputs file
    FiniteVolume_Parm.Initialize();

    //Read Phase-Field parameters from the inputs file
    PhaseField_Parm.Initialize();

    //Read linear solver parameters for Grad_QPF from the inputs file
    LinearSystem_Parm.Initialize();

    CheckPhaseLists();
    if (EOS == 0 || EOS == 1 || EOS == 3){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            CheckGamma(gamma[i_Phase], "prob.EOS.gamma", EOS, i_Phase);
        }
    }
#if (DIFFUSION == true)
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        CheckCv(cv[i_Phase], "prob.EOS.cv", "DIFFUSION = true", i_Phase);
    }
#endif

    switch(EOS){
    case 0:  // Ideal gas
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = 0.0;
            C[i_Phase] = 0.0;
            D[i_Phase] = 0.0;
        }
        break;
    case 1: // SG
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = gamma[i_Phase]*pinf[i_Phase]/(gamma[i_Phase]-1.);
            C[i_Phase] = 0.0;
            D[i_Phase] = 0.0;
        }
        break;
    case 2: // Mie-Gruneisen
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (gGamma[i_Phase] <= 0.0){
                amrex::Abort("prob.EOS.EOS = 2 (Mie-Gruneisen) requires prob.EOS.gGamma > 0 for every phase.");
            }
            A[i_Phase] = 1./gGamma[i_Phase];
            B[i_Phase] = -pref[i_Phase]/gGamma[i_Phase];
            C[i_Phase] = 0.0;
            D[i_Phase] = eref[i_Phase];
        }
        break;
    case 3: // NASG: the covolume b in C, the heat of formation q in D
        ReadNASG(C, D);
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = gamma[i_Phase]*pinf[i_Phase]/(gamma[i_Phase]-1.);
        }
        break;
    default:
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " is not recognized. "
                     "Valid values are 0 (ideal gas), 1 (stiffened gas), 2 (Mie-Gruneisen), "
                     "3 (Noble-Abel stiffened gas).");
        break;
    }

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        if (Physics_Parm.source_term == 0){
            p_min[i_Phase] = -B[i_Phase]/A[i_Phase];
        }
        else if (Physics_Parm.source_term == 1){
            p_min[i_Phase] = -B[i_Phase]/(A[i_Phase] + 1.0);
        }
        else {
            amrex::Abort("Physics.source_term must be 0 or 1.");
        }
    }

    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    UpperBound = {D_LOOP(NPHASE,inf), AMREX_D_DECL( inf,  inf,  inf),   inf, D_LOOP(NPHASE,1.0)};
    LowerBound = {D_LOOP(NPHASE,0.0), AMREX_D_DECL(-inf, -inf, -inf),  -inf, D_LOOP(NPHASE,0.0)};
    ID_BoundType   = {D_LOOP(NPHASE,-1),  AMREX_D_DECL(   2,    2,    2),    2,   D_LOOP(NPHASE,0)};

}

#endif


#if (PHYSICS == SIXEQS)

void Parm::Initialize ()
{
    //Read Physics parameters from the inputs file
    Physics_Parm.Initialize();

    //Read Finite Volume parameters from the inputs file
    FiniteVolume_Parm.Initialize();

    //Read Phase-Field parameters from the inputs file
    PhaseField_Parm.Initialize();

    //Read linear solver parameters for Grad_QPF from the inputs file
    LinearSystem_Parm.Initialize();

    if (EOS == 0 || EOS == 1 || EOS == 3){
        CheckGamma(gamma_1, "prob.EOS.gamma_1", EOS);
        CheckGamma(gamma_2, "prob.EOS.gamma_2", EOS);
    }
    if (Physics_Parm.pressure_temperature_relaxation){
        CheckCv(cv_1, "prob.EOS.cv_1", "Physics.pressure_temperature_relaxation");
        CheckCv(cv_2, "prob.EOS.cv_2", "Physics.pressure_temperature_relaxation");
    }
#if (DIFFUSION == true)
    CheckCv(cv_1, "prob.EOS.cv_1", "DIFFUSION = true");
    CheckCv(cv_2, "prob.EOS.cv_2", "DIFFUSION = true");
    // Each phase conducts only through itself, and only the temperature relaxation passes heat
    // from one phase to the other (once: the cases call Initialize again from DynamicInit)
    static bool Warned = false;
    if (!Warned && !Physics_Parm.pressure_temperature_relaxation && (kappa_1 > 0.0 || kappa_2 > 0.0)){
        Warned = true;
        amrex::Print() << "Warning: SIXEQS with DIFFUSION = true and a nonzero conductivity, but without "
                       << "Physics.pressure_temperature_relaxation = 1: no heat passes from one phase to the other, "
                       << "so the heat flux across a material interface is missing, and more so as the grid is "
                       << "refined (see Known limitations in docs/user-guide/models.md).\n";
    }
#endif

    switch(EOS){
    case 0:  // Ideal gas
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = 0.; B2 = 0.;
        C1 = 0.; C2 = 0.; 
        D1 = 0.; D2 = 0.;
        break;
    case 1: // SG
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = gamma_1*pinf_1/(gamma_1 - 1.); B2 = gamma_2*pinf_2/(gamma_2 - 1.);
        C1 = 0.; C2 = 0.; 
        D1 = 0.; D2 = 0.;
        break;
     case 2: // Mie-Gruneisen
        if (gGamma_1 <= 0.0 || gGamma_2 <= 0.0){
            amrex::Abort("prob.EOS.EOS = 2 (Mie-Gruneisen) requires prob.EOS.gGamma_1 > 0 and prob.EOS.gGamma_2 > 0.");
        }
        A1 = 1./gGamma_1; A2 = 1./gGamma_2;
        B1 = -pref_1/gGamma_1; B2 = -pref_2/gGamma_2;
        C1 = 0.; C2 = 0.; 
        D1 = eref_1; D2 = eref_2;
        break;
     case 3: // NASG: the covolume b in C, the heat of formation q in D
        ReadNASG(C1, C2, D1, D2);
        A1 = 1./(gamma_1 - 1.); A2 = 1./(gamma_2 - 1.);
        B1 = gamma_1*pinf_1/(gamma_1 - 1.); B2 = gamma_2*pinf_2/(gamma_2 - 1.);
        break;
     default:
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " is not recognized. "
                     "Valid values are 0 (ideal gas), 1 (stiffened gas), 2 (Mie-Gruneisen), "
                     "3 (Noble-Abel stiffened gas).");
        break;
    }


    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    p_min[0] = -B1/(A1 + 1.0);
    p_min[1] = -B2/(A2 + 1.0);
    UpperBound = {inf, inf, AMREX_D_DECL( inf,  inf,  inf),     inf,     inf,  1.0};
    LowerBound = {0.0, 0.0, AMREX_D_DECL(-inf, -inf, -inf),    -inf,    -inf,  0.0};
    ID_BoundType   = { -1,  -1, AMREX_D_DECL(   2,    2,    2),      2,      2,    0};

}

#endif



#if (PHYSICS == SIXEQS_IE_NPHASE)

void Parm::Initialize ()
{
    //Read Physics parameters from the inputs file
    Physics_Parm.Initialize();

    //Read Finite Volume parameters from the inputs file
    FiniteVolume_Parm.Initialize();

    //Read Phase-Field parameters from the inputs file
    PhaseField_Parm.Initialize();

    if (PhaseField_Parm.phase_field){
        if (PhaseField_Parm.ID_FluxPhase2Flux != 1)
            amrex::Abort("SixEq IE phase field requires PhaseField.ID_FluxPhase2Flux=1: the phasic mapping transports total rather than internal energies");
        if ((PhaseField_Parm.ID_Flux != 0 && PhaseField_Parm.ID_Flux != 1) ||
            (PhaseField_Parm.ID_Formulation != 0 && PhaseField_Parm.ID_Formulation != 1))
            amrex::Abort("SixEq IE phase-field flux and formulation selectors must be 0 or 1");
    }

    //Read linear solver parameters for Grad_QPF from the inputs file
    LinearSystem_Parm.Initialize();

    CheckPhaseLists();
    if(Physics_Parm.PressureTemperatureRelaxation()){
        amrex::ParmParse eos_inputs("prob.EOS");
        if(!eos_inputs.contains("cv") || eos_inputs.countval("cv")!=NPHASE)
            amrex::Abort("Physics.Relaxation=PressureTemperature requires explicit prob.EOS.cv for every phase");
    }
    static bool AnnouncedPT=false;
    if(Physics_Parm.PressureTemperatureRelaxation() && !AnnouncedPT){
        AnnouncedPT=true;
        amrex::Print()<<"SIXEQS_IE_NPHASE: instantaneous joint pressure-temperature relaxation enabled.\n";
    }
    static bool AnnouncedCompatibleDiffusion=false;
    if(Physics_Parm.compatible_diffusion && !AnnouncedCompatibleDiffusion){
        AnnouncedCompatibleDiffusion=true;
        amrex::Print()<<"SIXEQS_IE_NPHASE: compatible viscous diffusion and perfect thermal contact enabled.\n";
    }

    if (EOS == 0 || EOS == 1 || EOS == 3){
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            CheckGamma(gamma[i_Phase], "prob.EOS.gamma", EOS, i_Phase);
        }
    }
#if (DIFFUSION == true)
    bool Conducts = false;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        CheckCv(cv[i_Phase], "prob.EOS.cv", "DIFFUSION = true", i_Phase);
        Conducts = Conducts || (kappa[i_Phase] > 0.0);
    }
    // Each phase conducts only through itself, and the model has no temperature relaxation to
    // pass heat from one phase to another (once: the cases call Initialize again from DynamicInit)
    static bool Warned = false;
    if (!Warned && Conducts && !Physics_Parm.compatible_diffusion && !Physics_Parm.PressureTemperatureRelaxation()){
        Warned = true;
        amrex::Print() << "Warning: SIXEQS_IE_NPHASE with DIFFUSION = true and a nonzero conductivity: the model "
                       << "has no temperature relaxation, so no heat passes from one phase to another, and the heat "
                       << "flux across a material interface is missing, and more so as the grid is refined "
                       << "(see Known limitations in docs/user-guide/models.md).\n";
    }
#endif

    switch(EOS){
    case 0:  // Ideal gas
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = 0.0;
            C[i_Phase] = 0.0;
            D[i_Phase] = 0.0;
        }
        break;
    case 1: // SG
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = gamma[i_Phase]*pinf[i_Phase]/(gamma[i_Phase]-1.);
            C[i_Phase] = 0.0;
            D[i_Phase] = 0.0;
        }
        break;
    case 2: // Mie-Gruneisen
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (gGamma[i_Phase] <= 0.0){
                amrex::Abort("prob.EOS.EOS = 2 (Mie-Gruneisen) requires prob.EOS.gGamma > 0 for every phase.");
            }
            A[i_Phase] = 1./gGamma[i_Phase];
            B[i_Phase] = -pref[i_Phase]/gGamma[i_Phase];
            C[i_Phase] = 0.0;
            D[i_Phase] = eref[i_Phase];
        }
        break;
    case 3: // NASG: the covolume b in C, the heat of formation q in D
        ReadNASG(C, D);
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            A[i_Phase] = 1./(gamma[i_Phase]-1.);
            B[i_Phase] = gamma[i_Phase]*pinf[i_Phase]/(gamma[i_Phase]-1.);
        }
        break;
    default:
        amrex::Abort("prob.EOS.EOS = " + std::to_string(EOS) + " is not recognized. "
                     "Valid values are 0 (ideal gas), 1 (stiffened gas), 2 (Mie-Gruneisen), "
                     "3 (Noble-Abel stiffened gas).");
        break;
    }

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
       
        p_min[i_Phase] = -B[i_Phase]/(A[i_Phase] + 1.0);
    
    }

#ifdef PARM_POST_EOS_INITIALIZE
    // A testcase may convert raw EOS data to its solver units. This must run
    // after every Initialize call, including the framework's second call.
    PostEOSInitialize();
#endif

#ifdef PARM_EOS_SI_UNITS
    Physics_Parm.InitializeChemicalEOS(P0,P0/rho0);
#else
    Physics_Parm.InitializeChemicalEOS();
#endif
    if(Physics_Parm.chemical_relaxation){
        for(int l=0;l<NPHASE;++l){
            const auto pi=B[l]/(A[l]+1.);
            if(!(Physics_Parm.chemical_pref+pi>0.))
                amrex::Abort("ChemicalRelaxation: EOS p_ref+p_infinity must be positive for every phase");
        }
        static bool announced=false;
        if(!announced){announced=true;amrex::Print()<<"SIXEQS_IE_NPHASE: equilibrium chemical relaxation enabled (immiscible phases; nucleation permitted).\n";}
    }

    // Joint equilibrium uses caloric data even in inviscid builds. Validate
    // coefficients after testcase nondimensionalization, including repeats.
    if(Physics_Parm.PressureTemperatureRelaxation())for(int l=0;l<NPHASE;++l){
        if(!std::isfinite(cv[l]) || !(cv[l]>0.) ||
           !std::isfinite(A[l]) || !(A[l]>0.) || !std::isfinite(B[l]) ||
           !std::isfinite(C[l]) || C[l]<0. || !std::isfinite(D[l]))
            amrex::Abort("PressureTemperature relaxation requires finite A>0, cv>0, b>=0 and finite B,D for every phase");
    }

#if (DIFFUSION == true)
    // Check the final coefficients, including testcase unit conversions.
    for (int l=0; l<NPHASE; ++l){
        CheckTransport(mu[l], "prob.EOS.mu", l);
        CheckTransport(muB[l], "prob.EOS.muB", l);
        CheckTransport(kappa[l], "prob.EOS.kappa", l);
        if (!std::isfinite(cv[l]) || !(cv[l]>0.))
            amrex::Abort("SIXEQS_IE_NPHASE diffusion requires finite positive " +
                         EOSValueName("prob.EOS.cv", l));
    }
#endif

    amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    UpperBound = {D_LOOP(NPHASE,1.0), D_LOOP(NPHASE,inf), D_LOOP(NPHASE,inf), AMREX_D_DECL( inf,  inf,  inf),   inf };
    LowerBound = {D_LOOP(NPHASE,0.0), D_LOOP(NPHASE,0.0), D_LOOP(NPHASE,0.0), AMREX_D_DECL(-inf, -inf, -inf),  -inf};
    ID_BoundType   = {D_LOOP(NPHASE,0), D_LOOP(NPHASE,-1), D_LOOP(NPHASE,-1),  AMREX_D_DECL(   2,    2,    2),    2};

}

#endif