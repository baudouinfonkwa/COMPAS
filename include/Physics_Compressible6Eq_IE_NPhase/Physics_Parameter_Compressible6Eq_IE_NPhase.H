// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Physics_Parameter_Compressible6Eq_IE_NPhase_H_
#define Physics_Parameter_Compressible6Eq_IE_NPhase_H_

#include<AMReX_ParmParse.H>
#include<Macros.H>
#include<string>
#include<vector>
#include<cmath>
#if (PHYSICS==SIXEQS_IE_NPHASE)

enum class Riemann : int {
    LLF,
    HLL,
    HLLC
};

enum class HLLCImplementation : int {
    ConservativePart = 0,
    Consistent = 1
};

enum class Relaxation : int { Pressure = 0, PressureTemperature = 1 };

struct Physics_Parameter
{
    int pressure_relaxation = 0;
    int compatible_diffusion = 0; // compatible viscous diffusion and perfect thermal contact
    Relaxation relaxation = Relaxation::Pressure;
    int chemical_relaxation = 0;
    int chemical_max_sweeps = 128;
    amrex::Real chemical_tolerance = 1.e-10;
    amrex::Real chemical_Tref = 1., chemical_pref = 1.;
    amrex::GpuArray<int,NPHASE> chemical_group{};
    amrex::GpuArray<amrex::Real,NPHASE> chemical_sref{};

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool PressureTemperatureRelaxation() const
    { return pressure_relaxation && relaxation == Relaxation::PressureTemperature; }

    Riemann RiemannSolver = Riemann::HLLC;
    int ID_Riemann = 2;
    HLLCImplementation hllc_implementation = HLLCImplementation::ConservativePart;

    // Reference entropy is anchored at (p_ref,T_ref), in the same raw
    // units as the testcase cv/pressure inputs. GPU storage uses solver units.
    void InitializeChemicalEOS(amrex::Real pressure_unit=1.,amrex::Real energy_unit=1.)
    {
        if(!chemical_relaxation)return;
        amrex::ParmParse eos("prob.EOS");
        if(!eos.contains("sref") || eos.countval("sref")!=NPHASE ||
           !eos.contains("T_ref") || !eos.contains("p_ref"))
            amrex::Abort("ChemicalRelaxation requires explicit prob.EOS.sref (NPHASE entries), T_ref and p_ref");
        std::vector<amrex::Real> s;
        eos.getarr("sref",s);eos.get("T_ref",chemical_Tref);eos.get("p_ref",chemical_pref);
        if(!(chemical_Tref>0.) || !(chemical_pref>0.) || !std::isfinite(chemical_Tref) || !std::isfinite(chemical_pref))
            amrex::Abort("ChemicalRelaxation requires finite positive EOS T_ref and p_ref");
        chemical_pref/=pressure_unit;
        for(int l=0;l<NPHASE;++l){
            if(!std::isfinite(s[l]))amrex::Abort("ChemicalRelaxation requires finite EOS sref entries");
            chemical_sref[l]=s[l]/energy_unit;
        }
    }

    //Read parameters from the inputs file
    void Initialize ()
    {
        amrex::ParmParse pp("Physics");
        pp.query("pressure_relaxation",pressure_relaxation);
        std::string hllc = "ConservativePart";
        pp.query("HLLCImplementation",hllc);
        if (hllc == "ConservativePart" || hllc == "0")
            hllc_implementation = HLLCImplementation::ConservativePart;
        else if (hllc == "Consistent" || hllc == "1")
            hllc_implementation = HLLCImplementation::Consistent;
        else
            amrex::Abort("Physics.HLLCImplementation must be ConservativePart (0) or Consistent (1)");
        pp.query("CompatibleDiffusion",compatible_diffusion);
        if(compatible_diffusion!=0 && compatible_diffusion!=1)
            amrex::Abort("Physics.CompatibleDiffusion must be 0 or 1");
        if(compatible_diffusion && !pressure_relaxation)
            amrex::Abort("Physics.CompatibleDiffusion requires Physics.pressure_relaxation=1");
#if (DIFFUSION != true || NONCONSERVATIVE != true)
        if(compatible_diffusion)
            amrex::Abort("Physics.CompatibleDiffusion requires DIFFUSION=true and NONCONSERVATIVE=true");
#endif
        // Preserve decks from before the independent relaxation selector:
        // CompatibleDiffusion=1 implied joint P-T equilibrium when unspecified.
        std::string choice = compatible_diffusion ? "PressureTemperature" : "Pressure";
        pp.query("Relaxation",choice);
        if(choice=="Pressure" || choice=="0")relaxation=Relaxation::Pressure;
        else if(choice=="PressureTemperature" || choice=="1")relaxation=Relaxation::PressureTemperature;
        else amrex::Abort("Physics.Relaxation must be Pressure (0) or PressureTemperature (1)");
        if(relaxation==Relaxation::PressureTemperature && !pressure_relaxation)
            amrex::Abort("Physics.Relaxation=PressureTemperature requires Physics.pressure_relaxation=1");
#if (DIFFUSION != true)
        if(relaxation==Relaxation::PressureTemperature)
            amrex::Abort("Physics.Relaxation=PressureTemperature requires DIFFUSION=true");
#endif
        pp.query("ChemicalRelaxation",chemical_relaxation);
        if(chemical_relaxation!=0 && chemical_relaxation!=1)
            amrex::Abort("Physics.ChemicalRelaxation must be 0 or 1");
        for(int l=0;l<NPHASE;++l)chemical_group[l]=l;
        if(chemical_relaxation){
            if(!PressureTemperatureRelaxation())
                amrex::Abort("Physics.ChemicalRelaxation requires Physics.Relaxation=PressureTemperature and pressure_relaxation=1");
            amrex::ParmParse chemistry("Physics.ChemicalRelaxation");
            chemistry.query("Tolerance",chemical_tolerance);
            chemistry.query("MaxSweeps",chemical_max_sweeps);
            if(!(chemical_tolerance>0.) || chemical_tolerance>1.e-3 || !std::isfinite(chemical_tolerance) || chemical_max_sweeps<1)
                amrex::Abort("ChemicalRelaxation requires 0<Tolerance<=1e-3 and MaxSweeps>=1");
            std::vector<int> pairs;
            chemistry.queryarr("Pairs",pairs);
            if(pairs.empty() || pairs.size()%2)
                amrex::Abort("Physics.ChemicalRelaxation.Pairs requires one or more zero-based phase-index pairs");
            amrex::ParmParse eos("prob.EOS");
            std::vector<std::string> substance;
            eos.queryarr("Substance",substance);
            if(substance.size()!=NPHASE)
                amrex::Abort("ChemicalRelaxation requires prob.EOS.Substance for every phase");
            for(std::size_t n=0;n<pairs.size();n+=2){
                int i=pairs[n],j=pairs[n+1];
                if(i<0 || j<0 || i>=NPHASE || j>=NPHASE || i==j)
                    amrex::Abort("ChemicalRelaxation.Pairs contains an invalid phase index or self-pair");
                if(substance[i]!=substance[j])
                    amrex::Abort("ChemicalRelaxation.Pairs must connect phases of the same Substance; chemical reactions are not supported");
                for(std::size_t m=0;m<n;m+=2)
                    if((pairs[m]==i && pairs[m+1]==j) || (pairs[m]==j && pairs[m+1]==i))
                        amrex::Abort("ChemicalRelaxation.Pairs contains a duplicate pair");
                int a=chemical_group[i],b=chemical_group[j];
                int root=std::min(a,b);
                for(int l=0;l<NPHASE;++l)if(chemical_group[l]==a || chemical_group[l]==b)chemical_group[l]=root;
            }
        }
        amrex::Real legacy_trace = 0.;
        pp.query("pressure_relaxation_trace",legacy_trace);
        if (legacy_trace != 0.){
            amrex::Abort("SIXEQS_IE_NPHASE COMPAS pressure relaxation does not freeze trace phases; remove Physics.pressure_relaxation_trace");
        }
    };
};

#endif
#endif
