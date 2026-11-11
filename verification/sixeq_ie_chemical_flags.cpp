// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#include <Physics_Parameter_Compressible6Eq_IE_NPhase.H>
#include <iostream>
int main(int argc,char**argv){
    amrex::Initialize(argc,argv);
    {
        Physics_Parameter p;p.Initialize();
        amrex::Real pu=1.,eu=1.;amrex::ParmParse pp("test");pp.query("pressure_unit",pu);pp.query("energy_unit",eu);
        p.InitializeChemicalEOS(pu,eu);const auto entropy=p.chemical_sref;
        const auto pref=p.chemical_pref;p.InitializeChemicalEOS(pu,eu);
        bool changed=p.chemical_pref!=pref;
        for(int l=0;l<NPHASE;++l)changed=changed || p.chemical_sref[l]!=entropy[l];
        if(changed)amrex::Abort("Chemical reference units changed on repeated initialization");
        std::cout<<"CHEMICAL_FLAGS {\"groups\":[";
        for(int l=0;l<NPHASE;++l){if(l)std::cout<<",";std::cout<<p.chemical_group[l];}
        std::cout<<"],\"p_ref\":"<<p.chemical_pref<<",\"sref0\":"<<p.chemical_sref[0]<<"}\n";
    }
    amrex::Finalize();
}
