// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#include <Physics_Parameter_Compressible6Eq_IE_NPhase.H>
#include <iostream>
int main(int argc,char**argv){amrex::Initialize(argc,argv);{Physics_Parameter p;p.Initialize();std::cout<<"DIFFUSION_FLAG "<<p.compatible_diffusion<<" RELAXATION "<<static_cast<int>(p.relaxation)<<"\n";}amrex::Finalize();}
