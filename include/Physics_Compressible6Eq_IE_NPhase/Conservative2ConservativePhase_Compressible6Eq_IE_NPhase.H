// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef Conservative2ConservativePhase_Compressible6Eq_IE_NPhase_H_
#define Conservative2ConservativePhase_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void Conservative2ConservativePhase(amrex::Array<amrex::Real,(3+AMREX_SPACEDIM)*NPHASE>& ConservativePhase,
                                    amrex::Array<amrex::Real, NSTATE> const& Conservative,
                                    Parm const& EOS_Parameter)
{
    amrex::Array<amrex::Real, NPHASE> Alpha;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho;
    AMREX_D_TERM(amrex::Real rho_u;,
                 amrex::Real rho_v;,
                 amrex::Real rho_w );
    amrex::Real p, rhoE;
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_e;
    

    Conservative2State(Alpha,
                    Alpha_rho,
                    Alpha_rho_e,
       AMREX_D_DECL(rho_u, rho_v, rho_w),
                    rhoE,
                    Conservative);

    amrex::Real rho = 0.0;

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        rho += Alpha_rho[i_Phase];
    }

    AMREX_D_TERM(amrex::Real u = rho_u/rho;,
                amrex::Real v = rho_v/rho;,
                amrex::Real w = rho_w/rho;)
    
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_u;
#if (AMREX_SPACEDIM > 1)
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_v;
#endif
#if (AMREX_SPACEDIM > 2)
    amrex::Array<amrex::Real, NPHASE> Alpha_rho_w;
#endif
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_e_min;
    int isNegative = 0;
    amrex::Real rho_e = 0.0;

    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho_u[i_Phase] = Alpha_rho[i_Phase]*u;
#if (AMREX_SPACEDIM > 1)
        Alpha_rho_v[i_Phase] = Alpha_rho[i_Phase]*v;
#endif
#if (AMREX_SPACEDIM > 2)
        Alpha_rho_w[i_Phase] = Alpha_rho[i_Phase]*w;
#endif
        rho_e += Alpha_rho_e[i_Phase];
        bool PressureBase = true;//Pressure
        bool OutputInternalEnergy = true;//From Pressure to Alpha_rho_e
        
        amrex::Real p_mini = EOS_Parameter.p_min[i_Phase];
        amrex::Real Alphai_rhoi_ei_min;
        EquationOfState_Phase(Alphai_rhoi_ei_min,
                        Alpha[i_Phase]*p_mini,
                        Alpha[i_Phase],
                        Alpha_rho[i_Phase],
                        EOS_Parameter,
                        i_Phase,
                        PressureBase,
                        OutputInternalEnergy);
        Alpha_rho_e_min[i_Phase] = Alphai_rhoi_ei_min;
        
        if (Alpha_rho_e[i_Phase] < Alphai_rhoi_ei_min){
            isNegative = 1;
        }
    }
    
    amrex::Array<int,NPHASE> Phase_Indicator = {0};
    int projection_iterations=0;
    while (isNegative){
        if (++projection_iterations>NPHASE+1)
            amrex::Abort("SixEq IE phase conversion: energy bound redistribution did not converge");
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (Alpha_rho_e[i_Phase] < Alpha_rho_e_min[i_Phase]){
                Alpha_rho_e[i_Phase] = Alpha_rho_e_min[i_Phase];
                Phase_Indicator[i_Phase] = 1;
            }
        }
        
        amrex::Real Sum_Alpha_rho_e = 0.0;
        amrex::Real Sum_Alpha = 0.0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            Sum_Alpha_rho_e += Alpha_rho_e[i_Phase];
            if (Phase_Indicator[i_Phase] == 0){
                Sum_Alpha += Alpha[i_Phase];
            }
        }
        
        if (!(Sum_Alpha>0.))
            amrex::Abort("SixEq IE phase conversion: insufficient thermal energy for phase bounds");
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (Phase_Indicator[i_Phase] == 0){
                Alpha_rho_e[i_Phase] += Alpha[i_Phase]/Sum_Alpha*(rho_e - Sum_Alpha_rho_e);
            }
        }
        
        isNegative = 0;
        for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
            if (Alpha_rho_e[i_Phase] < Alpha_rho_e_min[i_Phase]){
                isNegative = 1;
            }
        }
    }
    
    amrex::Array<amrex::Real,NPHASE> Alpha_rho_E;
    for (int i_Phase = 0; i_Phase < NPHASE; i_Phase++){
        Alpha_rho_E[i_Phase] = (0.5)*Alpha_rho[i_Phase]*( AMREX_D_TERM(u*u,+v*v,+w*w) ) + Alpha_rho_e[i_Phase];
    }
    State2ConservativePhase(ConservativePhase,
                            Alpha_rho,
                            Alpha_rho_u,
#if (AMREX_SPACEDIM > 1)
                            Alpha_rho_v,
#endif
#if (AMREX_SPACEDIM > 2)
                            Alpha_rho_w,
#endif
                            Alpha_rho_E,
                            Alpha
                        );
    
}

#endif
#endif