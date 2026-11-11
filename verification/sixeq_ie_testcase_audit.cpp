// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Check saved testcase states using the production EOS and conversion kernels.
#include <COMPAS.H>
#include <AMReX_VisMF.H>
#include <AMReX_BoxIterator.H>
#include <iostream>
#include <iomanip>
#include <stdexcept>
int main(int argc,char** argv)
{
    amrex::Initialize(argc,argv);
    {
        Parm parm{};parm.DynamicInit();parm.Initialize();
        amrex::ParmParse pp("audit");bool parameters_only=false;pp.query("parameters_only",parameters_only);
        if(parameters_only){
            std::cout<<std::setprecision(17)<<"CASE_PARAMETERS {";
            for(int c=0;c<7;++c){
                if(c)std::cout<<",";
                std::cout<<"\""<<(c==0?"B":c==1?"C":c==2?"D":c==3?"cv":c==4?"mu":c==5?"muB":"kappa")<<"\":[";
                for(int l=0;l<NPHASE;++l){if(l)std::cout<<",";std::cout<<(c==0?parm.B[l]:c==1?parm.C[l]:c==2?parm.D[l]:c==3?parm.cv[l]:c==4?parm.mu[l]:c==5?parm.muB[l]:parm.kappa[l]);}
                std::cout<<"]";
            }
            std::cout<<"}\n";amrex::Finalize();return 0;
        }
        std::string checkpoint;pp.get("checkpoint",checkpoint);
        int levels=1;pp.query("levels",levels);
        bool initial=false;pp.query("check_initial_floors",initial);
        double expected_pressure=0.,expected_temperature=0.,expected_velocity=0.;
        const bool check_pressure=pp.query("expected_pressure",expected_pressure);
        const bool check_temperature=pp.query("expected_temperature",expected_temperature);
        const bool check_velocity=pp.query("expected_velocity_x",expected_velocity);
        double initial_pressure_error=0.,initial_temperature_error=0.,initial_velocity_error=0.;
        amrex::Vector<amrex::MultiFab> state(levels);
        for(int lev=0;lev<levels;++lev)
            amrex::VisMF::Read(state[lev],checkpoint+"/Level_"+std::to_string(lev)+"/");
        double energy_error=0.,volume_error=0.,pressure_error=0.,eos_error=0.;
        double max_temperature_error=0.,max_chemical_error=0.,total_energy_sum=0.;
        amrex::Array<double,NPHASE> phase_mass_sums{};
        double min_alpha=1.,min_mass=1.e300,min_xi=1.;long cells=0;
        for(int lev=0;lev<levels;++lev){
            amrex::BoxArray covered;
            if(lev+1<levels){covered=state[lev+1].boxArray();covered.coarsen(2);}
            for(amrex::MFIter mfi(state[lev]);mfi.isValid();++mfi){
                auto const a=state[lev].const_array(mfi);
                for(amrex::BoxIterator it(mfi.validbox());it.ok();++it){
                    auto const iv=it();if(lev+1<levels && covered.contains(iv))continue;
                    amrex::Array<amrex::Real,NSTATE> u{},q{};
                    for(int n=0;n<NSTATE;++n){u[n]=a(iv,n);if(!std::isfinite(u[n]))throw std::runtime_error("Nonfinite saved state");}
                    Conservative2Primitive(q,u,parm);
                    for(auto value:q)if(!std::isfinite(value))throw std::runtime_error("Nonfinite primitive/EOS state");
                    double rho=0.,sumH=0.,sumAlpha=0.,P=0.;
                    for(int l=0;l<NPHASE;++l){rho+=u[NPHASE+l];sumH+=u[2*NPHASE+l];sumAlpha+=u[l];P+=q[2*NPHASE+l];}
                    if(!(rho>0.))throw std::runtime_error("Nonpositive mixture density");
                    double common_temperature=0.;
                    if(parm.Physics_Parm.PressureTemperatureRelaxation()){
                        double reference=0.;bool have_reference=false;
                        for(int l=0;l<NPHASE;++l)if(u[NPHASE+l]>0.){
                            double mt;EquationOfState_Phase(mt,u[2*NPHASE+l],u[l],u[NPHASE+l],parm,l,false,false);
                            double T=mt/u[NPHASE+l];
                            if(have_reference)max_temperature_error=std::max(max_temperature_error,std::abs(T-reference)/std::max(std::abs(T),std::abs(reference)));
                            else{reference=T;common_temperature=T;have_reference=true;}
                        }
                    }
                    if(check_pressure)initial_pressure_error=std::max(initial_pressure_error,std::abs(P-expected_pressure));
                    if(check_temperature){
                        double T;Conservative2Temperature(T,u,parm);
                        initial_temperature_error=std::max(initial_temperature_error,std::abs(T-expected_temperature));
                    }
                    if(check_velocity)initial_velocity_error=std::max(initial_velocity_error,std::abs(u[INDEX_MomentumX]/rho-expected_velocity));
                    if(parm.Physics_Parm.chemical_relaxation){
                        double entropy;
                        max_chemical_error=std::max(max_chemical_error,SixEqIE_ChemicalResidual(u,P,common_temperature,parm,entropy));
                        for(int l=0;l<NPHASE;++l)if(u[INDEX_Mass1+l]==0. && (u[l]!=0. || u[INDEX_Energy1+l]!=0.))
                            throw std::runtime_error("Absent chemical phase retains volume or energy");
                    }
                    const double weight=std::ldexp(1.,-AMREX_SPACEDIM*lev);
                    for(int l=0;l<NPHASE;++l)phase_mass_sums[l]+=weight*u[INDEX_Mass1+l];
                    total_energy_sum+=weight*u[INDEX_Total_Energy];
                    double K=0.;for(int d=0;d<AMREX_SPACEDIM;++d)K+=.5*u[3*NPHASE+d]*u[3*NPHASE+d]/rho;
                    energy_error=std::max(energy_error,std::abs(sumH+K-u[INDEX_Total_Energy])/std::max(1.,std::abs(u[INDEX_Total_Energy])));
                    volume_error=std::max(volume_error,std::abs(sumAlpha-1.));
                    for(int l=0;l<NPHASE;++l){
                        if(initial && !parm.Physics_Parm.PressureTemperatureRelaxation() && u[l]<parm.alpha_floor[l]*(1.-1.e-12))
                            throw std::runtime_error("Initial phase fraction is below its requested floor");
                        min_alpha=std::min(min_alpha,u[l]);min_mass=std::min(min_mass,u[NPHASE+l]);
                        min_xi=std::min(min_xi,NASG_Xi(u[l],u[NPHASE+l],parm.C[l]));
                        pressure_error=std::max(pressure_error,std::abs(q[2*NPHASE+l]-u[l]*P)/std::max(1.,std::abs(P)));
                        double H;EquationOfState_Phase(H,q[2*NPHASE+l],u[l],u[NPHASE+l],parm,l,true,true);
                        if(!std::isfinite(H))throw std::runtime_error("Nonfinite reconstructed phase energy");
                        eos_error=std::max(eos_error,std::abs(H-u[2*NPHASE+l])/std::max(1.,std::abs(u[2*NPHASE+l])));
                    }
                    ++cells;
                }
            }
        }
        if(!cells || volume_error>1.e-9 || min_alpha < -1.e-11 || min_mass < -1.e-9 || eos_error>1.e-10)
            throw std::runtime_error("Saved state violates fraction/mass/EOS checks");
        // Initial cell averaging can retain unresolved kinetic energy and
        // pressure disequilibrium. Require reset closure only after evolution.
        if(!initial && parm.Physics_Parm.pressure_relaxation && (energy_error>1.e-9 || pressure_error>1.e-9))
            throw std::runtime_error("Relaxed state violates energy/common-pressure closure");
        if(parm.Physics_Parm.PressureTemperatureRelaxation() && max_temperature_error>1.e-9)
            throw std::runtime_error("P-T relaxed state has unequal temperatures");
        if(parm.Physics_Parm.chemical_relaxation && max_chemical_error>5.*parm.Physics_Parm.chemical_tolerance)
            throw std::runtime_error("Chemically relaxed state violates active-set Gibbs equilibrium");
        if((check_pressure && initial_pressure_error>1.e-9*std::max(1.,std::abs(expected_pressure))) ||
           (check_temperature && initial_temperature_error>1.e-9*std::max(1.,std::abs(expected_temperature))) ||
           (check_velocity && initial_velocity_error>1.e-9*std::max(1.,std::abs(expected_velocity))))
            throw std::runtime_error("Saved initial state is not uniform at the requested p,T,u");
        std::cout<<std::setprecision(17)<<"CASE_AUDIT {\"cells\":"<<cells<<",\"max_energy_error\":"<<energy_error
            <<",\"max_volume_error\":"<<volume_error<<",\"max_common_pressure_error\":"<<pressure_error
            <<",\"max_eos_roundtrip_error\":"<<eos_error<<",\"max_common_temperature_error\":"<<max_temperature_error<<",\"minimum_alpha\":"<<min_alpha
            <<",\"minimum_mass\":"<<min_mass<<",\"minimum_xi\":"<<min_xi
            <<",\"initial_pressure_error\":"<<initial_pressure_error<<",\"initial_temperature_error\":"<<initial_temperature_error<<",\"initial_velocity_error\":"<<initial_velocity_error
            <<",\"max_chemical_error\":"<<max_chemical_error<<",\"total_energy_sum\":"<<total_energy_sum<<",\"phase_mass_sums\":[";
        for(int l=0;l<NPHASE;++l){if(l)std::cout<<",";std::cout<<phase_mass_sums[l];}
        std::cout<<"]}\n";
    }
    amrex::Finalize();
}
