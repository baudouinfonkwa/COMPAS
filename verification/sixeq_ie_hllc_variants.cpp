// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Production-kernel checks for the two SIXEQS_IE_NPHASE HLLC formulations.
#include <COMPAS.H>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using State=amrex::Array<amrex::Real,NSTATE>;
static_assert(NPHASE==3 && AMREX_SPACEDIM==2,"Build this audit for three phases in 2D");
void near(double x,double y,double scale=0.){
    if(!std::isfinite(x) || !std::isfinite(y) || std::abs(x-y)>3.e-11*std::max({std::abs(x),std::abs(y),scale,1.e-250})){
        std::cerr<<std::setprecision(17)<<"actual="<<x<<" expected="<<y<<" scale="<<scale<<"\n";
        throw std::runtime_error("HLLC audit identity failed");
    }
}
Parm eos(int id){
    Parm p{};p.EOS=id;
    for(int k=0;k<NPHASE;++k){
        p.A[k]=1.5+.5*k;p.B[k]=id==0?0.:(id==2?-.15*(k+1):.5*(p.A[k]+1.));
        p.C[k]=id==3?.02*(k+1):0.;p.D[k]=id>=2?.4*(k-1):0.;p.cv[k]=1.+k;
    }
    return p;
}
State make_state(Parm const& p,bool right,double trace,double ux=.3,double uy=-.1,bool same_pressure=false){
    State v{},u;
    for(int k=0;k<NPHASE;++k){
        const double alpha=k==0?(right?.4:.3):(k==1?(right?.6:.7)-trace:trace);
        const double density=(right?.8:1.1)*(k+1);
        v[k]=alpha;v[NPHASE+k]=alpha*density;
        v[2*NPHASE+k]=alpha*(same_pressure?2.:(right?1.7:2.)+.1*k);
    }
    v[3*NPHASE]=ux;v[3*NPHASE+1]=uy;
    Primitive2Conservative(u,v,p);return u;
}
int main(int argc,char** argv){
    amrex::Initialize(argc,argv);
    std::string mode;amrex::ParmParse().query("audit_mode",mode);
    if(mode=="parse"){
        Physics_Parameter physics;physics.Initialize();
        amrex::Print()<<"PARSED_IMPLEMENTATION "<<static_cast<int>(physics.hllc_implementation)<<"\n";
        amrex::Finalize();return 0;
    }
    if(mode=="pole" || mode=="covolume" || mode=="bad_phase"){
        Parm p=eos(3);p.A[0]=1.;p.C[0]=mode=="covolume"?.4:0.;
        SixEqIE_HLLC_AlphaPressureStar(.2,.2,.4,mode=="bad_phase"?0.:.2,mode=="pole"?4.:3.,0,p);
        throw std::runtime_error("Invalid star was accepted");
    }
    int jumps=0,flux_checks=0,offset_checks=0,trace_checks=0;
    double max_work_residual=0.,max_alpha_flux_error=0.;
    for(int id=0;id<4;++id){
        Parm p=eos(id);
        for(double alpha:{.3,1.e-16,1.e-100})for(double r:{.8,1.,1.25})for(double lambda:{1.,r,.95*r}){
            const double mass=1.2*alpha,ap=2.*alpha,astar=lambda*alpha,mstar=r*mass;
            const double aps=SixEqIE_HLLC_AlphaPressureStar(alpha,mass,ap,astar,r,0,p);
            double H,Hs;EquationOfState_Phase(H,ap,alpha,mass,p,0,true,true);
            EquationOfState_Phase(Hs,aps,astar,mstar,p,0,true,true);
            const double residual=Hs/mstar-H/mass+.5*(aps/astar+ap/alpha)*(astar/mstar-alpha/mass);
            near(residual,0.,1.);max_work_residual=std::max(max_work_residual,std::abs(residual));
            if(lambda==r){near(aps,r*ap);near(Hs,r*H);}
            ++jumps;
        }
        near(SixEqIE_HLLC_AlphaPressureStar(0.,0.,0.,0.,1.2,0,p),0.);
        near(SixEqIE_HLLC_AlphaPressureStar(0.,0.,1.e-18,0.,1.2,0,p,1.),0.);
        for(auto implementation:{HLLCImplementation::ConservativePart,HLLCImplementation::Consistent}){
            p.Physics_Parm.hllc_implementation=implementation;
            for(double trace:{.1,1.e-16,0.})for(int dir=0;dir<AMREX_SPACEDIM;++dir){
                const State left=make_state(p,false,trace),right=make_state(p,true,trace,-.2,.2);
                State flux;amrex::Array<double,AMREX_SPACEDIM> velocity;double speed;
                Riemann_HLLC(flux,velocity,speed,left,right,p,dir);
                double alpha_flux=0.;for(int k=0;k<NPHASE;++k)alpha_flux+=flux[k];
                near(alpha_flux,velocity[dir],1.);max_alpha_flux_error=std::max(max_alpha_flux_error,std::abs(alpha_flux-velocity[dir]));++trace_checks;
                // All equal states, including supersonic states, recover physical flux.
                for(double velocity_scale:{1.,100.,-100.}){
                    State equal=make_state(p,false,trace,velocity_scale*.3,velocity_scale*.2),expected,actual;
                    Conservative2Flux(expected,equal,p,dir);Riemann_HLLC(actual,velocity,speed,equal,equal,p,dir);
                    for(int c=0;c<NSTATE;++c)near(actual[c],expected[c],1.);
                    ++flux_checks;
                }
                // Moving/stationary material contacts preserve pressure and velocity.
                for(double advect:{-.3,0.,.3}){
                    State l=make_state(p,false,trace,advect,advect,true),rr=make_state(p,true,trace,advect,advect,true),expected,actual;
                    Conservative2Flux(expected,advect>=0.?l:rr,p,dir);
                    Riemann_HLLC(actual,velocity,speed,l,rr,p,dir);
                    for(int c=0;c<NSTATE;++c)near(actual[c],expected[c],1.);
                    near(velocity[dir],advect,1.);++flux_checks;
                }
                // Independent reference-energy changes must only add D*mass flux.
                Riemann_HLLC(flux,velocity,speed,left,right,p,dir);
                Parm shifted=p;State l=left,rr=right,changed;
                for(int k=0;k<NPHASE;++k){
                    const double delta=(k==1?-1.:1.)*(k+1)*3.;shifted.D[k]+=delta;
                    l[2*NPHASE+k]+=delta*l[NPHASE+k];rr[2*NPHASE+k]+=delta*rr[NPHASE+k];
                    l[INDEX_Total_Energy]+=delta*l[NPHASE+k];rr[INDEX_Total_Energy]+=delta*rr[NPHASE+k];
                }
                amrex::Array<double,AMREX_SPACEDIM> changed_velocity;double changed_speed;
                Riemann_HLLC(changed,changed_velocity,changed_speed,l,rr,shifted,dir);
                double expected_energy=flux[INDEX_Total_Energy];
                for(int k=0;k<NPHASE;++k){
                    const double delta=(k==1?-1.:1.)*(k+1)*3.;
                    near(changed[k],flux[k],1.);near(changed[NPHASE+k],flux[NPHASE+k],1.);
                    near(changed[2*NPHASE+k],flux[2*NPHASE+k]+delta*flux[NPHASE+k],1.);
                    expected_energy+=delta*flux[NPHASE+k];
                }
                near(changed[INDEX_Total_Energy],expected_energy,1.);near(changed_speed,speed);
                for(int d=0;d<AMREX_SPACEDIM;++d){near(changed_velocity[d],velocity[d],1.);near(changed[3*NPHASE+d],flux[3*NPHASE+d],1.);}
                ++offset_checks;
            }
        }
    }
    amrex::Print()<<std::setprecision(17)<<"HLLC_RESULTS {\"jump_cases\":"<<jumps<<",\"flux_cases\":"<<flux_checks
        <<",\"reference_shift_cases\":"<<offset_checks<<",\"simplex_flux_cases\":"<<trace_checks
        <<",\"maximum_work_residual\":"<<max_work_residual<<",\"maximum_simplex_flux_error\":"<<max_alpha_flux_error<<"}\n";
    amrex::Finalize();
}
