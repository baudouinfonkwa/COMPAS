// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
#ifndef RelaxationPressureTemperature_Compressible6Eq_IE_NPhase_H_
#define RelaxationPressureTemperature_Compressible6Eq_IE_NPhase_H_
#include <limits>
#include <cmath>
#if PHYSICS == SIXEQS_IE_NPHASE

enum class SixEqIE_PTStatus { Success, Nonfinite, NegativePhase, InvalidEOS,
    InvalidVolume, InvalidEnergy, NoConvergence };

// Instantaneous mechanical/thermal equilibrium at fixed partial masses,
// mixture momentum and conservative total energy. EOS 0--3 (constant A/B/C/D).
// alpha_k=b_k m_k + T w_k/(p+pi_k), w_k=m_k cv_k/A_k.
// V=1-sum b*m, Cp=sum m cv (1+1/A), pi=B/(A+1).
// With x=p+min(pi)>0, h=(H-sum mD)/V-min(pi)>0, solve
// Cp/S(x)-x=h, S=sum w/(x+pi-min(pi)). This function is strictly increasing.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
SixEqIE_PTStatus SixEqIE_TryPressureTemperatureCell(amrex::Array<amrex::Real,NSTATE>& out,
    amrex::Array<amrex::Real,NSTATE> const& state,Parm const& p,
    amrex::Real& pressure_out,amrex::Real& temperature_out)
{
    using amrex::Real;
    out=state;
    Real rho=0.,capacity=0.,sum_w=0.,free_volume=1.,offset=0.;
    Real pi_min=std::numeric_limits<Real>::max();
    amrex::Array<Real,NPHASE> w{},pi{};
    for(int n=0;n<NSTATE;++n)if(!std::isfinite(state[n]))
        return SixEqIE_PTStatus::Nonfinite;
    for(int l=0;l<NPHASE;++l){
        Real m=state[INDEX_Mass1+l];
        if(m<0. || state[l]<-64.*std::numeric_limits<Real>::epsilon())
            return SixEqIE_PTStatus::NegativePhase;
        if(!(p.A[l]>0.) || !(p.cv[l]>0.) || p.C[l]<0.)
            return SixEqIE_PTStatus::InvalidEOS;
        rho+=m;capacity+=m*p.cv[l];w[l]=m*p.cv[l]/p.A[l];sum_w+=w[l];
        free_volume-=p.C[l]*m;offset+=p.D[l]*m;
        pi[l]=p.B[l]/(1.+p.A[l]);if(m>0.)pi_min=std::min(pi_min,pi[l]);
    }
    if(!(rho>0.) || !(capacity>0.) || !(free_volume>0.) || !(sum_w>0.))
        return SixEqIE_PTStatus::InvalidVolume;
    Real kinetic=0.;for(int d=0;d<AMREX_SPACEDIM;++d){Real mom=state[INDEX_MomentumX+d];kinetic+=.5*mom*mom/rho;}
    const Real H=state[INDEX_Total_Energy]-kinetic;
    const Real h=(H-offset)/free_volume-pi_min,cp=capacity+sum_w;
    if(!(h>0.) || !std::isfinite(h))
        return SixEqIE_PTStatus::InvalidEnergy;
    Real lo=0.,hi=h*sum_w/capacity,x=hi,S=0.;
    bool converged=false;
    for(int it=0;it<192;++it){
        S=0.;Real S2=0.;
        for(int l=0;l<NPHASE;++l)if(w[l]>0.){Real y=x+(pi[l]-pi_min),v=w[l]/y;S+=v;S2+=v/y;}
        Real value=cp/S-x-h;
        if(std::abs(value)<=64.*std::numeric_limits<Real>::epsilon()*std::max(h,cp/S)){
            converged=true;break;
        }
        if(value>0.)hi=x;else lo=x;
        Real derivative=cp*S2/(S*S)-1.;
        Real next=x-value/derivative;
        if(!(next>lo && next<hi) || !std::isfinite(next))next=.5*(lo+hi);
        x=next;
    }
    if(!converged)return SixEqIE_PTStatus::NoConvergence;
    Real T=free_volume/S,pressure=x-pi_min;
    for(int l=0;l<NPHASE;++l){
        Real m=state[INDEX_Mass1+l];
        out[l]=m>0. ? p.C[l]*m+T*w[l]/(x+(pi[l]-pi_min)) : 0.;
        if(m>0. && !(out[l]-p.C[l]*m>0.))return SixEqIE_PTStatus::InvalidVolume;
        // Every pressure/thermal-energy conversion uses the canonical EOS.
        EquationOfState_Phase(out[INDEX_Energy1+l],out[l]*pressure,out[l],m,p,l,true,true);
    }
    pressure_out=pressure;temperature_out=T;
    return SixEqIE_PTStatus::Success;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void SixEqIE_RelaxPressureTemperatureCell(amrex::Array<amrex::Real,NSTATE>& out,
    amrex::Array<amrex::Real,NSTATE> const& state,Parm const& p)
{
    amrex::Real pressure,T;
    const auto status=SixEqIE_TryPressureTemperatureCell(out,state,p,pressure,T);
    switch(status){
    case SixEqIE_PTStatus::Success:return;
    case SixEqIE_PTStatus::Nonfinite:amrex::Abort("SixEq IE P-T relaxation: nonfinite state");break;
    case SixEqIE_PTStatus::NegativePhase:amrex::Abort("SixEq IE P-T relaxation: negative phase mass or fraction");break;
    case SixEqIE_PTStatus::InvalidEOS:amrex::Abort("SixEq IE P-T relaxation: invalid EOS/heat capacity");break;
    case SixEqIE_PTStatus::InvalidVolume:amrex::Abort("SixEq IE P-T relaxation: invalid density or free volume");break;
    case SixEqIE_PTStatus::InvalidEnergy:amrex::Abort("SixEq IE P-T relaxation: energy below equilibrium temperature bound");break;
    case SixEqIE_PTStatus::NoConvergence:amrex::Abort("SixEq IE P-T relaxation: equilibrium solve did not converge");break;
    }
}

AMREX_GPU_HOST inline
void SixEqIE_RelaxPressureTemperature(amrex::MultiFab& state,Parm const* p)
{
    for(amrex::MFIter mfi(state,amrex::TilingIfNotGPU());mfi.isValid();++mfi){
        auto a=state.array(mfi);auto box=mfi.tilebox();
        amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE(int i,int j,int k){
            amrex::Array<amrex::Real,NSTATE> in{},out{};
            for(int n=0;n<NSTATE;++n)in[n]=a(i,j,k,n);
            SixEqIE_RelaxPressureTemperatureCell(out,in,*p);
            for(int n=0;n<NSTATE;++n)a(i,j,k,n)=out[n];
        });
    }
}
#endif
#endif
