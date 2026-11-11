// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
#ifndef RelaxationPressure_Instantaneous_Compressible6Eq_IE_NPhase_H_
#define RelaxationPressure_Instantaneous_Compressible6Eq_IE_NPhase_H_
#include <Macros.H>
#include <cmath>
#include <cstdio>
#include <limits>
#if PHYSICS == SIXEQS_IE_NPHASE

struct SixEqIERelaxationResult {
    amrex::Real work_pressure=0.,reset_pressure=0.,volume_residual=0.;
    int iterations=0;
    amrex::Real initial_work_pressure=0.,initial_residual=0.;
    int initial_halvings=0;
};

// g(x)=g_inf+sum c_k/(x+shift_k), x=P_e+min(pi_k)>0.
// All c_k are nonnegative. Evaluate g' without squaring tiny denominators.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void SixEqIE_RelaxationResidual(amrex::Real x,amrex::Real g_inf,
    amrex::Array<amrex::Real,NPHASE> const& c,
    amrex::Array<amrex::Real,NPHASE> const& shift,
    amrex::Real& g,amrex::Real& derivative,amrex::Real& scale)
{
    g=g_inf;derivative=0.;scale=std::abs(g_inf);
    for(int k=0;k<NPHASE;++k)if(c[k]>0.){
        const amrex::Real y=x+shift[k],term=c[k]/y;
        g+=term;derivative-=term/y;scale+=term;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void PressureRelaxation_Cell(amrex::Array<amrex::Real,NSTATE>& out,
    amrex::Array<amrex::Real,NSTATE> const& state,Parm const& parm,
    int i=-1,int j=-1,int k=-1,SixEqIERelaxationResult* details=nullptr)
{
    using amrex::Real;
    if(parm.EOS<0 || parm.EOS>3)amrex::Abort("SixEq IE relaxation: unsupported EOS identifier");
    constexpr Real eps=std::numeric_limits<Real>::epsilon();
    constexpr Real tolerance=32.*eps;
    out=state;
    Real rho=0.,mass_scale=1.;
    for(int n=0;n<NSTATE;++n)if(!std::isfinite(state[n]))
        amrex::Abort("SixEq IE pressure relaxation: invalid incoming conservative state");
    for(int l=0;l<NPHASE;++l){rho+=state[NPHASE+l];mass_scale=std::max(mass_scale,std::abs(state[NPHASE+l]));}
    if(!(rho>0.) || !std::isfinite(rho))
        amrex::Abort("SixEq IE pressure relaxation: invalid incoming conservative state");
    amrex::Array<Real,NSTATE> primitive;
    // Every incoming H -> alpha*p conversion goes through the canonical EOS.
    Conservative2Primitive(primitive,state,parm);
    amrex::Array<Real,NPHASE> c{},d{},pi{},shift{};
    Real sum_alpha=0.,sum_d=0.,sum_c=0.,coefficient_scale=0.;
    Real pi_min=std::numeric_limits<Real>::max();
    for(int l=0;l<NPHASE;++l){
        const Real alpha=state[l],mass=state[NPHASE+l],A=parm.A[l],B=parm.B[l],b=parm.C[l];
        if(alpha < -64.*eps || mass < -64.*eps*mass_scale){
            printf("SixEq IE incoming invalid cell=(%d,%d,%d) phase=%d alpha=%.17g mass=%.17g\n",i,j,k,l,alpha,mass);
            amrex::Abort("SixEq IE relaxation: invalid incoming volume fraction or partial mass");
        }
        if(!(A>0.) || !(b>=0.) || !std::isfinite(A) || !std::isfinite(B)
           || !std::isfinite(b) || !std::isfinite(parm.D[l]))
            amrex::Abort("SixEq IE relaxation: invalid EOS coefficients");
        sum_alpha+=alpha;pi[l]=B/(A+1.);
        if(b>0. && alpha<=0. && mass>64.*eps*mass_scale)
            amrex::Abort("SixEq IE relaxation: zero-volume NASG phase has finite mass");
        // No positive trace cutoff. Nonpositive roundoff/absent volumes do not
        // participate in volume exchange, as in the existing weighted-state policy.
        if(alpha<=0.)continue;
        const Real xi=NASG_Xi(alpha,mass,b),phi=xi*alpha;
        const Real Pi=xi*primitive[2*NPHASE+l]; // Pi=phi*p, not an EOS conversion.
        const Real gamma=(A+1.)/A;
        c[l]=(Pi+phi*pi[l])/gamma;d[l]=phi/gamma;
        coefficient_scale+=(std::abs(Pi)+std::abs(phi*pi[l]))/gamma;
        if(!std::isfinite(c[l]) || !std::isfinite(d[l]))
            amrex::Abort("SixEq IE relaxation: nonfinite pressure-work coefficients");
        sum_d+=d[l];pi_min=std::min(pi_min,pi[l]);
    }
    // Transport/reflux can leave signed roundoff in nominally absent phases.
    // Round only negative coefficients at mixture floating-point scale; this
    // is not a volume-fraction cutoff and no positive coefficient is changed.
    if(!std::isfinite(coefficient_scale))
        amrex::Abort("SixEq IE relaxation: nonfinite pressure-work scale");
    const Real coefficient_error=64.*eps*coefficient_scale;
    for(int l=0;l<NPHASE;++l){
        if(c[l]<-coefficient_error){
            printf("SixEq IE relaxation invalid work state cell=(%d,%d,%d) phase=%d alpha=%.17g mass=%.17g alpha_p=%.17g c=%.17g tolerance=%.17g\n",
                i,j,k,l,state[l],state[NPHASE+l],primitive[2*NPHASE+l],c[l],coefficient_error);
            amrex::Abort("SixEq IE relaxation: incoming phase is below its EOS temperature bound");
        }
        if(c[l]<0.)c[l]=0.;
        sum_c+=c[l];
    }
    const Real g_inf=(sum_alpha-1.)-sum_d;
    if(!(g_inf<0.) || !(sum_c>0.) || !std::isfinite(sum_c) || !std::isfinite(g_inf))
        amrex::Abort("SixEq IE relaxation: no admissible pressure-volume root");
    Real g_zero=g_inf;bool pole=false;
    for(int l=0;l<NPHASE;++l)if(state[l]>0.){
        shift[l]=pi[l]-pi_min;
        if(c[l]>0.){
            if(shift[l]==0.)pole=true;
            else g_zero+=c[l]/shift[l];
        }
    }
    if(!pole && !(g_zero>0.))
        amrex::Abort("SixEq IE relaxation: volume root is outside the common EOS pressure domain");
    // Start at the previous mixture pressure, keeping the right-most branch
    // x=P_e+pi_min>0. If p0 is outside that branch, choose a positive gap
    // using the incoming pressure scale (no fixed dimensional pressure floor).
    Real p0=0.;
    for(int l=0;l<NPHASE;++l)p0+=primitive[2*NPHASE+l];
    if(!std::isfinite(p0))amrex::Abort("SixEq IE relaxation: nonfinite initial mixture pressure");
    Real x=p0+pi_min;
    if(!(x>0.) || !std::isfinite(x)){
        x=std::max(std::abs(p0),std::abs(pi_min));
        if(!(x>0.))x=sum_c;
    }
    if(!(x>0.) || !std::isfinite(x))
        amrex::Abort("SixEq IE relaxation: invalid initial pressure gap");
    Real g,derivative,scale;
    SixEqIE_RelaxationResidual(x,g_inf,c,shift,g,derivative,scale);
    bool converged=std::isfinite(g) && std::abs(g)<=tolerance*scale;
    int halvings=0;
    // g decreases with pressure: halve the distance to the EOS bound until
    // g>0, unless the initial/search point already meets the root tolerance.
    // 2200 halvings cover the exponent range of finite double-precision gaps.
    for(;!converged && g<0. && halvings<2200;++halvings){
        const Real next=.5*x;
        if(!(next>0.) || next==x)
            amrex::Abort("SixEq IE relaxation: no representable positive-residual initial guess");
        x=next;SixEqIE_RelaxationResidual(x,g_inf,c,shift,g,derivative,scale);
        converged=std::isfinite(g) && std::abs(g)<=tolerance*scale;
    }
    if(!converged && (!(g>0.) || !std::isfinite(g) || !std::isfinite(scale)))
        amrex::Abort("SixEq IE relaxation: could not find a positive-residual initial guess");
    const Real initial_pressure=x-pi_min,initial_residual=g;
    int iterations=0;
    // Plain Newton-Raphson. For this decreasing convex residual, a positive
    // start approaches the root from below in pressure; a sign crossing is
    // not required. There is no bracket or midpoint/bisection fallback.
    for(;!converged && iterations<160;++iterations){
        if(!(derivative<0.) || !std::isfinite(derivative))
            amrex::Abort("SixEq IE relaxation: invalid Newton derivative");
        const Real candidate=x-g/derivative;
        if(!(candidate>0.) || !std::isfinite(candidate) || candidate==x)
            amrex::Abort("SixEq IE relaxation: Newton step left the pressure domain or stagnated");
        x=candidate;SixEqIE_RelaxationResidual(x,g_inf,c,shift,g,derivative,scale);
        if(!std::isfinite(g) || !std::isfinite(scale))
            amrex::Abort("SixEq IE relaxation: nonfinite Newton residual");
        converged=std::abs(g)<=tolerance*scale;
    }
    if(!converged)amrex::Abort("SixEq IE relaxation: pressure-volume root did not converge");
    const Real Pe=x-pi_min;
    Real alpha_sum=0.;
    for(int l=0;l<NPHASE;++l){
        if(state[l]>0.)out[l]=state[l]-d[l]+(c[l]>0.?c[l]/(x+shift[l]):0.);
        alpha_sum+=out[l];
        if(!std::isfinite(out[l]) || out[l]<-64.*eps)
        amrex::Abort("SixEq IE relaxation: invalid relaxed volume fraction");
    }
    if(std::abs(alpha_sum-1.)>256.*eps)
    amrex::Abort("SixEq IE relaxation: relaxed volumes do not sum to one");

    Real kinetic=0.;
    for(int n=0;n<AMREX_SPACEDIM;++n)kinetic+=.5*state[3*NPHASE+n]*state[3*NPHASE+n]/rho;
    const Real thermal=state[3*NPHASE+AMREX_SPACEDIM]-kinetic;
    Real offset=0.,compliance=0.;
    for(int l=0;l<NPHASE;++l){
        Real intercept;
        // p=0 evaluates the affine EOS intercept; it is not a stored phase state.
        EquationOfState_Phase(intercept,0.,out[l],state[NPHASE+l],parm,l,true,true);
        offset+=intercept;
        compliance+=parm.A[l]*NASG_Xi(out[l],state[NPHASE+l],parm.C[l])*out[l];
    }
    if(!(compliance>0.) || !std::isfinite(compliance) || !std::isfinite(thermal) || !std::isfinite(offset))
        amrex::Abort("SixEq IE relaxation: invalid energy-reset compliance or thermal budget");
    const Real pressure=(thermal-offset)/compliance;
    // ALWAYS reset common pressure and phasic energies from conserved mixture E,
    // using the conservative energy budget.
    for(int l=0;l<NPHASE;++l){
        if(!std::isfinite(pressure) || (out[l]>0. && !(pressure+pi[l]>0.))){
            printf("SixEq IE relaxation invalid reset cell=(%d,%d,%d) Pe=%.17g P_reset=%.17g thermal=%.17g\n",i,j,k,Pe,pressure,thermal);
            amrex::Abort("SixEq IE pressure relaxation produced an inadmissible state");
        }
        const Real ap=out[l]*pressure;
        EquationOfState_Phase(out[2*NPHASE+l],ap,out[l],state[NPHASE+l],parm,l,true,true);
        if(!std::isfinite(out[2*NPHASE+l]))amrex::Abort("SixEq IE relaxation: nonfinite reset phase energy");
    }
    if(details){details->work_pressure=Pe;details->reset_pressure=pressure;details->volume_residual=alpha_sum-1.;details->iterations=iterations;
        details->initial_work_pressure=initial_pressure;details->initial_residual=initial_residual;
        details->initial_halvings=halvings;}
}

AMREX_GPU_HOST AMREX_FORCE_INLINE
void PressureRelaxationInstantaneous(amrex::MultiFab& mf,Parm const parm)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    for(MFIter mfi(mf,TilingIfNotGPU());mfi.isValid();++mfi){
        const auto box=mfi.tilebox();const auto state=mf.array(mfi);
        amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE(int i,int j,int k){
            amrex::Array<amrex::Real,NSTATE> input,output;
            for(int c=0;c<NSTATE;++c)input[c]=state(i,j,k,c);
            PressureRelaxation_Cell(output,input,parm,i,j,k);
            for(int l=0;l<NPHASE;++l){state(i,j,k,l)=output[l];state(i,j,k,2*NPHASE+l)=output[2*NPHASE+l];}
        });
    }
}
#endif
#endif
