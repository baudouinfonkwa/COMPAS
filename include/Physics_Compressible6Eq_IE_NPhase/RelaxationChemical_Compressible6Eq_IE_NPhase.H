// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
#ifndef RelaxationChemical_Compressible6Eq_IE_NPhase_H_
#define RelaxationChemical_Compressible6Eq_IE_NPhase_H_
#include <RelaxationPressureTemperature_Compressible6Eq_IE_NPhase.H>
#include <cstdio>
#if PHYSICS == SIXEQS_IE_NPHASE

struct SixEqIE_ChemicalThermo {
    amrex::Real v=0.,e=0.,s=0.,g=0.,scale=1.;
};

// Pure-phase chemical potentials: gases are immiscible; no mixing entropy.
// sref is entropy at (p_ref,T_ref), not the unnormalised logarithmic constant.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool SixEqIE_ChemicalThermodynamics(SixEqIE_ChemicalThermo& t,
    amrex::Real pressure,amrex::Real T,int l,Parm const& p)
{
    using amrex::Real;
    const Real pi=p.B[l]/(p.A[l]+1.),R=p.cv[l]/p.A[l],cp=p.cv[l]+R;
    if(!(T>0.) || !(pressure+pi>0.))return false;
    t.v=p.C[l]+R*T/(pressure+pi);
    t.s=p.Physics_Parm.chemical_sref[l]
        +cp*std::log(T/p.Physics_Parm.chemical_Tref)
        -R*std::log((pressure+pi)/(p.Physics_Parm.chemical_pref+pi));
    const Real h=cp*T+p.C[l]*pressure+p.D[l];
    t.e=h-pressure*t.v;t.g=h-T*t.s;
    t.scale=std::max({Real(1.),std::abs(h),T*std::abs(t.s),p.cv[l]*T});
    return std::isfinite(t.v) && std::isfinite(t.e) && std::isfinite(t.s) && std::isfinite(t.g);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real SixEqIE_ChemicalViolation(amrex::Real f,amrex::Real scale,
    amrex::Real mi,amrex::Real mj)
{
    if(mi>0. && mj>0.)return std::abs(f)/scale;
    if(mi>0.)return std::max(f,amrex::Real(0.))/scale;
    if(mj>0.)return std::max(-f,amrex::Real(0.))/scale;
    return 0.;
}

// Entropy and active-set KKT residual. If an absent phase has no pure state
// at (p,T), its potential cannot certify the endpoint. Force an interior
// flash search instead; notably, ideal vapor can nucleate in stretched liquid.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real SixEqIE_ChemicalResidual(amrex::Array<amrex::Real,NSTATE> const& u,
    amrex::Real pressure,amrex::Real T,Parm const& p,amrex::Real& entropy)
{
    using amrex::Real;
    amrex::Array<SixEqIE_ChemicalThermo,NPHASE> thermo{};
    amrex::Array<bool,NPHASE> valid{};
    entropy=0.;Real residual=0.;
    for(int l=0;l<NPHASE;++l){
        valid[l]=SixEqIE_ChemicalThermodynamics(thermo[l],pressure,T,l,p);
        if(u[INDEX_Mass1+l]>0.){
            if(!valid[l])amrex::Abort("Chemical relaxation: inadmissible active phase thermodynamics");
            entropy+=u[INDEX_Mass1+l]*thermo[l].s;
        }
    }
    for(int i=0;i<NPHASE;++i)for(int j=i+1;j<NPHASE;++j)
        if(p.Physics_Parm.chemical_group[i]==p.Physics_Parm.chemical_group[j]){
            const Real mi=u[INDEX_Mass1+i],mj=u[INDEX_Mass1+j];
            if(mi==0. && mj==0.)continue;
            if(valid[i] && valid[j])
                residual=std::max(residual,SixEqIE_ChemicalViolation(thermo[i].g-thermo[j].g,
                    std::max(thermo[i].scale,thermo[j].scale),mi,mj));
            else residual=std::max(residual,Real(1.));
        }
    return residual;
}

struct SixEqIE_ChemicalPoint {
    amrex::Array<amrex::Real,NSTATE> state{};
    amrex::Real y=0.,pressure=0.,T=0.,entropy=0.,f=0.,scale=1.,derivative=0.;
    bool valid=false;
};

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
SixEqIE_ChemicalPoint SixEqIE_ChemicalTrial(amrex::Real y,int i,int j,
    amrex::Real total,amrex::Array<amrex::Real,NSTATE> const& base,Parm const& p)
{
    using amrex::Real;
    SixEqIE_ChemicalPoint point;point.y=y;
    auto trial=base;
    trial[INDEX_Mass1+j]=total*y;
    trial[INDEX_Mass1+i]=total-trial[INDEX_Mass1+j];
    if(SixEqIE_TryPressureTemperatureCell(point.state,trial,p,point.pressure,point.T)!=SixEqIE_PTStatus::Success)return point;
    SixEqIE_ChemicalThermo a,b;
    bool va=SixEqIE_ChemicalThermodynamics(a,point.pressure,point.T,i,p);
    bool vb=SixEqIE_ChemicalThermodynamics(b,point.pressure,point.T,j,p);
    if((y<1. && !va) || (y>0. && !vb))return point;
    // An undefined absent-phase potential forces exploration of the interior.
    point.f=va && vb ? a.g-b.g : (va ? std::numeric_limits<Real>::max() : -std::numeric_limits<Real>::max());
    point.scale=std::max(a.scale,b.scale);
    Real unused=SixEqIE_ChemicalResidual(point.state,point.pressure,point.T,p,point.entropy);
    amrex::ignore_unused(unused);
    if(va && vb){
        Real vp=0.,vt=0.,ep=0.,et=0.;
        for(int l=0;l<NPHASE;++l)if(point.state[INDEX_Mass1+l]>0.){
            const Real m=point.state[INDEX_Mass1+l],pi=p.B[l]/(p.A[l]+1.);
            const Real dvt=p.cv[l]/p.A[l]/(point.pressure+pi),dvp=-dvt*point.T/(point.pressure+pi);
            vp+=m*dvp;vt+=m*dvt;ep+=m*pi*dvp;et+=m*(p.cv[l]+pi*dvt);
        }
        const Real determinant=vp*et-vt*ep;
        if(determinant!=0.){
            const Real rhs_v=-total*(b.v-a.v),rhs_e=-total*(b.e-a.e);
            const Real dp=(rhs_v*et-vt*rhs_e)/determinant,dT=(vp*rhs_e-rhs_v*ep)/determinant;
            point.derivative=(a.v-b.v)*dp+(b.s-a.s)*dT;
        }
    }
    point.valid=true;return point;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void SixEqIE_ChemicalConsider(SixEqIE_ChemicalPoint const& candidate,
    SixEqIE_ChemicalPoint& best,int i,int j)
{
    using amrex::Real;
    if(!candidate.valid)return;
    const Real roundoff=256.*std::numeric_limits<Real>::epsilon()*std::max({Real(1.),std::abs(best.entropy),std::abs(candidate.entropy)});
    const Real old_error=SixEqIE_ChemicalViolation(best.f,best.scale,best.state[INDEX_Mass1+i],best.state[INDEX_Mass1+j]);
    const Real new_error=SixEqIE_ChemicalViolation(candidate.f,candidate.scale,candidate.state[INDEX_Mass1+i],candidate.state[INDEX_Mass1+j]);
    if(candidate.entropy>best.entropy+roundoff ||
       (candidate.entropy>=best.entropy-roundoff && new_error<old_error))best=candidate;
}

// Feasible interior compositions have positive free volume and positive
// H-offset-pi_min*V. Both are affine in y while both pair members are active.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool SixEqIE_ChemicalPositiveInterval(amrex::Real f0,amrex::Real f1,
    amrex::Real& lo,amrex::Real& hi)
{
    if(f0<=0. && f1<=0.)return false;
    if(f0<=0.)lo=std::max(lo,f0/(f0-f1));
    if(f1<=0.)hi=std::min(hi,f0/(f0-f1));
    return hi>lo;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
SixEqIE_ChemicalPoint SixEqIE_ChemicalOptimizePair(int i,int j,
    amrex::Array<amrex::Real,NSTATE> const& base,Parm const& p)
{
    using amrex::Real;
    const Real total=base[INDEX_Mass1+i]+base[INDEX_Mass1+j];
    auto best=SixEqIE_ChemicalTrial(base[INDEX_Mass1+j]/total,i,j,total,base,p);
    if(!best.valid)amrex::Abort("Chemical relaxation: invalid starting pair state");
    const Real tolerance=p.Physics_Parm.chemical_tolerance;
    if(SixEqIE_ChemicalViolation(best.f,best.scale,best.state[INDEX_Mass1+i],best.state[INDEX_Mass1+j])<=tolerance)return best;
    auto left=SixEqIE_ChemicalTrial(0.,i,j,total,base,p);
    auto right=SixEqIE_ChemicalTrial(1.,i,j,total,base,p);
    if(left.valid && left.f<=tolerance*left.scale)SixEqIE_ChemicalConsider(left,best,i,j);
    if(right.valid && right.f>=-tolerance*right.scale)SixEqIE_ChemicalConsider(right,best,i,j);
    if(SixEqIE_ChemicalViolation(best.f,best.scale,best.state[INDEX_Mass1+i],best.state[INDEX_Mass1+j])<=tolerance)return best;

    Real rho=0.,offset0=0.,V0=1.,pi_min=std::min(p.B[i]/(p.A[i]+1.),p.B[j]/(p.A[j]+1.));
    for(int l=0;l<NPHASE;++l){
        rho+=base[INDEX_Mass1+l];
        const Real m=l==i?total:l==j?0.:base[INDEX_Mass1+l];
        V0-=m*p.C[l];offset0+=m*p.D[l];
        if(m>0.)pi_min=std::min(pi_min,p.B[l]/(p.A[l]+1.));
    }
    Real H=base[INDEX_Total_Energy];
    for(int d=0;d<AMREX_SPACEDIM;++d){Real mom=base[INDEX_MomentumX+d];H-=.5*mom*mom/rho;}
    const Real V1=V0-total*(p.C[j]-p.C[i]),offset1=offset0+total*(p.D[j]-p.D[i]);
    Real lo=0.,hi=1.;
    if(!SixEqIE_ChemicalPositiveInterval(V0,V1,lo,hi) ||
       !SixEqIE_ChemicalPositiveInterval(H-offset0-pi_min*V0,H-offset1-pi_min*V1,lo,hi))return best;
    // Interior limits exclude only floating-point endpoints, never a mass floor.
    Real yl=lo==0. && left.valid ? 0. : std::nextafter(lo,hi);
    Real yr=hi==1. && right.valid ? 1. : std::nextafter(hi,lo);
    auto a=yl==0.?left:SixEqIE_ChemicalTrial(yl,i,j,total,base,p);
    auto b=yr==1.?right:SixEqIE_ChemicalTrial(yr,i,j,total,base,p);
    const Real middle=.5*(lo+hi);
    // Near a thermodynamic boundary, roundoff may make the limit inadmissible.
    for(int n=0;n<64 && !a.valid;++n){yl=.5*(yl+middle);a=SixEqIE_ChemicalTrial(yl,i,j,total,base,p);}
    for(int n=0;n<64 && !b.valid;++n){yr=.5*(yr+middle);b=SixEqIE_ChemicalTrial(yr,i,j,total,base,p);}
    if(!a.valid || !b.valid || !(a.y<b.y) || !(a.f>=0.) || !(b.f<=0.))return best;
    auto x=SixEqIE_ChemicalTrial(.5*(a.y+b.y),i,j,total,base,p);
    for(int iteration=0;iteration<128;++iteration){
        if(!x.valid)amrex::Abort("Chemical relaxation: inadmissible interior mass-transfer trial");
        if(std::abs(x.f)<=tolerance*x.scale){SixEqIE_ChemicalConsider(x,best,i,j);return best;}
        if(x.f>0.)a=x;else b=x;
        Real next=x.y-x.f/x.derivative;
        if(!(x.derivative<0.) || !std::isfinite(next) || !(next>a.y && next<b.y))next=.5*(a.y+b.y);
        if(next==a.y || next==b.y)break;
        x=SixEqIE_ChemicalTrial(next,i,j,total,base,p);
    }
    amrex::Abort("Chemical relaxation: bounded pair solve did not converge");
    return best;
}

// A chemical equilibrium flash conserves substance-group masses, not the
// incoming phase partition. Reflux can invalidate that partition after exact
// disappearance. Construct a feasible solver seed from the conserved groups;
// this is not an output clamp and does not change any group mass or energy.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
bool SixEqIE_ChemicalGroupSeed(amrex::Array<amrex::Real,NSTATE>& out,
    amrex::Array<amrex::Real,NSTATE> const& input,Parm const& p,
    amrex::Real& pressure,amrex::Real& T)
{
    using Real=amrex::Real;
    amrex::Array<Real,NPHASE> mass{},uniform{};
    amrex::Array<int,NPHASE> count{};
    for(int n=0;n<NSTATE;++n)if(!std::isfinite(input[n]))return false;
    Real rho=0.,pi_min=std::numeric_limits<Real>::max();
    for(int l=0;l<NPHASE;++l){int g=p.Physics_Parm.chemical_group[l];mass[g]+=input[INDEX_Mass1+l];++count[g];}
    for(int g=0;g<NPHASE;++g)if(mass[g]<0.)return false;else rho+=mass[g];
    if(!(rho>0.))return false;
    Real H=input[INDEX_Total_Energy],scale=std::max(Real(1.),std::abs(H));
    for(int d=0;d<AMREX_SPACEDIM;++d){Real mom=input[INDEX_MomentumX+d];H-=.5*mom*mom/rho;}
    for(int l=0;l<NPHASE;++l){
        int g=p.Physics_Parm.chemical_group[l];uniform[l]=mass[g]/count[g];
        if(mass[g]>0.)pi_min=std::min(pi_min,p.B[l]/(p.A[l]+1.));
        scale=std::max(scale,std::abs(p.D[l])*mass[g]);
    }
    scale=std::max(scale,std::abs(pi_min));
    auto trial=input;
    auto project=[&](amrex::Array<Real,NPHASE> const& m)->bool {
        for(int l=0;l<NPHASE;++l){trial[INDEX_Mass1+l]=m[l];trial[l]=m[l]>0.?1./NPHASE:0.;}
        return SixEqIE_TryPressureTemperatureCell(out,trial,p,pressure,T)==SixEqIE_PTStatus::Success;
    };
    // A fully active composition has two affine admissibility constraints:
    // f=1-sum b*m>0, g=(H-sum q*m-pi_min*f)/scale>0.
    auto constraints=[&](amrex::Array<Real,NPHASE> const& m,Real& f,Real& g){
        f=1.;Real offset=0.;for(int l=0;l<NPHASE;++l){f-=p.C[l]*m[l];offset+=p.D[l]*m[l];}
        g=(H-offset-pi_min*f)/scale;
    };
    if(project(uniform))return true;
    amrex::Array<Real,NPHASE> a{},b{},candidate{};
    Real lo=0.,hi=1.;bool have_a=false,have_b=false;
    // Dual of max min(f,g) over the product of group simplices. At each
    // dual weight its maximizer assigns a group to one of its pure phases.
    // Opposite one-sided maximizers bracket a feasible convex combination.
    for(int iteration=0;iteration<96;++iteration){
        Real weight=iteration==0?0.:iteration==1?1.:.5*(lo+hi);
        for(int l=0;l<NPHASE;++l)candidate[l]=0.;
        for(int group=0;group<NPHASE;++group)if(mass[group]>0.){
            Real best=-std::numeric_limits<Real>::max();int phase=-1;
            for(int l=0;l<NPHASE;++l)if(p.Physics_Parm.chemical_group[l]==group){
                Real score=-weight*p.C[l]+(1.-weight)*(pi_min*p.C[l]-p.D[l])/scale;
                if(score>best){best=score;phase=l;}
            }
            candidate[phase]=mass[group];
        }
        if(project(candidate))return true;
        Real f,g;constraints(candidate,f,g);
        if(f>0. && g>0.){
            // Retain a small, feasible interior contribution when a pure
            // vertex excludes the phase imposing the pressure-domain bound.
            Real mix=.5;amrex::Array<Real,NPHASE> interior{};
            for(int n=0;n<96;++n){
                for(int l=0;l<NPHASE;++l)interior[l]=(1.-mix)*candidate[l]+mix*uniform[l];
                Real fi,gi;constraints(interior,fi,gi);
                if(fi>0. && gi>0. && project(interior))return true;
                mix*=.5;
            }
        }
        if(f-g<0.){lo=weight;a=candidate;have_a=true;}else{hi=weight;b=candidate;have_b=true;}
    }
    if(!have_a || !have_b)return false;
    Real fa,ga,fb,gb;constraints(a,fa,ga);constraints(b,fb,gb);
    Real denominator=(fb-fa)-(gb-ga);
    if(denominator==0.)return false;
    Real weight=(ga-fa)/denominator;
    if(!(weight>=0. && weight<=1.))return false;
    for(int l=0;l<NPHASE;++l)candidate[l]=(1.-weight)*a[l]+weight*b[l];
    if(project(candidate))return true;
    Real f,g;constraints(candidate,f,g);
    if(f>0. && g>0.){
        Real mix=.5;amrex::Array<Real,NPHASE> interior{};
        for(int n=0;n<96;++n){
            for(int l=0;l<NPHASE;++l)interior[l]=(1.-mix)*candidate[l]+mix*uniform[l];
            Real fi,gi;constraints(interior,fi,gi);
            if(fi>0. && gi>0. && project(interior))return true;
            mix*=.5;
        }
    }
    return false;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void SixEqIE_RelaxChemicalCell(amrex::Array<amrex::Real,NSTATE>& out,
    amrex::Array<amrex::Real,NSTATE> const& input,Parm const& p)
{
    using amrex::Real;
    Real pressure,T;
    const auto status=SixEqIE_TryPressureTemperatureCell(out,input,p,pressure,T);
    if(status!=SixEqIE_PTStatus::Success && !SixEqIE_ChemicalGroupSeed(out,input,p,pressure,T))
        amrex::Abort("Chemical relaxation: conserved substance masses/energy have no admissible equilibrium seed");
    Real initial_entropy;
    if(SixEqIE_ChemicalResidual(out,pressure,T,p,initial_entropy)<=p.Physics_Parm.chemical_tolerance)return;
    for(int sweep=0;sweep<p.Physics_Parm.chemical_max_sweeps;++sweep){
        // Use all members of each connected group: an absent intermediate
        // phase must not block exchange between two other group members.
        for(int i=0;i<NPHASE;++i)for(int j=i+1;j<NPHASE;++j)
            if(p.Physics_Parm.chemical_group[i]==p.Physics_Parm.chemical_group[j] &&
               out[INDEX_Mass1+i]+out[INDEX_Mass1+j]>0.){
                auto point=SixEqIE_ChemicalOptimizePair(i,j,out,p);
                out=point.state;pressure=point.pressure;T=point.T;
            }
        Real entropy;
        const Real error=SixEqIE_ChemicalResidual(out,pressure,T,p,entropy);
        if(entropy<initial_entropy-1024.*std::numeric_limits<Real>::epsilon()*std::max(Real(1.),std::abs(initial_entropy)))
            amrex::Abort("Chemical relaxation: mixture entropy decreased");
        if(error<=p.Physics_Parm.chemical_tolerance)return;
    }
    // Report the stalled state before aborting; never accept a failed flash.
    Real diagnostic_entropy;
    const Real diagnostic_residual=SixEqIE_ChemicalResidual(out,pressure,T,p,diagnostic_entropy);
    printf("Chemical relaxation stalled: p=%.17g T=%.17g residual=%.17g tolerance=%.17g sweeps=%d\n",
           double(pressure),double(T),double(diagnostic_residual),double(p.Physics_Parm.chemical_tolerance),p.Physics_Parm.chemical_max_sweeps);
    for(int l=0;l<NPHASE;++l){
        SixEqIE_ChemicalThermo thermo;
        const bool valid=SixEqIE_ChemicalThermodynamics(thermo,pressure,T,l,p);
        printf("Chemical phase %d: input_alpha=%.17g input_mass=%.17g input_H=%.17g alpha=%.17g mass=%.17g g=%.17g valid=%d\n",
               l,double(input[l]),double(input[INDEX_Mass1+l]),double(input[INDEX_Energy1+l]),double(out[l]),double(out[INDEX_Mass1+l]),double(thermo.g),int(valid));
    }
    for(int n=0;n<NSTATE;++n)printf("Chemical input[%d]=%.17g\n",n,double(input[n]));
    amrex::Abort("Chemical relaxation: coupled substance groups did not converge; see stalled-state diagnostics");
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void SixEqIE_RelaxEquilibriumCell(amrex::Array<amrex::Real,NSTATE>& out,
    amrex::Array<amrex::Real,NSTATE> const& input,Parm const& p)
{
    if(p.Physics_Parm.chemical_relaxation)SixEqIE_RelaxChemicalCell(out,input,p);
    else SixEqIE_RelaxPressureTemperatureCell(out,input,p);
}

AMREX_GPU_HOST inline
void SixEqIE_RelaxEquilibrium(amrex::MultiFab& state,Parm const* p,bool chemical)
{
    // Keep the existing smaller GPU kernel when chemistry is disabled.
    if(!chemical){SixEqIE_RelaxPressureTemperature(state,p);return;}
    for(amrex::MFIter mfi(state,amrex::TilingIfNotGPU());mfi.isValid();++mfi){
        auto a=state.array(mfi);auto box=mfi.tilebox();
        amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE(int i,int j,int k){
            amrex::Array<amrex::Real,NSTATE> in{},out{};
            for(int n=0;n<NSTATE;++n)in[n]=a(i,j,k,n);
            SixEqIE_RelaxChemicalCell(out,in,*p);
            for(int n=0;n<NSTATE;++n)a(i,j,k,n)=out[n];
        });
    }
}
#endif
#endif
