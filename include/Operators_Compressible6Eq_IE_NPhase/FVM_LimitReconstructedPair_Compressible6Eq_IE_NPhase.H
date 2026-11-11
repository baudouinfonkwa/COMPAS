// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_LimitReconstructedPair_Compressible6Eq_IE_NPhase_H_
#define FVM_LimitReconstructedPair_Compressible6Eq_IE_NPhase_H_
#include <limits>
#if PHYSICS == SIXEQS_IE_NPHASE

// A single convex factor for both faces preserves the cell average and the
// phase partition. Include the interior quadrature state, not just the faces.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void SixEqIE_ReconstructionFactor(amrex::Real& theta,amrex::Real center,
    amrex::Real plus,amrex::Real minus,amrex::Real weight)
{
    const auto interior=(center-weight*(plus+minus))/(1.-2.*weight);
    const auto minimum=std::min({plus,minus,interior});
    if(!std::isfinite(center) || !std::isfinite(minimum) || center<0.)
        amrex::Abort("SixEq IE reconstruction: inadmissible cell average");
    // The ratio can round to exactly one when a negative undershoot is
    // tiny relative to the average. Still retreat into the admissible set.
    if(minimum<0.)theta=std::min(theta,(1.-64.*std::numeric_limits<amrex::Real>::epsilon())*center/(center-minimum));
}

// The primitive energy slots here hold INTRINSIC phasic pressures. Reconstruct
// p_k, then form alpha_k*p_k with the limited alpha, so absent phases carry no
// spurious pressure energy. This is reconstruction limiting, not cell clipping.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void SixEqIE_LimitReconstructedPair(amrex::Array<amrex::Real,NSTATE>& plus,
    amrex::Array<amrex::Real,NSTATE>& minus,
    amrex::Array<amrex::Real,NSTATE> vp,
    amrex::Array<amrex::Real,NSTATE> vm,
    amrex::Array<amrex::Real,NSTATE> const& center,Parm const& p,amrex::Real weight)
{
    using amrex::Real;
    amrex::Array<Real,NSTATE> vc;
    Conservative2Primitive(vc,center,p);
    // Intrinsic pressure is unobservable at exactly zero volume. Extend it
    // inside that phase's EOS domain; its stored alpha*p stays exactly zero.
    for(int l=0;l<NPHASE;++l)
        vc[2*NPHASE+l]=vc[l]>0.?vc[2*NPHASE+l]/vc[l]:std::max(vc[NSTATE-1],-p.B[l]/(p.A[l]+1.));
    Real sp=0.,sm=0.;for(int l=0;l<NPHASE;++l){sp+=vp[l];sm+=vm[l];}
    if(!(sp>0.) || !(sm>0.) || !std::isfinite(sp) || !std::isfinite(sm)){
        plus=center;minus=center;return;
    }
    for(int l=0;l<2*NPHASE;++l){vp[l]/=sp;vm[l]/=sm;}
    Real theta=1.;
    for(int l=0;l<NPHASE;++l){
        SixEqIE_ReconstructionFactor(theta,vc[l],vp[l],vm[l],weight);
        SixEqIE_ReconstructionFactor(theta,vc[NPHASE+l],vp[NPHASE+l],vm[NPHASE+l],weight);
        const auto b=p.C[l],pi=p.B[l]/(p.A[l]+1.);
        SixEqIE_ReconstructionFactor(theta,vc[l]-b*vc[NPHASE+l],
            vp[l]-b*vp[NPHASE+l],vm[l]-b*vm[NPHASE+l],weight);
        // Pressure is not a conserved polynomial; only its face domain is needed.
        const auto pc=vc[2*NPHASE+l]+pi;
        for(Real pf:{vp[2*NPHASE+l]+pi,vm[2*NPHASE+l]+pi}){
            if(!std::isfinite(pc) || !std::isfinite(pf) || pc<0.)
                amrex::Abort("SixEq IE reconstruction: invalid cell pressure domain");
            if(pf<0.)theta=std::min(theta,(1.-64.*std::numeric_limits<Real>::epsilon())*pc/(pc-pf));
        }
    }
    if(theta<1.)theta*=1.-64.*std::numeric_limits<Real>::epsilon();
    for(int n=0;n<NSTATE;++n){vp[n]=(1.-theta)*vc[n]+theta*vp[n];vm[n]=(1.-theta)*vc[n]+theta*vm[n];}
    for(int l=0;l<NPHASE;++l){vp[2*NPHASE+l]*=vp[l];vm[2*NPHASE+l]*=vm[l];}
    Primitive2Conservative(plus,vp,p);Primitive2Conservative(minus,vm,p);
    // The EOS thermal domain is linear in the weighted conservative variables:
    // H - m D - (alpha-b*m)*pi = (alpha-b*m)*A*(p+pi) >= 0.
    theta=1.;
    for(int l=0;l<NPHASE;++l){
        const auto b=p.C[l],pi=p.B[l]/(p.A[l]+1.);
        const auto thermal=[&](amrex::Array<Real,NSTATE> const& u){
            return u[2*NPHASE+l]-u[NPHASE+l]*p.D[l]-(u[l]-b*u[NPHASE+l])*pi;
        };
        SixEqIE_ReconstructionFactor(theta,thermal(center),thermal(plus),thermal(minus),weight);
    }
    if(theta<1.)theta*=1.-64.*std::numeric_limits<Real>::epsilon();
    for(int n=0;n<NSTATE;++n){plus[n]=(1.-theta)*center[n]+theta*plus[n];minus[n]=(1.-theta)*center[n]+theta*minus[n];}
}
#endif
#endif
