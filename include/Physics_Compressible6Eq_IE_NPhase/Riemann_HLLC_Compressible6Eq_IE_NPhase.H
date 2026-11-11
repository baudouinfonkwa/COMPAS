// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
#ifndef Riemann_HLLC_Compressible6Eq_IE_NPhase_H_
#define Riemann_HLLC_Compressible6Eq_IE_NPhase_H_
#include <Macros.H>
#include <HLLC_ThermalStar_Compressible6Eq_IE_NPhase.H>
#if PHYSICS == SIXEQS_IE_NPHASE

template<typename TF,typename TV,typename TU>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void Riemann_HLLC(TF& flux,TV& velocity,amrex::Real& cmax,
    TU const& left,TU const& right,Parm const& parm,int dir)
{
    using amrex::Real;
    constexpr int momentum=3*NPHASE;
    constexpr int energy=momentum+AMREX_SPACEDIM;
    const auto implementation=parm.Physics_Parm.hllc_implementation;
    const bool physical=implementation==HLLCImplementation::Consistent;
    if(!physical && implementation!=HLLCImplementation::ConservativePart)
        amrex::Abort("SixEq IE HLLC: unknown implementation");
    amrex::Array<Real,NSTATE> VL,VR,FL,FR;
    Conservative2Primitive(VL,left,parm);
    Conservative2Primitive(VR,right,parm);
    Conservative2Flux(FL,left,parm,dir);
    Conservative2Flux(FR,right,parm,dir);
    amrex::Array<Real,NPHASE> alphaL,massL,apL,alphaR,massR,apR;
    Real rhoL=0.,rhoR=0.,pL=0.,pR=0.,absPL=0.,absPR=0.;
    for(int k=0;k<NPHASE;++k){
        alphaL[k]=VL[k];alphaR[k]=VR[k];
        massL[k]=VL[NPHASE+k];massR[k]=VR[NPHASE+k];
        apL[k]=VL[2*NPHASE+k];apR[k]=VR[2*NPHASE+k];
        rhoL+=massL[k];rhoR+=massR[k];pL+=apL[k];pR+=apR[k];
        absPL+=std::abs(apL[k]);absPR+=std::abs(apR[k]);
    }
    Real cL,cR;SoundSpeed(cL,alphaL,massL,apL,parm);SoundSpeed(cR,alphaR,massR,apR,parm);
    const Real uL=VL[momentum+dir],uR=VR[momentum+dir];
    const Real average_u=.5*(uL+uR),average_c=.5*(cL+cR);
    // Keep the pre-existing wave-speed estimates to isolate the star/flux choice.
    // const Real SL=std::min(average_u-average_c,uL-cL);
    // const Real SR=std::max(average_u+average_c,uR+cR);
    const Real SL=std::min(uR-cR,uL-cL);
    const Real SR=std::max(uL+cL,uR+cR);
    Real normL=0.,normR=0.;
    for(int d=0;d<AMREX_SPACEDIM;++d){normL+=VL[momentum+d]*VL[momentum+d];normR+=VR[momentum+d]*VR[momentum+d];}
    cmax=std::max(cL+std::sqrt(normL),cR+std::sqrt(normR));
    if(!std::isfinite(cmax) || !(rhoL>0.) || !(rhoR>0.) || !(SL<SR))
        amrex::Abort("SixEq IE HLLC: invalid density or wave-speed bounds");
    // Supersonic faces do not require a star-state construction.
    if(SL>=0. || SR<=0.){
        const bool choose_left=SL>=0.;
        for(int c=0;c<NSTATE;++c)flux[c]=choose_left?FL[c]:FR[c];
        for(int d=0;d<AMREX_SPACEDIM;++d)velocity[d]=choose_left?VL[momentum+d]:VR[momentum+d];
        return;
    }
    const Real qL=rhoL*(uL-SL),qR=rhoR*(uR-SR);
    if(qL==qR || !std::isfinite(qL) || !std::isfinite(qR))
        amrex::Abort("SixEq IE HLLC: degenerate contact-speed denominator");
    const Real contact=(qL*uL-qR*uR+pL-pR)/(qL-qR);
    if(!(contact>SL && contact<SR) || !std::isfinite(contact))
        amrex::Abort("SixEq IE HLLC: contact speed outside acoustic fan");
    const Real pressure_star=pL+qL*(uL-contact);
    amrex::Array<Real,NSTATE> starL=VL,starR=VR;
    const Real rL=(uL-SL)/(contact-SL),rR=(uR-SR)/(contact-SR);
    for(int k=0;k<NPHASE;++k){
        starL[k]=physical?alphaL[k]:rL*alphaL[k];
        starR[k]=physical?alphaR[k]:rR*alphaR[k];
        starL[NPHASE+k]=rL*massL[k];starR[NPHASE+k]=rR*massR[k];
        starL[2*NPHASE+k]=SixEqIE_HLLC_AlphaPressureStar(alphaL[k],massL[k],apL[k],starL[k],rL,k,parm,std::max(absPL,absPR));
        starR[2*NPHASE+k]=SixEqIE_HLLC_AlphaPressureStar(alphaR[k],massR[k],apR[k],starR[k],rR,k,parm,std::max(absPL,absPR));
    }
    for(int d=0;d<AMREX_SPACEDIM;++d){
        starL[momentum+d]=d==dir?contact:VL[momentum+d];
        starR[momentum+d]=d==dir?contact:VR[momentum+d];
    }
    // Both implementations first produce primitive alpha*p*. All phasic H*
    // values come from the existing EOS conversion, including NASG and D offsets.
    amrex::Array<Real,NSTATE> UL,UR;
    Primitive2Conservative(UL,starL,parm);
    Primitive2Conservative(UR,starR,parm);
    // Preserve the independent conservative total-energy jump. It need not
    // equal sum(H*)+K* for this approximate Riemann fan.
    UL[energy]=((uL-SL)*left[energy]+pL*uL-pressure_star*contact)/(contact-SL);
    UR[energy]=((uR-SR)*right[energy]+pR*uR-pressure_star*contact)/(contact-SR);
    const Real weightL=.5*(1.+Sign(contact)),weightR=1.-weightL;
    for(int c=0;c<NSTATE;++c){
        // Use the same convective flux for nonconservative variables in both
        // modes. Conservative-part stars recover their conservative RH flux.
        if(c<NPHASE || (c>=2*NPHASE && c<3*NPHASE))
            flux[c]=contact*(weightL*UL[c]+weightR*UR[c]);
        else
            flux[c]=weightL*(FL[c]+SL*(UL[c]-left[c]))
                   +weightR*(FR[c]+SR*(UR[c]-right[c]));
    }
    for(int d=0;d<AMREX_SPACEDIM;++d){
        if(d!=dir)velocity[d]=weightL*VL[momentum+d]+weightR*VR[momentum+d];
        else if(physical)velocity[d]=contact;
        else{
            // Conservative-part alpha flux = alpha_donor*u_hat; use the
            // same trace for +alpha div(u_hat) to preserve sum(alpha)=1.
            velocity[d]=weightL*(uL+SL*(rL-1.))+weightR*(uR+SR*(rR-1.));
        }
    }
}
#endif
#endif
