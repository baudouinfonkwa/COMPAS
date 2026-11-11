// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
#ifndef HLLC_ThermalStar_Compressible6Eq_IE_NPhase_H_
#define HLLC_ThermalStar_Compressible6Eq_IE_NPhase_H_
#include <Macros.H>
#include <cmath>
#include <limits>
#if PHYSICS == SIXEQS_IE_NPHASE

// Returns alpha* p*. H* is deliberately recovered separately by the canonical
// EOS function. The current EOS family is e=(A p+B)(v-b)+D, with b=C[k].
// For physical acoustic stars alpha*=alpha and m*=r m, chi=(v*-b)/(v-b).
// (p*+pi)/(p+pi) = [A-(chi-1)/2]/[A*chi+(chi-1)/2].
// Thus alpha*p*=(alpha*/alpha)*[(alpha*p+alpha*pi)*f(chi)-alpha*pi].
// The caller supplies alpha*. This function contains no implementation switch.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real SixEqIE_HLLC_AlphaPressureStar(
    amrex::Real alpha,amrex::Real mass,amrex::Real alpha_p,
    amrex::Real alpha_star,amrex::Real compression,int phase,Parm const& parm,
    amrex::Real pressure_scale=0.)
{
    if(!(compression>0.) || !std::isfinite(compression))
        amrex::Abort("SixEq IE HLLC: invalid star compression ratio");
    if(!std::isfinite(alpha) || !std::isfinite(mass) || !std::isfinite(alpha_p)
       || !std::isfinite(alpha_star))
        amrex::Abort("SixEq IE HLLC: invalid incoming phase state");
    const amrex::Real A=parm.A[phase],B=parm.B[phase],b=parm.C[phase];
    if(!(A>0.) || !(b>=0.) || !std::isfinite(A) || !std::isfinite(B) || !std::isfinite(b))
        amrex::Abort("SixEq IE HLLC: invalid EOS coefficients");
    if(alpha==0.){
        // Reconstruction/reflux can leave an ulp-scale alpha*p at exactly
        // zero alpha. Use the finite-pressure absent-phase limit only within
        // roundoff of the supplied face pressure scale, not an alpha cutoff.
        const amrex::Real tolerance=64.*std::numeric_limits<amrex::Real>::epsilon()*pressure_scale;
        if(std::abs(alpha_p)>tolerance || alpha_star!=0. || (b>0. && mass!=0.))
            amrex::Abort("SixEq IE HLLC: inconsistent zero-volume phase");
        return 0.;
    }
    const amrex::Real phi=alpha-b*mass;
    const amrex::Real shifted_pressure=alpha_p+alpha*(B/(A+1.));
    // For b=0 retain the existing EOS's algebraic extension to signed roundoff
    // traces; do not clip alpha, mass, or alpha*p. NASG requires positive
    // available volume for every nonzero phase. Zero temperature is a limit.
    if(!std::isfinite(phi) || !std::isfinite(shifted_pressure)
       || (b>0. && (!(phi>0.) || !(alpha>0.) || mass<0.))
       || (alpha>0. && shifted_pressure<0.) || (alpha<0. && shifted_pressure>0.))
        amrex::Abort("SixEq IE HLLC: incoming phase violates available-volume or temperature bound");
    const amrex::Real phi_star=alpha_star-b*(compression*mass);
    if(!std::isfinite(phi_star) || (b>0. && !(phi_star>0.))
       || !(alpha_star/alpha>0.))
        amrex::Abort("SixEq IE HLLC: invalid star available volume or fraction ratio");
    // chi=(phi*/phi)/r, lambda=alpha*/alpha. No intrinsic phase density
    // is recovered, and no products of two trace fractions are required.
    const amrex::Real d=(phi_star/phi)/compression-1.;
    const amrex::Real denominator=A+(A+.5)*d;
    const amrex::Real numerator=A-.5*d;
    if(!(denominator>0.) || !(numerator>0.) || !std::isfinite(denominator)
       || !std::isfinite(numerator))
        amrex::Abort("SixEq IE HLLC: thermal Hugoniot pole or inadmissible expansion");
    const amrex::Real result=(alpha_star/alpha)*
        (alpha_p-((A+1.)*d/denominator)*shifted_pressure);
    if(!std::isfinite(result))amrex::Abort("SixEq IE HLLC: nonfinite thermal star pressure");
    return result;
}
#endif
#endif
