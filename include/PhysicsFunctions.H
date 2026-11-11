// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef PHYSICSFUNCTIONS_
#define PHYSICSFUNCTIONS_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

// Constant-coefficient EOS 0--3: rho*e=(1-b*rho)*(A*p+B)+rho*D.
// All arguments use one consistent unit system. B is the canonical EOS
// intercept, while b is the covolume (the C coefficient).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real AffineEOSSoundSpeed(amrex::Real density,amrex::Real pressure,
    amrex::Real A,amrex::Real B,amrex::Real b)
{
    const auto pi=B/(A+1.),xi=1.-b*density;
    if(!(density>0.) || !(A>0.) || !(b>=0.) || !(xi>0.) || !(pressure+pi>0.) ||
       !std::isfinite(density) || !std::isfinite(pressure) || !std::isfinite(A) ||
       !std::isfinite(B) || !std::isfinite(b))
        amrex::Abort("Shock initialization: inadmissible upstream EOS state");
    return std::sqrt((1.+1./A)*(pressure+pi)/(density*xi));
}

// Generalized overload for ideal gas, SG, constant-reference MG and NASG.
// Mach is measured relative to the upstream fluid. The returned entries are
// density, lab velocity, pressure, total energy density and lab shock speed.
// Reference energy D cancels from the jump but is retained in total energy.
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void StiffenedShockConditions(amrex::Real Ma,amrex::Real u0,
    amrex::Real p0,amrex::Real rho0,amrex::Real A,amrex::Real B,
    amrex::Real b,amrex::Real D,T* ret_param)
{
    const auto c=AffineEOSSoundSpeed(rho0,p0,A,B,b);
    if(!(Ma>0.) || !std::isfinite(Ma) || !std::isfinite(u0) || !std::isfinite(D))
        amrex::Abort("Shock initialization: invalid Mach, velocity or reference energy");
    const auto g=1.+1./A,pi=B/(A+1.),v=1./rho0,M2=Ma*Ma;
    const auto post=b+(v-b)*((g-1.)+2./M2)/(g+1.);
    const auto rho=1./post,u=u0+Ma*c*(1.-post/v);
    const auto p=p0+2.*g*(p0+pi)*(M2-1.)/(g+1.);
    const auto energy=(1.-b*rho)*(A*p+B)+rho*D+.5*rho*u*u;
    const auto speed=u0+Ma*c;
    if(!(rho>0.) || !(1.-b*rho>0.) || !(p+pi>0.) ||
       !std::isfinite(rho) || !std::isfinite(u) || !std::isfinite(p) ||
       !std::isfinite(energy) || !std::isfinite(speed))
        amrex::Abort("Shock initialization: inadmissible downstream EOS state");
    ret_param[0]=rho;ret_param[1]=u;ret_param[2]=p;
    ret_param[3]=energy;ret_param[4]=speed;
}

template <typename T>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
StiffenedShockConditions(amrex::Real Ma, 
                         amrex::Real u0, 
                         amrex::Real p0, 
                         amrex::Real rho0, 
                         amrex::Real n, 
                         amrex::Real B,
                         T * ret_param)
{
    if(!(n>1.) || !std::isfinite(n))amrex::Abort("Shock initialization requires gamma>1");
    // Backward-compatible gamma/p_infinity interface for existing callers.
    StiffenedShockConditions(Ma,u0,p0,rho0,1./(n-1.),n*B/(n-1.),0.,0.,ret_param);

}

#endif /* PHYSICSFUNCTIONS_ */
