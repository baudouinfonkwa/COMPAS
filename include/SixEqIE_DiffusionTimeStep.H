// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef SixEqIE_DiffusionTimeStep_H_
#define SixEqIE_DiffusionTimeStep_H_
#include <AMReX_REAL.H>
#include <AMReX.H>
#include <cmath>
#include <limits>

// run.vnn is a fraction (0,1] of the multidimensional Forward Euler scalar
// diffusion limit. SSP-RK2/RK3 retain that limit. This is the diffusion bound;
// it does not prove stability/positivity of the combined nonlinear operators.
inline amrex::Real SixEqIE_DiffusionTimeStep(amrex::Real diffusivity,
    amrex::Real safety, amrex::Real const* dx, int dimensions)
{
    if (!std::isfinite(safety) || !(safety>0.) || safety>1.)
        amrex::Abort("SIXEQS_IE_NPHASE diffusion requires 0 < run.vnn <= 1 (finite)");
    if (!std::isfinite(diffusivity) || diffusivity<0.)
        amrex::Abort("SIXEQS_IE_NPHASE: invalid diffusion timestep coefficient");
    amrex::Real inverse_square_sum=0.;
    for(int d=0;d<dimensions;++d)inverse_square_sum+=1./(dx[d]*dx[d]);
    return diffusivity>0. ? safety/(2.*diffusivity*inverse_square_sum)
                         : std::numeric_limits<amrex::Real>::max();
}

// Positive WENO reconstruction includes an interior quadrature state. The
// unsplit FE update must respect its smallest weight in the sum of directions.
inline amrex::Real SixEqIE_WENOTimeStep(amrex::Real speed,amrex::Real const* dx,
    int dimensions,int weno,int limiter,int polynomial)
{
    if(!weno || limiter<4 || limiter>6 || !(speed>0.))
        return std::numeric_limits<amrex::Real>::max();
    const amrex::Real weight=(limiter==4 || polynomial==2)?1./6.:1./12.;
    amrex::Real inverse_sum=0.;for(int d=0;d<dimensions;++d)inverse_sum+=1./dx[d];
    return .9*weight/(speed*inverse_sum);
}
#endif
