// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef IC_H_
#define IC_H_

#include <AMReX_ParmParse.H>
#include <Parm.H>
#include "PhysicsFunctions.H"
#include <Tools.H>


// =======================================
// Initial condition
//
// Shock-bubble interaction: NPHASE=2 water/air, or NPHASE=3 water/vapor/air.
// A radius-0.5 cavity is centered on the symmetry plane. Three-phase raw gas
// volumes are split with prob.bubble_vapor_fraction. With a positive
// prob.initial_temperature, all upstream densities come from EOS(P0,T).
// Otherwise the two-phase case retains its prescribed rho0/rho_air densities.
// A frozen-water RH shock of Mach prob.Mach starts at x=-0.85. The reference
// density, pressure, velocity and length are rho0, P0, sqrt(P0/rho0) and the
// bubble diameter. Per-phase initialization floors do not constrain transport
// or chemical equilibrium; chemistry can set phase mass/volume exactly to zero.
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void initial_condition(amrex::Array<amrex::Real, NSTATE>& Conservative,
                       bool& sharp,
                       amrex::Real x,
                       amrex::Real y,
                       amrex::Real z,
                       AMREX_D_DECL(amrex::Real dx,
                                    amrex::Real dy,
                                    amrex::Real dz),
                       AMREX_D_DECL(amrex::Real dx_FinestLevel,
                                    amrex::Real dy_FinestLevel,
                                    amrex::Real dz_FinestLevel),
                       Parm const& parm)
{
    // Prevent unused-variable warnings
    sharp = false;
    amrex::ignore_unused(AMREX_D_DECL(dx, dy, dz));

    // ---------------------------------------
    // Interface: tanh profile of width eta about the bubble surface
    // ---------------------------------------
    const amrex::Real eta = parm.PhaseField_Parm.Eta_Multiplier * dx_FinestLevel;
    const amrex::Real Psi = 0.5 - std::sqrt(x*x + y*y + z*z);

    const amrex::Real indicator=eta>0. ? std::tanh(-Psi/(std::sqrt(2.)*eta)) : (Psi<0.?1.:-1.);
    amrex::Array<amrex::Real,NPHASE> alpha{0.5*(1.+indicator),0.5*(1.-indicator)};
    if constexpr(NPHASE==3){
        const amrex::Real bubble=alpha[1];
        alpha[1]=bubble*parm.bubble_vapor_fraction;
        alpha[2]=bubble*(1.-parm.bubble_vapor_fraction);
    }
    parm.ApplyInitialAlphaFloor(alpha);

    // ---------------------------------------
    // Fluids at rest, and the post-shock state of the water behind the shock (SI units)
    // ---------------------------------------
    const amrex::Real x_Shock = -0.85;

    amrex::Real rho1 = parm.upstream_rho*parm.rho0;
    amrex::Real P1   = parm.P0;
    amrex::Real P2   = parm.P0;
    amrex::Real u    = 0.0;
    amrex::Real v    = 0.0;
    amrex::Real w    = 0.0;

    if (x < x_Shock) {
        rho1 = parm.shock_rho * parm.rho0;
        u    = parm.shock_u * parm.u0;
        P1   = parm.shock_p * parm.P0;
    }

    // ---------------------------------------
    // Scale, convert to primitive and then conservative
    // ---------------------------------------
    amrex::Array<amrex::Real, NPHASE> Alphas, Alpha_rhos, Alpha_Ps;
    for(int l=0;l<NPHASE;++l){
        Alphas[l]=alpha[l];
        const amrex::Real density=l==0 ? rho1 :
            (NPHASE==2 && parm.initial_temperature==0. ? parm.rho_air : parm.upstream_phase_rho[l]*parm.rho0);
        Alpha_rhos[l]=alpha[l]*density/parm.rho0;
        Alpha_Ps[l]=alpha[l]*(l==0 ? P1 : P2)/parm.P0;
    }
    AMREX_D_TERM(u /= parm.u0;, v /= parm.u0;, w /= parm.u0;)
    amrex::Real P=0.;for(int l=0;l<NPHASE;++l)P+=Alpha_Ps[l];

    amrex::Array<amrex::Real, NSTATE> Primitive;
    State2Primitive(Primitive,
                    Alphas,
                    Alpha_rhos,
                    Alpha_Ps,
       AMREX_D_DECL(u, v, w),
                    P);
    Primitive2Conservative(Conservative,
                           Primitive,
                           parm);
}


// =======================================
// Spatio-temporal refinement tagging
//
// Tags every cell within a distance 2 of the bubble center for refinement to the finest
// level at all times, in addition to the gradient criteria in the inputs file.
// =======================================

#define SPATIO_TEMPORAL_TAG_USER_FUNC_
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void spatio_temporal_tag(int& TagVal,
                         int lev,
                         amrex::Real x,
                         amrex::Real y,
                         amrex::Real z,
                         amrex::Real t,
                         Parm const& parm)
{
    amrex::ignore_unused(lev, t, parm);

    const amrex::Real r = std::sqrt(x*x + y*y + z*z);
    if (r < 2.0) {
        TagVal = TagBox::SET;
    }
}


// =======================================
// In-situ diagnostics
//
// With run.user_output_int = N, COMPAS calls this every N coarse steps. Each call appends
// one line to <plot directory>/bubble.csv with
//   V        gas-cavity volume, int sum_{k>0} alpha_k dV (upper-half area in 2D)
//   x_c      x-centroid of the air, int alpha_2 x dV / V
//   u_c      alpha_2-weighted mean x-velocity of the air, int alpha_2 u dV / V
//   x_up     upstream interface on the axis, the smallest x where alpha_2 crosses 0.5
//   x_down   downstream interface on the axis, the largest x where alpha_2 crosses 0.5
// The integrals are over the composite grid; the axial crossings are taken on the finest
// level, in the row of cells next to y = 0 (and z = 0 in 3D).
// Lengths, times and velocities are in the scaled units of the inputs file.
// =======================================

#define USER_OUTPUT_FUNC
template <class T>
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void UserOutputFunction(T& o)
{
    using namespace UserOutput;
    const Real big = 1.0e30;

    // One ghost cell, for the +x neighbor in the crossings
    const auto v = Measure(o, {Sum, Sum, Sum, Min, Max}, 1,
    [=] AMREX_GPU_DEVICE (Cell const& c, State const& U, Parm const&) -> Values<5>
    {
        const int  i = c.i, j = c.j, k = c.k;
        // Whole gas cavity: vapor plus air for the three-phase build.
        Real a2=0.,a2e=0.;
        for(int l=1;l<NPHASE;++l){a2+=U(i,j,k,INDEX_VolumeFraction1+l);a2e+=U(i+1,j,k,INDEX_VolumeFraction1+l);}
        Real rho = 0.0;
        for (int n = 0; n < NPHASE; ++n) { rho += U(i,j,k,INDEX_Mass1+n); }
        const Real u   = U(i,j,k,INDEX_MomentumX) / rho;

        // alpha_2 = 0.5 crossing toward the +x neighbor, on the axis
        bool axis = c.finest && (j == c.domain.smallEnd(1));
#if (AMREX_SPACEDIM == 3)
        axis = axis && (c.x[2] - 0.5 * c.dx[2] <= 0.0) && (c.x[2] + 0.5 * c.dx[2] > 0.0);
#endif
        Real x_up = big, x_down = -big;
        if (axis && (a2 - 0.5) * (a2e - 0.5) <= 0.0 && a2 != a2e) {
            x_up = x_down = c.x[0] + c.dx[0] * (0.5 - a2) / (a2e - a2);
        }

        return {a2 * c.dV, a2 * c.x[0] * c.dV, a2 * u * c.dV, x_up, x_down};
    });

    const Real V = v[0];
    WriteRow(o, "bubble.csv", {{"V", V}, {"x_c", (V > 0.0) ? v[1] / V : 0.0}, {"u_c", (V > 0.0) ? v[2] / V : 0.0},
                                {"x_up", v[3]}, {"x_down", v[4]}});
}

#endif  // IC_H_
