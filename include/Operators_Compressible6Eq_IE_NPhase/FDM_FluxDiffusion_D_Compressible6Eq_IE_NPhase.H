// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FDM_FluxDiffusion_D_Compressible6Eq_IE_NPhase_H_
#define FDM_FluxDiffusion_D_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE && DIFFUSION == true)

// =======================================
// Face quantities of FluxDiffusion from the cell values of FDM_Conservative2FluxDiffusion_D:
//   - the face velocity and velocity gradient (FDM_FaceVelocityGradient_D, V[b][s] the
//     velocity component b at the stencil point s);
//   - the face viscosities, the harmonic mean of the mixture viscosities of the cells m and c
//     (PhaseDiffusion_FaceViscosity_D);
//   - per phase, the temperature difference across the face, the harmonic face weight and
//     the thermal diffusivity bound for the time step (where the phase conducts).
// =======================================

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FDM_FluxDiffusion_D(amrex::Array<amrex::Real,NSTATE>& Flux,
                         amrex::Real& dmax,
                         amrex::Array<amrex::Array<amrex::Real,FaceStencil::N>,AMREX_SPACEDIM> const& V,
                         amrex::Real const& mu_m,  amrex::Real const& mu_c,
                         amrex::Real const& mub_m, amrex::Real const& mub_c,
                         amrex::Real const& rho_m, amrex::Real const& rho_c,
                         amrex::Array<amrex::Real,NPHASE> const& Alpha_m,
                         amrex::Array<amrex::Real,NPHASE> const& Alpha_c,
                         amrex::Array<amrex::Real,NPHASE> const& Alpha_rho_m,
                         amrex::Array<amrex::Real,NPHASE> const& Alpha_rho_c,
                         amrex::Array<amrex::Real,NPHASE> const& T_m,
                         amrex::Array<amrex::Real,NPHASE> const& T_c,
                         amrex::Array<bool,NPHASE> const& Present_m,
                         amrex::Array<bool,NPHASE> const& Present_c,
            AMREX_D_DECL(amrex::Real const& dx,
                         amrex::Real const& dy,
                         amrex::Real const& dz),
                         Parm const& EOS_Parameter,
                         int ID_Dimension,
                         SixEqIEDiffusionFace* face_data=nullptr)
{
    amrex::Array<amrex::Real,AMREX_SPACEDIM> Vel;
    amrex::Array<amrex::Array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> Grad_Vel;
    FDM_FaceVelocityGradient_D(Vel,
                               Grad_Vel,
                               V,
                  AMREX_D_DECL(dx, dy, dz),
                               ID_Dimension);

    if(EOS_Parameter.Physics_Parm.compatible_diffusion){
        using amrex::Real;
        const int n=ID_Dimension;
        const amrex::Array<Real,AMREX_SPACEDIM> widths={AMREX_D_DECL(dx,dy,dz)};
        const Real h=widths[n];
        const int minus=FaceStencil::Minus(n);
        const Real ml=std::max(mu_m,Real(0.)),mr=std::max(mu_c,Real(0.));
        const Real bl=std::max(mub_m,Real(0.)),br=std::max(mub_c,Real(0.));
        const Real etaL=4.*ml/3.+bl,etaR=4.*mr/3.+br;
        SixEqIEDiffusionFace data;
        for(int b=0;b<AMREX_SPACEDIM;++b){
            const Real kl=b==n?etaL:ml,kr=b==n?etaR:mr;
            Real crossL=ml*Grad_Vel[b][n],crossR=mr*Grad_Vel[b][n];
            if(b==n){
                Real transverse_div=0.;
                for(int a=0;a<AMREX_SPACEDIM;++a)if(a!=n)transverse_div+=Grad_Vel[a][a];
                crossL=(bl-2.*ml/3.)*transverse_div;
                crossR=(br-2.*mr/3.)*transverse_div;
            }
            const Real ul=V[b][minus],ur=V[b][FaceStencil::c],sum=kl+kr;
            Real uf=.5*(ul+ur),traction=0.;
            if(sum>0.){
                const Real wl=kl/sum,wr=kr/sum;
                uf=wl*ul+wr*ur+.5*h*(crossR-crossL)/sum;
                traction=(2.*wl*kr)*(ur-ul)/h+wr*crossL+wl*crossR;
            }
            data.traction[b]=traction;data.work+=uf*traction;
        }
        // Perfect thermal contact: shared cell temperature (heat-capacity
        // weighted for a supplied nonequilibrium initial/AMR state), with
        // harmonic mixture conductivity across two half-cell resistances.
        Real capL=0.,capR=0.,TL=0.,TR=0.,kL=0.,kR=0.;
        for(int l=0;l<NPHASE;++l){
            const Real cl=Alpha_rho_m[l]*EOS_Parameter.cv[l];
            const Real cr=Alpha_rho_c[l]*EOS_Parameter.cv[l];
            capL+=cl;capR+=cr;TL+=cl*T_m[l];TR+=cr*T_c[l];
            kL+=Alpha_m[l]*EOS_Parameter.kappa[l];
            kR+=Alpha_c[l]*EOS_Parameter.kappa[l];
        }
        if(!(capL>0.) || !(capR>0.))amrex::Abort("SixEq IE compatible diffusion: invalid heat capacity");
        TL/=capL;TR/=capR;
        Real conductivity;PhaseDiffusion_FaceViscosity_D(conductivity,kL,kR);
        data.heat=conductivity*(TR-TL)/h;
        for(int c=0;c<NSTATE;++c)Flux[c]=0.;
        for(int b=0;b<AMREX_SPACEDIM;++b)Flux[INDEX_MomentumX+b]=data.traction[b];
        Flux[INDEX_Total_Energy]=data.work+data.heat;
        // Phasic thermal diffusion is a cell heat-capacity partition of the
        // shared energy flux divergence, not a flux into a massless phase.
        dmax=std::max(std::max(etaL,etaR)/std::min(rho_m,rho_c),
                      conductivity/std::min(capL,capR));
        if(face_data)*face_data=data;
        return;
    }

    amrex::Real mu, mub;
    PhaseDiffusion_FaceViscosity_D(mu, mu_m, mu_c);
    PhaseDiffusion_FaceViscosity_D(mub, mub_m, mub_c);

    amrex::Array<amrex::Real,AMREX_SPACEDIM> const dX = {AMREX_D_DECL(dx, dy, dz)};
    amrex::Array<amrex::Real,NPHASE> Weight, Grad_T;
    amrex::Real Chi_T = 0.0;
    for (int i_Phase = 0; i_Phase < NPHASE; ++i_Phase){
        PhaseDiffusion_FaceWeight_D(Weight[i_Phase], Alpha_m[i_Phase], Alpha_c[i_Phase], Present_m[i_Phase], Present_c[i_Phase]);
        Grad_T[i_Phase] = 0.0;
        if (Weight[i_Phase] > 0.0){
            FDM_Gradient_D(Grad_T[i_Phase], T_m[i_Phase], T_c[i_Phase], dX[ID_Dimension]);
            // The face weight over the smaller volume fraction (between 1 and 2) bounds the
            // explicit limit of both cells
            amrex::Real rho_k = std::min(Alpha_rho_m[i_Phase]/Alpha_m[i_Phase], Alpha_rho_c[i_Phase]/Alpha_c[i_Phase]);
            amrex::Real Ratio = Weight[i_Phase]/std::min(Alpha_m[i_Phase], Alpha_c[i_Phase]);
            Chi_T = std::max(Chi_T, Ratio*EOS_Parameter.kappa[i_Phase]/(rho_k*EOS_Parameter.cv[i_Phase]));
        }
    }

    amrex::Real rho = std::min(rho_m, rho_c);

    FluxDiffusion(Flux,
                  dmax,
                  Vel,
                  Grad_Vel,
                  Weight,
                  Grad_T,
                  mu,
                  mub,
                  rho,
                  Chi_T,
                  EOS_Parameter,
                  ID_Dimension);
    if(face_data){
        face_data->heat=0.;
        for(int l=0;l<NPHASE;++l)face_data->heat+=Flux[INDEX_Energy1+l];
        face_data->work=Flux[INDEX_Total_Energy]-face_data->heat;
        for(int b=0;b<AMREX_SPACEDIM;++b)face_data->traction[b]=Flux[INDEX_MomentumX+b];
    }
}

#endif
#endif
