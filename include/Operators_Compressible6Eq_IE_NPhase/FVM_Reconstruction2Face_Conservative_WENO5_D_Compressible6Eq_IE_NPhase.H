// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#ifndef FVM_Reconstruction2Face_Conservative_WENO5_D_Compressible6Eq_IE_NPhase_H_
#define FVM_Reconstruction2Face_Conservative_WENO5_D_Compressible6Eq_IE_NPhase_H_

#include<Macros.H>
#if (PHYSICS==SIXEQS_IE_NPHASE)

// Component-wise fifth-order WENO with an admissible two-face polynomial.
// The convex limiter is inactive on smooth states inside the EOS domain.
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void FVM_Reconstruction2Face_Conservative_WENO5_D(amrex::Array<amrex::Real,NSTATE>& Conservative_Face,
                                                  amrex::Array<amrex::Real,NSTATE> const& Conservative_mm,
                                                  amrex::Array<amrex::Real,NSTATE> const& Conservative_m,
                                                  amrex::Array<amrex::Real,NSTATE> const& Conservative_c,
                                                  amrex::Array<amrex::Real,NSTATE> const& Conservative_p,
                                                  amrex::Array<amrex::Real,NSTATE> const& Conservative_pp,
                                                  Parm const& EOS_Parameter,
                                                  amrex::Real const& dx,
                                                  int ID_Face)
{
    amrex::Array<amrex::Real,NSTATE> stencil[5];
    amrex::Array<amrex::Real,NSTATE> const* u[5]={&Conservative_mm,&Conservative_m,&Conservative_c,&Conservative_p,&Conservative_pp};
    for(int j=0;j<5;++j){
        Conservative2Primitive(stencil[j],*u[j],EOS_Parameter);
        for(int l=0;l<NPHASE;++l)
            stencil[j][2*NPHASE+l]=stencil[j][l]>0.?stencil[j][2*NPHASE+l]/stencil[j][l]:std::max(stencil[j][NSTATE-1],-EOS_Parameter.B[l]/(EOS_Parameter.A[l]+1.));
    }
    amrex::Array<amrex::Real,NSTATE> face[2],limited[2];
    for(int side=0;side<2;++side)for(int n=0;n<NSTATE;++n)
        Reconstruction2Face_WENO5_D(face[side][n],stencil[0][n],stencil[1][n],stencil[2][n],stencil[3][n],stencil[4][n],side==0?1:-1);
    SixEqIE_LimitReconstructedPair(limited[0],limited[1],face[0],face[1],Conservative_c,EOS_Parameter,1./12.);
    Conservative_Face=limited[ID_Face==1?0:1];
}

#endif
#endif