// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.
#include <COMPAS.H>
#include <random>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using Real=amrex::Real;
using State=amrex::Array<Real,NSTATE>;
static_assert(NPHASE==5 && AMREX_SPACEDIM==2,"Use the five-phase 2D IE advection build");
void require(bool ok,char const* text){if(!ok)throw std::runtime_error(text);}
void reconstruct(State& face,State const* s,Parm const& p,int scheme,int side){
    if(scheme==3)FVM_Reconstruction2Face_Conservative_WENO3B_D(face,s[1],s[2],s[3],p,1.,side);
    else if(scheme==5)FVM_Reconstruction2Face_Conservative_WENO5_D(face,s[0],s[1],s[2],s[3],s[4],p,1.,side);
    else FVM_Reconstruction2Face_Conservative_WENO5B_D(face,s[0],s[1],s[2],s[3],s[4],p,1.,side,4);
}
int main(int argc,char**argv){amrex::Initialize(argc,argv);{
    Parm p{};std::mt19937 rng(98241);std::uniform_real_distribution<Real> random(0.,1.);
    long faces=0;Real worstP=0.,worstA=0.;
    for(int eos=0;eos<4;++eos)for(int sample=0;sample<120;++sample){
        p.EOS=eos;
        for(int l=0;l<NPHASE;++l){p.A[l]=1.5+.3*l;p.B[l]=eos==0?0.:(eos==2?-.2:2.)*(p.A[l]+1.);p.C[l]=eos==3?.02*(l+1):0.;p.D[l]=eos>=2?.3*(l-1):0.;p.cv[l]=1.+.1*l;}
        State s[5];
        for(int j=0;j<5;++j){
            State v{};Real sum=0.;
            for(int l=0;l<NPHASE;++l){v[l]=sample%3==0?(l==(j+sample)%NPHASE?1.:0.):(random(rng)<.3?1e-100:random(rng));sum+=v[l];}
            for(int l=0;l<NPHASE;++l){v[l]/=sum;v[NPHASE+l]=v[l]*(1.+.2*l);v[2*NPHASE+l]=2.*v[l];}
            v[INDEX_MomentumX]=.3;v[INDEX_MomentumY]=-.2;Primitive2Conservative(s[j],v,p);
        }
        for(int scheme:{3,5,6}){
            State face[2];for(int side=0;side<2;++side){
                reconstruct(face[side],s,p,scheme,side==0?1:-1);State v;Conservative2Primitive(v,face[side],p);
                Real sum=0.;for(int l=0;l<NPHASE;++l){
                    if(face[side][l]<0. || face[side][NPHASE+l]<0.)std::cerr<<"negative eos="<<eos<<" sample="<<sample<<" scheme="<<scheme<<" phase="<<l<<" alpha="<<std::setprecision(17)<<face[side][l]<<" mass="<<face[side][NPHASE+l]<<"\n";require(face[side][l]>=0. && face[side][NPHASE+l]>=0.,"negative reconstructed phase");sum+=face[side][l];
                    const auto phi=face[side][l]-p.C[l]*face[side][NPHASE+l];
                    require(phi>=0.,"negative reconstructed available volume");
                    if(face[side][l]>0.){auto err=std::abs(v[2*NPHASE+l]/face[side][l]-2.);worstP=std::max(worstP,err);if(err>=1e-11)std::cerr<<"eos="<<eos<<" sample="<<sample<<" scheme="<<scheme<<" phase="<<l<<" alpha="<<std::setprecision(17)<<face[side][l]<<" pressure="<<v[2*NPHASE+l]/face[side][l]<<"\n";require(err<1e-11,"constant phasic pressure not preserved");}
                    else require(v[2*NPHASE+l]==0.,"absent phase has pressure energy");
                }
                worstA=std::max(worstA,std::abs(sum-1.));require(std::abs(sum-1.)<3e-14,"face partition not closed");++faces;
            }
            // The independent EOS temperature numerator of the interior
            // quadrature state must also be nonnegative, including absent phases.
            const Real w=scheme==3?1./6.:1./12.;
            for(int l=0;l<NPHASE;++l){
                const auto a=(s[2][l]-w*(face[0][l]+face[1][l]))/(1.-2.*w);
                const auto m=(s[2][NPHASE+l]-w*(face[0][NPHASE+l]+face[1][NPHASE+l]))/(1.-2.*w);
                const auto H=(s[2][2*NPHASE+l]-w*(face[0][2*NPHASE+l]+face[1][2*NPHASE+l]))/(1.-2.*w);
                const auto phi=a-p.C[l]*m;
                const auto thermal=H-m*p.D[l]-phi*p.B[l]/(p.A[l]+1.);
                require(a>=-2e-15 && m>=-2e-15 && phi>=-2e-15 && thermal>=-1e-13,"interior quadrature state outside EOS domain");
            }
        }
    }
    // A stretched pure liquid is admissible although an absent gas would
    // not admit its pressure. The absent phase must stay exactly absent.
    p.EOS=3;
    for(int l=0;l<NPHASE;++l){p.A[l]=2.5;p.B[l]=l==0?10.5:0.;p.C[l]=l==0?.02:0.;p.D[l]=0.;p.cv[l]=1.;}
    State liquid{},v{};v[0]=1.;v[NPHASE]=1.;v[2*NPHASE]=-.5;
    Primitive2Conservative(liquid,v,p);State liquid_stencil[5]={liquid,liquid,liquid,liquid,liquid};
    for(int scheme:{3,5,6}){State face;reconstruct(face,liquid_stencil,p,scheme,1);State primitive;Conservative2Primitive(primitive,face,p);
        require(std::abs(primitive[2*NPHASE]+.5)<1e-13,"stretched liquid pressure changed");
        for(int l=1;l<NPHASE;++l)require(face[l]==0. && face[NPHASE+l]==0. && face[2*NPHASE+l]==0.,"absent gas populated in stretched liquid");}
    // Smooth finite-volume cell averages: fifth-order face reconstruction must
    // survive the limiter, rather than passing by falling back to first order.
    Real ratio3=0.,ratio5=0.,ratio5b=0.;p.EOS=0;
    for(int l=0;l<NPHASE;++l){p.A[l]=2.5;p.B[l]=p.C[l]=p.D[l]=0.;p.cv[l]=1.;}
    for(int scheme:{3,5,6}){
      Real previous=0.,ratio=0.;
      for(int n:{32,64,128,256}){
        const Real h=2.*std::acos(-1.)/n,sinc=std::sin(.5*h)/(.5*h);Real error=0.;
        for(int i=0;i<n;++i){State s[5];
            for(int j=0;j<5;++j){State v{};const Real a=.2+.05*sinc*std::sin((i+j-2+.5)*h);v[0]=a;v[1]=.4-a;
                for(int l=2;l<NPHASE;++l)v[l]=.6/(NPHASE-2);
                for(int l=0;l<NPHASE;++l){v[NPHASE+l]=v[l]*(1.+.2*l);v[2*NPHASE+l]=2.*v[l];}
                v[INDEX_MomentumX]=.3;v[INDEX_MomentumY]=-.2;Primitive2Conservative(s[j],v,p);}
            State face;reconstruct(face,s,p,scheme,1);error+=std::abs(face[0]-(.2+.05*std::sin((i+1)*h)))/n;
        }
        if(previous)ratio=previous/error;previous=error;
    }
      std::cerr<<"accuracy scheme="<<scheme<<" ratio="<<ratio<<"\n";require(ratio>(scheme==3?6.:20.),"WENO lost high-order accuracy in smooth flow");
      if(scheme==3)ratio3=ratio;else if(scheme==5)ratio5=ratio;else ratio5b=ratio;
    }
    std::cout<<std::setprecision(17)<<"RECONSTRUCTION {\"faces\":"<<faces<<",\"max_pressure_error\":"<<worstP<<",\"max_partition_error\":"<<worstA<<",\"smooth_refinement_ratio\":"<<ratio5<<",\"smooth_ratio_weno3b\":"<<ratio3<<",\"smooth_ratio_weno5b\":"<<ratio5b<<"}\n";
}amrex::Finalize();}
