// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Independent jump/energy/EOS checks of the opt-in production kernels.
#include <COMPAS.H>
#include <AMReX_FArrayBox.H>
#include <random>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using Real=amrex::Real;
using State=amrex::Array<Real,NSTATE>;
void require(bool c,char const* text){if(!c)throw std::runtime_error(text);}
void near(Real a,Real b,char const* text,Real tolerance=3.e-11){require(std::abs(a-b)<=tolerance*std::max({Real(1.),std::abs(a),std::abs(b)}),text);}
State state(Parm const& p,Real a,Real u,Real v,Real T)
{
    State q{};Real rho=0.,H=0.;
    for(int l=0;l<NPHASE;++l){q[l]=l==0?a:l==1?1.-a:0.;q[NPHASE+l]=q[l]*(l+1.);
        EquationOfState_Phase(q[2*NPHASE+l],q[NPHASE+l]*T,q[l],q[NPHASE+l],p,l,false,true);
        rho+=q[NPHASE+l];H+=q[2*NPHASE+l];}
    q[INDEX_MomentumX]=rho*u;q[INDEX_MomentumY]=rho*v;q[INDEX_Total_Energy]=H+.5*rho*(u*u+v*v);return q;
}
SixEqIEDiffusionFace face(Parm const& p,State const& l,State const& r,State& flux)
{
    amrex::Array<State,FaceStencil::N> stencil;
    for(int s=0;s<FaceStencil::N;++s){int i,j,k;FaceStencil::OffsetOnFace(s,0,i,j,k);stencil[s]=i<0?l:r;}
    SixEqIEDiffusionFace data;Real d;
    FDM_Conservative2FluxDiffusion_D(flux,d,stencil,1.,1.,p,0,&data);return data;
}
int main(int argc,char**argv)
{
    amrex::Initialize(argc,argv);
    {
        Parm p{};p.DynamicInit();p.Initialize();p.Physics_Parm.compatible_diffusion=1;
        for(int l=0;l<NPHASE;++l){p.A[l]=2.+l;p.B[l]=0.;p.C[l]=0.;p.D[l]=0.;p.cv[l]=2.+l;p.mu[l]=l==0?.7:.2;p.muB[l]=l==0?.1:.4;p.kappa[l]=l==0?.6:.9;}
        State flux{};auto data=face(p,state(p,1.,0.,-.5/.7,2.),state(p,0.,0.,.5/.2,2.),flux);
        near(data.traction[1],1.,"layered shear traction");near(data.work,0.,"layered shear work");
        Real etaL=4.*.7/3.+.1,etaR=4.*.2/3.+.4;
        data=face(p,state(p,1.,-.5/etaL,0.,2.),state(p,0.,.5/etaR,0.,2.),flux);
        near(data.traction[0],1.,"layered normal traction");near(data.work,0.,"layered normal work");
        data=face(p,state(p,1.,0.,0.,2.),state(p,0.,0.,0.,4.),flux);
        near(data.heat,2.*(.6*.9/(.6+.9))*2.,"sharp material heat flux");
        for(int l=0;l<NPHASE;++l)near(flux[2*NPHASE+l],0.,"heat flux directly populated a massless phase");
        // Uniform rigid rotation must have zero viscous traction despite a mu jump.
        for(int dir=0;dir<2;++dir){
            amrex::Array<State,FaceStencil::N> s;
            for(int q=0;q<FaceStencil::N;++q){int i,j,k;FaceStencil::OffsetOnFace(q,dir,i,j,k);
                Real x=i+(dir==0?.5:0.),y=j+(dir==1?.5:0.);s[q]=state(p,(dir==0?x:y)<0.?1.:0.,-y,x,2.);}
            Real bound;FDM_Conservative2FluxDiffusion_D(flux,bound,s,1.,1.,p,dir,&data);
            near(data.traction[0],0.,"rigid rotation x stress");near(data.traction[1],0.,"rigid rotation y stress");
        }
        std::mt19937 generator(731);std::uniform_real_distribution<Real> rand(0.,1.);
        Real worstP=0.,worstT=0.,worstEnergy=0.;int projections=0;
        for(int eos=0;eos<4;++eos)for(int test=0;test<200;++test){
            p.EOS=eos;
            Real P=std::pow(10.,-6.+12.*rand(generator)),T=std::pow(10.,-3.+6.*rand(generator));
            State q{};Real thermal=0.,rho=0.;
            amrex::Array<Real,NPHASE> alphas{};
            Real alpha=.05+.9*rand(generator);alphas[0]=alpha;alphas[1]=1.-alpha;
            for(int l=0;l<NPHASE;++l){
                p.A[l]=.5+5.*rand(generator);p.cv[l]=.5+5.*rand(generator);
                Real pi=eos==0?0.:P*(eos==2?(-.5+2.*rand(generator)):1000.*rand(generator));
                p.B[l]=(1.+p.A[l])*pi;
                Real rhoT=p.A[l]*(P+pi)/p.cv[l];
                p.C[l]=eos==3?.05*T/rhoT:0.;p.D[l]=eos>=2?(.5-rand(generator))*T*p.cv[l]:0.;
                Real density=rhoT/(T+p.C[l]*rhoT);
                q[l]=alphas[l];q[NPHASE+l]=alphas[l]*density;
                EquationOfState_Phase(q[2*NPHASE+l],q[l]*P,q[l],q[NPHASE+l],p,l,true,true);
                thermal+=q[2*NPHASE+l];rho+=q[NPHASE+l];
            }
            q[INDEX_MomentumX]=rho*.3;q[INDEX_MomentumY]=-rho*.2;
            q[INDEX_Total_Energy]=thermal+.5*rho*(.3*.3+.2*.2);
            // Nonequilibrium alpha/H predictor, conservative invariants fixed.
            q[0]=alpha*.8;q[1]=1.-q[0];
            q[2*NPHASE]*=1.2;
            State relaxed;SixEqIE_RelaxPressureTemperatureCell(relaxed,q,p);
            amrex::Array<Real,NPHASE> temperatures;amrex::Array<bool,NPHASE> present;
            Conservative2PhaseTemperature(temperatures,present,relaxed,p);
            State primitive;Conservative2Primitive(primitive,relaxed,p);
            Real sumAlpha=0.,sumH=0.;
            for(int l=0;l<NPHASE;++l){
                sumAlpha+=relaxed[l];sumH+=relaxed[2*NPHASE+l];
                require(relaxed[NPHASE+l]==q[NPHASE+l],"P-T changed partial mass");
                if(alphas[l]>0.){
                    Real pressure=primitive[2*NPHASE+l]/relaxed[l];
                    worstP=std::max(worstP,std::abs(pressure-P)/P);
                    worstT=std::max(worstT,std::abs(temperatures[l]-T)/T);
                    near(relaxed[l],alphas[l],"P-T failed known equilibrium fractions",3.e-10);
                }else require(relaxed[l]==0. && relaxed[2*NPHASE+l]==0.,"P-T populated absent phase");
            }
            near(sumAlpha,1.,"P-T volume sum");
            near(sumH,thermal,"P-T energy closure",3.e-10);
            worstEnergy=std::max(worstEnergy,std::abs(sumH-thermal)/std::max(Real(1.),std::abs(thermal)));
            require(relaxed[INDEX_Total_Energy]==q[INDEX_Total_Energy],"P-T changed total energy");++projections;
        }
        require(worstP<5.e-9 && worstT<5.e-10,"P-T inaccurate pressure/temperature");
        // Real production NC source versus real face fluxes, periodic shear.
        p.A.fill(2.5);p.B.fill(0.);p.C.fill(0.);p.D.fill(0.);p.cv.fill(2.);
        p.mu.fill(.4);p.muB.fill(0.);p.kappa.fill(0.);
        Real maxCompatibility=0.,checkerboardHeating=0.;
        for(int mode=0;mode<2;++mode){
            int N=16;Real h=1./N;
            amrex::Box box(amrex::IntVect(AMREX_D_DECL(-1,-1,0)),amrex::IntVect(AMREX_D_DECL(N,1,0)));
            amrex::FArrayBox states(box,NSTATE),rhs(box,NSTATE),ncx(box,NC_TERMS),ncy(box,NC_TERMS),fx(box,NSTATE),fy(box,NSTATE);
            auto u=states.array();auto ax=ncx.array();auto ay=ncy.array();auto f=fx.array();
            auto sample=[&](int i){Real v=mode==1?(i%2==0?1.:-1.):std::sin(2.*std::acos(-1.)*(i+.5)*h);return state(p,.4,0.,v,2.);};
            for(int i=-1;i<=N;++i)for(int j=-1;j<=1;++j){auto q=sample(i);for(int c=0;c<NSTATE;++c)u(i,j,0,c)=q[c];}
            rhs.setVal(0.);ncx.setVal(0.);ncy.setVal(0.);fx.setVal(0.);fy.setVal(0.);
            for(int i=0;i<=N;++i){
                amrex::Array<State,FaceStencil::N> stencil;
                for(int q=0;q<FaceStencil::N;++q){int di,dj,dk;FaceStencil::OffsetOnFace(q,0,di,dj,dk);stencil[q]=sample(i+di);}
                Real bound;FDM_Conservative2FluxDiffusion_D(flux,bound,stencil,h,h,p,0,&data);
                for(int c=0;c<NSTATE;++c)f(i,0,0,c)=flux[c];
                for(int b=0;b<2;++b)ax(i,0,0,IE_NC_TRACTION+b)=data.traction[b];
                ax(i,0,0,IE_NC_WORK)=data.work;ax(i,0,0,IE_NC_HEAT)=data.heat;
            }
            Real integrated=0.;
            for(int i=0;i<N;++i){
                FVM_SurfaceIntegral_NC(i,0,0,rhs.array(),states.const_array(),fx.const_array(),fy.const_array(),ncx.const_array(),ncy.const_array(),h,h,p,p.Physics_Parm);
                Real heat=0.;for(int l=0;l<NPHASE;++l)heat+=rhs.array()(i,0,0,2*NPHASE+l);
                Real v=sample(i)[INDEX_MomentumY]/1.6;
                Real expected=(f(i+1,0,0,INDEX_Total_Energy)-f(i,0,0,INDEX_Total_Energy)
                            -v*(f(i+1,0,0,INDEX_MomentumY)-f(i,0,0,INDEX_MomentumY)))/h;
                maxCompatibility=std::max(maxCompatibility,std::abs(heat-expected)/std::max(Real(1.),std::abs(expected)));
                integrated+=heat*h;
            }
            if(mode==1){checkerboardHeating=integrated;near(integrated,409.6,"checkerboard lost viscous heating");}
        }
        require(maxCompatibility<2.e-13,"diffusion source not compatible with face work/momentum");
        // Full 2D variable-mixture tensor stress: global dissipation check.
        // Local compatible heating is recorded, not assumed positive.
        Real minimumGlobalDissipation=std::numeric_limits<Real>::max(),minimumLocalHeating=0.;
        for(int trial=0;trial<50;++trial){
            constexpr int N=8;Real h=1./N;
            for(int l=0;l<NPHASE;++l){p.mu[l]=std::pow(10.,-2.+4.*rand(generator));p.muB[l]=std::pow(10.,-2.+4.*rand(generator));}
            amrex::Array<State,N*N> grid;
            for(int i=0;i<N*N;++i)grid[i]=state(p,.05+.9*rand(generator),2.*rand(generator)-1.,2.*rand(generator)-1.,2.);
            auto sample=[&](int i,int j)->State const&{return grid[((j+N)%N)*N+(i+N)%N];};
            amrex::Array<amrex::Array<SixEqIEDiffusionFace,N*N>,2> faces;
            for(int dir=0;dir<2;++dir)for(int j=0;j<N;++j)for(int i=0;i<N;++i){
                amrex::Array<State,FaceStencil::N> stencil;
                for(int q=0;q<FaceStencil::N;++q){int di,dj,dk;FaceStencil::OffsetOnFace(q,dir,di,dj,dk);stencil[q]=sample(i+di,j+dj);}
                Real bound;FDM_Conservative2FluxDiffusion_D(flux,bound,stencil,h,h,p,dir,&faces[dir][j*N+i]);
            }
            Real loss=0.;
            for(int j=0;j<N;++j)for(int i=0;i<N;++i){
                State const& cell=sample(i,j);Real rho=cell[NPHASE]+cell[NPHASE+1];
                Real Q=0.;
                for(int dir=0;dir<2;++dir){
                    auto const& left=faces[dir][j*N+i];
                    auto const& right=faces[dir][((j+(dir==1))%N)*N+(i+(dir==0))%N];
                    Q+=(right.work-left.work)/h;
                    for(int b=0;b<2;++b){
                        Real power=cell[INDEX_MomentumX+b]/rho*(right.traction[b]-left.traction[b])/h;
                        loss-=power*h*h;Q-=power;
                    }
                }
                minimumLocalHeating=std::min(minimumLocalHeating,Q);
            }
            minimumGlobalDissipation=std::min(minimumGlobalDissipation,loss);
            require(loss>=-1.e-10,"variable-viscosity tensor operator creates global kinetic energy");
        }
        std::cout<<std::setprecision(17)<<"COMPATIBLE_DIFFUSION {\"projections\":"<<projections
            <<",\"max_pressure_error\":"<<worstP<<",\"max_temperature_error\":"<<worstT
            <<",\"max_energy_error\":"<<worstEnergy<<",\"max_source_error\":"<<maxCompatibility
            <<",\"checkerboard_heating\":"<<checkerboardHeating
            <<",\"variable_tensor_trials\":50,\"minimum_global_dissipation\":"<<minimumGlobalDissipation
            <<",\"minimum_local_heating\":"<<minimumLocalHeating<<"}\n";
    }
    amrex::Finalize();
}
