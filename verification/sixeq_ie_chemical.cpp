// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Native verification of equilibrium chemistry; synthetic EOS calibrations
// isolate conservation and active-set algebra from empirical material fitting.
#include <COMPAS.H>
#include <random>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using Real=amrex::Real;
using State=amrex::Array<Real,NSTATE>;
void require(bool condition,char const* text){if(!condition)throw std::runtime_error(text);}
void near(Real a,Real b,char const* text,Real tol=2.e-8){require(std::abs(a-b)<=tol*std::max({Real(1.),std::abs(a),std::abs(b)}),text);}
Parm parameters()
{
    Parm p{};
    p.Physics_Parm.pressure_relaxation=1;p.Physics_Parm.relaxation=Relaxation::PressureTemperature;
    p.Physics_Parm.chemical_relaxation=1;
    p.Physics_Parm.chemical_Tref=1.;p.Physics_Parm.chemical_pref=1.;
    for(int l=0;l<NPHASE;++l){p.A[l]=2.;p.B[l]=0.;p.C[l]=0.;p.D[l]=0.;p.cv[l]=1.;
        p.Physics_Parm.chemical_group[l]=l;p.Physics_Parm.chemical_sref[l]=0.;}
    p.A[0]=2.;p.cv[0]=2.;p.A[1]=1.;p.cv[1]=1.;p.D[1]=3.;
    p.Physics_Parm.chemical_sref[1]=2.;p.Physics_Parm.chemical_group[1]=0;
    return p;
}
State state(Parm const& p,amrex::Array<Real,NPHASE> const& mass,Real H)
{
    State q{},out{};Real rho=0.;
    for(int l=0;l<NPHASE;++l){q[l]=1./NPHASE;q[INDEX_Mass1+l]=mass[l];rho+=mass[l];}
    q[INDEX_MomentumX]=.2*rho;q[INDEX_MomentumY]=-.1*rho;
    q[INDEX_Total_Energy]=H+.025*rho;
    SixEqIE_RelaxPressureTemperatureCell(out,q,p);return out;
}
Real worstEnergy=0.,worstAffinity=0.,worstTemperature=0.;int checks=0;
void audit(State const& before,State const& after,Parm const& p)
{
    near(before[INDEX_Total_Energy],after[INDEX_Total_Energy],"total energy changed",0.);
    for(int d=0;d<AMREX_SPACEDIM;++d)near(before[INDEX_MomentumX+d],after[INDEX_MomentumX+d],"momentum changed",0.);
    Real rho=0.,H=0.,sumAlpha=0.,pressure=0.,T=0.;
    for(int l=0;l<NPHASE;++l){
        require(after[l]>=0. && after[INDEX_Mass1+l]>=0.,"negative phase");
        rho+=after[INDEX_Mass1+l];H+=after[INDEX_Energy1+l];sumAlpha+=after[l];
        Real ap;EquationOfState_Phase(ap,after[INDEX_Energy1+l],after[l],after[INDEX_Mass1+l],p,l,true,false);pressure+=ap;
        if(after[INDEX_Mass1+l]==0.){require(after[l]==0. && after[INDEX_Energy1+l]==0.,"absent phase has residual volume or energy");continue;}
        Real mt;EquationOfState_Phase(mt,after[INDEX_Energy1+l],after[l],after[INDEX_Mass1+l],p,l,false,false);
        Real tl=mt/after[INDEX_Mass1+l];if(T>0.)worstTemperature=std::max(worstTemperature,std::abs(tl-T)/T);else T=tl;
    }
    Real K=0.;for(int d=0;d<AMREX_SPACEDIM;++d)K+=.5*after[INDEX_MomentumX+d]*after[INDEX_MomentumX+d]/rho;
    worstEnergy=std::max(worstEnergy,std::abs(H+K-after[INDEX_Total_Energy])/std::max(Real(1.),std::abs(after[INDEX_Total_Energy])));
    near(sumAlpha,1.,"volume constraint");
    for(int group=0;group<NPHASE;++group){Real a=0.,b=0.;int members=0;
        for(int l=0;l<NPHASE;++l)if(p.Physics_Parm.chemical_group[l]==group){a+=before[INDEX_Mass1+l];b+=after[INDEX_Mass1+l];++members;}
        near(a,b,"substance mass changed",1.e-13);
        if(members==1)for(int l=0;l<NPHASE;++l)if(p.Physics_Parm.chemical_group[l]==group)
            near(before[INDEX_Mass1+l],after[INDEX_Mass1+l],"inert mass changed",0.);
    }
    Real entropy;Real residual=SixEqIE_ChemicalResidual(after,pressure,T,p,entropy);worstAffinity=std::max(worstAffinity,residual);
    require(residual<=1.01*p.Physics_Parm.chemical_tolerance,"chemical KKT residual");
    Real oldP=0.,oldT=0.;State old;
    if(SixEqIE_TryPressureTemperatureCell(old,before,p,oldP,oldT)==SixEqIE_PTStatus::Success){
        Real oldS;SixEqIE_ChemicalResidual(old,oldP,oldT,p,oldS);
        require(entropy>=oldS-2.e-11*std::max(Real(1.),std::abs(oldS)),"entropy decreased");
    }
    State again;SixEqIE_RelaxChemicalCell(again,after,p);
    for(int n=0;n<NSTATE;++n)near(after[n],again[n],"relaxation not idempotent",1.e-8);
    ++checks;
}
int main(int argc,char**argv)
{
    amrex::Initialize(argc,argv);
    {
        Parm p=parameters();
        amrex::Array<Real,NPHASE> m{};m[0]=.64;m[1]=.16;m[2]=.4;
        for(Real target:{Real(-.2),Real(.4),Real(1.2)}){
            auto q=state(p,m,2.+1.6*target);State out;SixEqIE_RelaxChemicalCell(out,q,p);audit(q,out,p);
            if(target<0.)require(out[INDEX_Mass1+1]==0.,"complete condensation did not select exact endpoint");
            else if(target>1.)require(out[INDEX_Mass1]==0.,"complete evaporation did not select exact endpoint");
            else{near(out[INDEX_Mass1+1],.8*target,"wrong coexistence mass");near(out[0],.8*(1.-target),"wrong coexistence volume");}
        }
        // Initially absent vapor can enter: equilibrium nucleation, no mass floor.
        m[0]=.8;m[1]=0.;auto q=state(p,m,2.64);State out;SixEqIE_RelaxChemicalCell(out,q,p);audit(q,out,p);near(out[INDEX_Mass1+1],.32,"absent phase did not reappear");
        // An invalid incoming phase partition still has valid conserved groups:
        // solve the group flash, without clipping or changing total water/air.
        q[INDEX_Mass1]=.81;q[INDEX_Mass1+1]=-.01;q[1]=-.01;
        SixEqIE_RelaxChemicalCell(out,q,p);audit(q,out,p);near(out[INDEX_Mass1+1],.32,"group flash failed to recover equilibrium partition");
        auto bad=q;bad[INDEX_Mass1+2]=-.4;Real ps,ts;State seed;
        require(!SixEqIE_ChemicalGroupSeed(seed,bad,p,ps,ts),"negative inert mass was repaired");
        bad=q;bad[INDEX_Total_Energy]=-100.;
        require(!SixEqIE_ChemicalGroupSeed(seed,bad,p,ps,ts),"inadmissible group energy was repaired");
        // Feasible interior with neither pure endpoint nor uniform split admissible.
        Parm constrained=parameters();constrained.C[0]=2.;constrained.D[1]=2.;
        State raw{};raw[INDEX_Mass1]=2.;raw[INDEX_Mass1+1]=-1.;raw[INDEX_Total_Energy]=1.5;
        require(SixEqIE_ChemicalGroupSeed(seed,raw,constrained,ps,ts),"affine feasibility seed failed");
        near(seed[INDEX_Mass1]+seed[INDEX_Mass1+1],1.,"seed changed group mass",1.e-13);
        require(seed[INDEX_Total_Energy]==raw[INDEX_Total_Energy],"seed changed energy");
        // Group absent entirely: inert state remains unchanged chemically.
        m.fill(0.);m[2]=2.;q=state(p,m,2.);SixEqIE_RelaxChemicalCell(out,q,p);audit(q,out,p);
        // Nucleation from stretched pure liquid: the absent ideal vapor has
        // no EOS state at the incoming negative pressure; do not certify that
        // endpoint merely because its Gibbs potential cannot be evaluated.
        p=parameters();p.B[0]=9.;m.fill(0.);m[0]=2.8;
        q=state(p,m,8.6);SixEqIE_RelaxChemicalCell(out,q,p);audit(q,out,p);
        require(out[INDEX_Mass1+1]>0.,"vapor did not nucleate from stretched liquid");
        // Disabled dispatch reproduces the P-T projection bit for bit.
        p.Physics_Parm.chemical_relaxation=0;State a,b;SixEqIE_RelaxEquilibriumCell(a,q,p);SixEqIE_RelaxPressureTemperatureCell(b,q,p);
        for(int n=0;n<NSTATE;++n)require(a[n]==b[n],"disabled chemistry changed PT arithmetic");
        std::mt19937 random(901);std::uniform_real_distribution<Real> uniform(0.,1.);
        // Known coexistence equilibria across EOS 0--3 coefficient families.
        for(int sample=0;sample<100;++sample){
            p=parameters();p.EOS=sample%4;
            for(int l=0;l<NPHASE;++l){
                p.A[l]=1.+2.*uniform(random);p.cv[l]=1.+uniform(random);
                p.B[l]=sample%4==0?0.:(p.A[l]+1.)*(l==0?3.:.2);
                p.C[l]=sample%4==3?.03*(l+1):0.;p.D[l]=sample%4>=2?.2*l:0.;
                p.Physics_Parm.chemical_sref[l]=0.;
            }
            // Calibrate g0=g1 at p=T=1 using the explicit reference entropies.
            SixEqIE_ChemicalThermo t0,t1;SixEqIE_ChemicalThermodynamics(t0,1.,1.,0,p);SixEqIE_ChemicalThermodynamics(t1,1.,1.,1,p);
            p.Physics_Parm.chemical_sref[1]=t1.g-t0.g;
            m.fill(0.);m[0]=.25+.2*uniform(random);m[1]=.15+.1*uniform(random);m[2]=.2;
            Real volume=0.,H=0.;for(int l=0;l<NPHASE;++l){SixEqIE_ChemicalThermo t;SixEqIE_ChemicalThermodynamics(t,1.,1.,l,p);volume+=m[l]*t.v;}
            for(int l=0;l<NPHASE;++l)m[l]/=volume;
            for(int l=0;l<NPHASE;++l){SixEqIE_ChemicalThermo t;SixEqIE_ChemicalThermodynamics(t,1.,1.,l,p);H+=m[l]*t.e;}
            const Real target=m[1];m[0]+=.05*m[1];m[1]*=.95;
            q=state(p,m,H);SixEqIE_RelaxChemicalCell(out,q,p);audit(q,out,p);near(out[INDEX_Mass1+1],target,"known EOS coexistence state",1.e-7);
        }
#if NPHASE >= 5
        // Two distinct groups exchange simultaneously through common p,T.
        p=parameters();p.A[2]=2.;p.cv[2]=1.;p.A[3]=2.;p.cv[3]=1.2;p.D[3]=.3;
        p.A[4]=1.;p.cv[4]=1.2;p.D[4]=2.2;p.Physics_Parm.chemical_sref[4]=2.5;
        p.Physics_Parm.chemical_group[4]=3;
        m={.32,.08,.36,.35,.15};q=state(p,m,2.53);SixEqIE_RelaxChemicalCell(out,q,p);audit(q,out,p);
        near(out[INDEX_Mass1+1],.12,"coupled first group",1.e-7);near(out[INDEX_Mass1+4],.2,"coupled second group",1.e-7);
        // Connected group with absent intermediate must permit direct 0<->2.
        p=parameters();p.A[2]=1.;p.cv[2]=1.;p.D[2]=3.;p.Physics_Parm.chemical_sref[2]=2.;
        p.D[1]=20.;p.Physics_Parm.chemical_sref[1]=0.;p.Physics_Parm.chemical_group[2]=0;
        m={.64,0.,.16,.4,0.};q=state(p,m,2.64);SixEqIE_RelaxChemicalCell(out,q,p);audit(q,out,p);
        require(out[INDEX_Mass1+1]==0.,"unstable intermediate phase was populated");near(out[INDEX_Mass1+2],.32,"absent intermediate blocked group transfer");
#endif
        require(worstEnergy<1.e-11 && worstTemperature<1.e-10,"EOS closure errors");
        std::cout<<std::setprecision(17)<<"CHEMICAL_RESULTS {\"checks\":"<<checks<<",\"max_energy_error\":"<<worstEnergy
            <<",\"max_affinity\":"<<worstAffinity<<",\"max_temperature_error\":"<<worstTemperature<<"}\n";
    }
    amrex::Finalize();
}
