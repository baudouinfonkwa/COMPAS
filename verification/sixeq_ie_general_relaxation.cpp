// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Native checks of production generalized instantaneous pressure relaxation.
#include <COMPAS.H>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using State=amrex::Array<amrex::Real,NSTATE>;
void require(bool ok,char const* message){if(!ok)throw std::runtime_error(message);}
int main(int argc,char** argv){
 amrex::Initialize(argc,argv);
 std::string invalid;amrex::ParmParse().query("invalid",invalid);
 int count=0,resets=0,halved_cases=0,direct_cases=0,max_halvings=0,max_newton_steps=0;double max_energy=0.,max_volume=0.,max_work=0.;
 for(int eos=0;eos<4;++eos)for(int trace_case=0;trace_case<4;++trace_case)
 for(int sample=0;sample<20;++sample){
  Parm p{};p.EOS=eos;
  State primitive{},u{},out{},q{};
  const double trace=trace_case==0?.2:(trace_case==1?1.e-30:(trace_case==2?0.:1.e-100));
  for(int l=0;l<NPHASE;++l){
   p.A[l]=1./(.4+.1*l);p.B[l]=eos==0?0.:(eos==2?-(.2+.1*l):(2.+l));
   p.C[l]=eos==3?.04*(l+1):0.;p.D[l]=eos>=2?-.2+.3*l:0.;p.cv[l]=1.;
   primitive[l]=l==0?.3:(l==1?.7-trace:trace);
   primitive[NPHASE+l]=primitive[l]*(1.+.2*l);
   primitive[2*NPHASE+l]=primitive[l]*(1.+.09*sample+.7*l);
  }
  primitive[3*NPHASE]=.3;primitive[3*NPHASE+1]=-.2;
  Primitive2Conservative(u,primitive,p);
  // Mixture energy is independently transported; force a nonzero reset.
  u[INDEX_Total_Energy]+=.13;
  if(!invalid.empty()){
   if(invalid=="energy")u[INDEX_Total_Energy]=-100.;
   if(invalid=="temperature")u[2*NPHASE]=-100.;
   if(invalid=="fraction")u[0]=-.1;
   if(invalid=="covolume"){p.EOS=3;p.C[0]=2.;}
  }
  State before=u;
  SixEqIERelaxationResult info;PressureRelaxation_Cell(out,u,p,-1,-1,-1,&info);
  require(invalid.empty(),"invalid input accepted");
  if(info.initial_halvings>0)++halved_cases;else ++direct_cases;
  max_halvings=std::max(max_halvings,info.initial_halvings);
  max_newton_steps=std::max(max_newton_steps,info.iterations);
  if(info.iterations>0)require(info.initial_residual>0.,"Newton did not start with positive residual");
  Conservative2Primitive(q,out,p);
  double rho=0.,sumH=0.,sumA=0.;
  for(int l=0;l<NPHASE;++l){
   rho+=out[NPHASE+l];sumH+=out[2*NPHASE+l];sumA+=out[l];
   require(out[l]>=0.,"negative relaxed volume");
   if(out[l]>0.)require(std::abs(q[2*NPHASE+l]/out[l]-info.reset_pressure)<1.e-11,"noncommon reset pressure");
   require(out[NPHASE+l]==before[NPHASE+l],"mass changed");
   {
    // Independent EOS-based endpoint-work identity for the provisional root.
    double H;EquationOfState_Phase(H,out[l]*info.work_pressure,out[l],out[NPHASE+l],p,l,true,true);
    double defect=H-u[2*NPHASE+l]+info.work_pressure*(out[l]-u[l]);
    max_work=std::max(max_work,std::abs(defect));
    require(std::abs(defect)<2.e-11,"endpoint work identity failed");
   }
  }
  double kinetic=0.;for(int d=0;d<AMREX_SPACEDIM;++d){
   require(out[3*NPHASE+d]==before[3*NPHASE+d],"momentum changed");
   kinetic+=.5*out[3*NPHASE+d]*out[3*NPHASE+d]/rho;
  }
  require(out[INDEX_Total_Energy]==before[INDEX_Total_Energy],"total energy changed");
  const double error=std::abs(sumH+kinetic-out[INDEX_Total_Energy])/std::max(1.,std::abs(out[INDEX_Total_Energy]));
  max_energy=std::max(max_energy,error);max_volume=std::max(max_volume,std::abs(sumA-1.));
  require(error<3.e-14,"mixture energy not closed");require(std::abs(sumA-1.)<6.e-14,"volume not closed");
  if(std::abs(info.work_pressure-info.reset_pressure)>1.e-5)++resets;
  // An equilibrium state must be a fixed point.
  State twice;PressureRelaxation_Cell(twice,out,p);
  for(int l=0;l<NSTATE;++l)require(std::abs(twice[l]-out[l])<2.e-12*std::max(1.,std::abs(out[l])),"relaxation not idempotent");
  ++count;
 }
 // Equilibrium fixed points at negative admissible pressure, near an EOS pole,
 // almost packed NASG volume, and pressure scales spanning 240 orders.
 int edge_cases=0;
 for(double scale : {1.e-120,1.,1.e120})for(int packed=0;packed<2;++packed)
 for(int near_pole=0;near_pole<2;++near_pole){
  Parm p{};p.EOS=packed?3:1;State v{},u{},out{},q{};
  const double P=scale*(near_pole?(-1.+1.e-10):-.5);
  for(int l=0;l<NPHASE;++l){
   p.A[l]=2.+l;p.B[l]=(p.A[l]+1.)*scale*(1.+l);
   p.C[l]=packed?.999/(1.+l):0.;p.D[l]=0.;p.cv[l]=1.;
   v[l]=l==0?.3:(l==1?.7-1.e-30:1.e-30);
   v[NPHASE+l]=v[l]*(1.+l);v[2*NPHASE+l]=v[l]*P;
  }
  Primitive2Conservative(u,v,p);SixEqIERelaxationResult info;
  PressureRelaxation_Cell(out,u,p,-1,-1,-1,&info);
  Conservative2Primitive(q,out,p);
  for(int l=0;l<NPHASE;++l){
   require(std::abs(out[l]-u[l])<3.e-12,"edge equilibrium volume changed");
   require(std::abs(q[2*NPHASE+l]/out[l]-P)<1.e-11*scale,"edge equilibrium pressure changed");
  }
  ++edge_cases;
 }
 // Signed trace perturbations of mixture-roundoff size must survive reflux
 // without freezing any positive coefficient or changing conserved quantities.
 int roundoff_cases=0;
 for(int eos=0;eos<4;++eos)for(int sign_case=0;sign_case<2;++sign_case){
  Parm p{};p.EOS=eos;State v{},u{},out{};
  for(int l=0;l<NPHASE;++l){
   p.A[l]=2.+l;p.B[l]=eos?3.:0.;p.C[l]=eos==3?.01:0.;p.D[l]=eos>=2?.2:0.;p.cv[l]=1.;
   v[l]=l==0?.3:(l==1?.7:0.);v[NPHASE+l]=v[l];v[2*NPHASE+l]=v[l];
  }
  Primitive2Conservative(u,v,p);
  u[2]=sign_case?-1.e-35:1.e-50;u[NPHASE+2]=sign_case?1.e-35:-1.e-35;
  u[2*NPHASE+2]=-1.e-35;
  PressureRelaxation_Cell(out,u,p);
  double sum=0.;for(int l=0;l<NPHASE;++l){
   sum+=out[2*NPHASE+l];require(out[NPHASE+l]==u[NPHASE+l],"roundoff mass changed");
  }
  require(std::abs(sum-out[INDEX_Total_Energy])<2.e-13,"roundoff reset did not close energy");
  ++roundoff_cases;
 }
 // Valid disequilibrium phases can have a mixture pressure outside the
 // common EOS branch, including either sign of the pressure bound.
 int outside_branch_cases=0;
 for(int positive_bound=0;positive_bound<2;++positive_bound){
  Parm p{};p.EOS=positive_bound?2:1;State v{},u{},out{};
  double mixture=0.;
  for(int l=0;l<NPHASE;++l){
   p.A[l]=2.;p.B[l]=l==0?3.*(positive_bound?-1.:10.):0.;p.C[l]=p.D[l]=0.;p.cv[l]=1.;
   v[l]=l==0?(positive_bound?.1:.9):(l==1?(positive_bound?.9:.1):0.);
   v[NPHASE+l]=v[l];
   const double pressure=l==0?(positive_bound?1.1:-9.):(positive_bound?.1:1.);
   v[2*NPHASE+l]=v[l]*pressure;mixture+=v[2*NPHASE+l];
  }
  Primitive2Conservative(u,v,p);SixEqIERelaxationResult info;
  PressureRelaxation_Cell(out,u,p,-1,-1,-1,&info);
  const double bound=positive_bound?1.:0.;
  require(mixture<=bound,"test did not place mixture outside common branch");
  require(info.initial_work_pressure>bound && info.initial_residual>0.,"branch initialization failed");
  require(info.work_pressure>bound && info.initial_halvings>0,"branch root/search failed");
  double total=0.;for(int l=0;l<NPHASE;++l)total+=out[l];
  require(std::abs(total-1.)<3.e-14,"outside-branch volumes not closed");
  ++outside_branch_cases;
 }
 // Two ideal-gas mixtures have an independent closed-form root. Reversing
 // which gamma carries the high pressure puts p0 on opposite sides of it.
 for(int direct=0;direct<2;++direct){
  Parm p{};p.EOS=0;State v{},u{},out{};double numerator=0.,denominator=0.;
  for(int l=0;l<NPHASE;++l){
   const double gamma=1.4+.2*l;
   p.A[l]=1./(gamma-1.);p.B[l]=p.C[l]=p.D[l]=0.;p.cv[l]=1.;
   v[l]=l==0?.3:(l==1?.7:0.);v[NPHASE+l]=v[l];
   const double pressure=(l==(direct?0:1))?3.:1.;v[2*NPHASE+l]=v[l]*pressure;
   numerator+=v[l]*pressure/gamma;denominator+=v[l]/gamma;
  }
  Primitive2Conservative(u,v,p);SixEqIERelaxationResult info;
  PressureRelaxation_Cell(out,u,p,-1,-1,-1,&info);
  require(std::abs(info.work_pressure-numerator/denominator)<2.e-13,"ideal-mixture analytic root mismatch");
  require(info.initial_residual>0.,"targeted Newton start was not positive");
  require(direct?info.initial_halvings==0:info.initial_halvings>0,"targeted initial search path mismatch");
  if(info.initial_halvings>0)++halved_cases;else ++direct_cases;
 }
 require(halved_cases>0 && direct_cases>0,"both initial-guess paths were not exercised");
 require(resets>0,"unconditional reset not exercised");
 std::cout<<std::setprecision(17)<<"RELAXATION_RESULTS {\"cases\":"<<count<<",\"halved_start_cases\":"<<halved_cases<<",\"direct_start_cases\":"<<direct_cases
 <<",\"max_initial_halvings\":"<<max_halvings<<",\"max_newton_steps\":"<<max_newton_steps
 <<",\"outside_branch_cases\":"<<outside_branch_cases<<",\"edge_cases\":"<<edge_cases<<",\"roundoff_cases\":"<<roundoff_cases<<",\"distinct_unconditional_resets\":"<<resets
 <<",\"max_energy_residual\":"<<max_energy<<",\"max_volume_residual\":"<<max_volume<<",\"max_endpoint_work_residual\":"<<max_work<<"}\n";
 amrex::Finalize();
}
