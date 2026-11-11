// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Actual class/flux-register check of the extra compatible diffusion NC slots.
#include <COMPAS.H>
#include <Class_Compressible6Eq_IE_NPhase_PhaseField.H>
#include <AMReX_BoxIterator.H>
#include <iostream>
#include <iomanip>
#include <stdexcept>
int main(int argc,char**argv){
 amrex::Initialize(argc,argv);
 {
  Compressible6Eq_IE_NPhase_PhaseField sim;
  auto domain=sim.Geom(0).Domain();amrex::BoxArray coarse(domain);coarse.maxSize(4);
  amrex::DistributionMapping cdm(coarse);auto half=domain;half.setBig(0,3);half.refine(sim.refRatio(0));
  amrex::BoxArray fine(half);fine.maxSize(4);amrex::DistributionMapping fdm(fine);
  sim.SetBoxArray(0,coarse);sim.SetDistributionMap(0,cdm);sim.SetBoxArray(1,fine);sim.SetDistributionMap(1,fdm);
  Parm& p=*sim.h_parm;
  amrex::Array<amrex::Real,NSTATE> q{},v{};
  for(int l=0;l<NPHASE;++l){v[l]=l==0?.4:.6;v[NPHASE+l]=v[l]*(l+1.);v[2*NPHASE+l]=v[l];}
  v[INDEX_MomentumX]=.3;v[INDEX_MomentumY]=-.2;Primitive2Conservative(q,v,p);
  amrex::MultiFab before(coarse,cdm,NSTATE,0),after(coarse,cdm,NSTATE,0);
  for(int n=0;n<NSTATE;++n){before.setVal(q[n],n,1);after.setVal(q[n],n,1);}
  amrex::Vector<std::unique_ptr<amrex::FluxRegister>> nc(2);
  nc[1]=std::make_unique<amrex::FluxRegister>(fine,fdm,sim.refRatio(0),1,NC_TERMS);nc[1]->setVal(0.);
  amrex::FluxRegister conservative(fine,fdm,sim.refRatio(0),1,NSTATE);conservative.setVal(0.);
  auto face_boxes=amrex::convert(coarse,amrex::IntVect::TheDimensionVector(0));
  amrex::MultiFab fn(face_boxes,cdm,NC_TERMS,0),fc(face_boxes,cdm,NSTATE,0);fn.setVal(0.);fc.setVal(0.);
  fn.setVal(.2,IE_NC_TRACTION,1);fn.setVal(-.1,IE_NC_TRACTION+1,1);
  fn.setVal(.3,IE_NC_WORK,1);fn.setVal(.7,IE_NC_HEAT,1);
  fc.setVal(-.2,INDEX_MomentumX,1);fc.setVal(.1,INDEX_MomentumY,1);fc.setVal(-1.,INDEX_Total_Energy,1);
  nc[1]->CrseInit(fn,0,0,0,NC_TERMS,.001,amrex::FluxRegister::COPY);
  conservative.CrseInit(fc,0,0,0,NSTATE,.001,amrex::FluxRegister::COPY);
  conservative.Reflux(after,1.,0,0,NSTATE,sim.Geom(0));
  sim.RefluxSourceTerms(after,before,0,nc);
  double error=0.;int changed=0;
  for(amrex::MFIter mfi(after);mfi.isValid();++mfi){auto a=after.const_array(mfi),b=before.const_array(mfi);
   for(amrex::BoxIterator it(mfi.validbox());it.ok();++it){auto iv=it();
    double H=0.,expected=a(iv,INDEX_Total_Energy)-b(iv,INDEX_Total_Energy);
    expected-=.3*(a(iv,INDEX_MomentumX)-b(iv,INDEX_MomentumX));
    expected-=-.2*(a(iv,INDEX_MomentumY)-b(iv,INDEX_MomentumY));
    for(int l=0;l<NPHASE;++l)H+=a(iv,2*NPHASE+l)-b(iv,2*NPHASE+l);
    if(expected!=0.)++changed;
    error=std::max(error,std::abs(H-expected));
   }
  }
  if(!changed || error>2.e-13)throw std::runtime_error("AMR compatible diffusion source has wrong sign or slots");
  std::cout<<std::setprecision(17)<<"COMPATIBLE_REFLUX {\"changed_cells\":"<<changed<<",\"max_source_error\":"<<error<<"}\n";
 }
 amrex::Finalize();
}
