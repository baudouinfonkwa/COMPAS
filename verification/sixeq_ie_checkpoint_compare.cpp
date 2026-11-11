// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Compare valid saved cells, excluding checkpoint ghost cells/uninitialized padding.
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BoxIterator.H>
#include <iostream>
#include <iomanip>
#include <cstring>
int main(int argc,char**argv){
 amrex::Initialize(argc,argv);
 {
  amrex::ParmParse pp;std::string left,right;int levels=1;
  pp.get("left",left);pp.get("right",right);pp.query("levels",levels);
  long changed=0,values=0;double error=0.;
  for(int lev=0;lev<levels;++lev){
   amrex::MultiFab a,b;amrex::VisMF::Read(a,left+"/Level_"+std::to_string(lev)+"/");
   amrex::VisMF::Read(b,right+"/Level_"+std::to_string(lev)+"/");
   if(a.boxArray()!=b.boxArray() || a.nComp()!=b.nComp())amrex::Abort("comparison layouts differ");
   for(amrex::MFIter mfi(a);mfi.isValid();++mfi){auto aa=a.const_array(mfi),bb=b.const_array(mfi);
    for(amrex::BoxIterator it(mfi.validbox());it.ok();++it)for(int n=0;n<a.nComp();++n){
     double x=aa(it(),n),y=bb(it(),n);++values;if(std::memcmp(&x,&y,sizeof(double)))++changed;
     if(!std::isfinite(x) || !std::isfinite(y))amrex::Abort("nonfinite comparison state");
     error=std::max(error,std::abs(x-y)/std::max({1.,std::abs(x),std::abs(y)}));
    }
   }
  }
  std::cout<<std::setprecision(17)<<"CHECKPOINT_COMPARE {\"changed_values\":"<<changed<<",\"values\":"<<values<<",\"max_relative_error\":"<<error<<"}\n";
 }
 amrex::Finalize();
}
