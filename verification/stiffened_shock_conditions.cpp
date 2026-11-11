// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

#include <PhysicsFunctions.H>
#include <iostream>
#include <stdexcept>
#include <algorithm>
void near(amrex::Real a,amrex::Real b,char const* message,amrex::Real tolerance=2.e-11){
    if(std::abs(a-b)>tolerance*std::max({amrex::Real(1.),std::abs(a),std::abs(b)}))throw std::runtime_error(message);
}
int main(int argc,char** argv){
    amrex::Initialize(argc,argv);
    int cases=0,legacy=0;
    for(int eos=0;eos<4;++eos)for(double A:{2.5,1./3.4})for(double mach:{1.,1.01,1.1,2.,10.})
    for(double density:{.1,1.,1000.})for(double velocity:{-10.,0.,7.}){
        const double B=eos==0?0.:(eos==2?-1.:20.),b=eos==3?.0005:0.,D=eos>=2?-123.:0.,p0=2.;
        double out[5];StiffenedShockConditions(mach,velocity,p0,density,A,B,b,D,out);
        double r=out[0],u=out[1],p=out[2],speed=out[4];
        double e0=(A*p0+B)*(1./density-b)+D,e1=(A*p+B)*(1./r-b)+D;
        near(density*(speed-velocity),r*(speed-u),"mass jump");
        near(p0+density*(speed-velocity)*(speed-velocity),p+r*(speed-u)*(speed-u),"momentum jump");
        near(e0+p0/density+.5*(speed-velocity)*(speed-velocity),e1+p/r+.5*(speed-u)*(speed-u),"energy jump");
        near(out[3],r*(e1+.5*u*u),"EOS total energy");
        if(eos<2){
            const double g=1.+1./A,pi=B/(A+1.),c=std::sqrt(g*(p0+pi)/density);
            double old_p=p0*(1.+2.*g/(g+1.)*(mach*mach-1.)*(1.+pi/p0));
            double old_r=density*((g+1.)*(old_p+pi)/((g-1.)*(p0+pi))+1.)/((g+1.)/(g-1.)+(old_p+pi)/(p0+pi));
            double old_u=velocity+c/g*(old_p/p0-1.)*p0/(p0+pi)/std::sqrt((g+1.)/(2.*g)*(old_p/p0-1.)*p0/(p0+pi)+1.);
            double compat[5];StiffenedShockConditions(mach,velocity,p0,density,g,pi,compat);
            near(compat[0],old_r,"legacy density");near(compat[1],old_u,"legacy velocity");near(compat[2],old_p,"legacy pressure");++legacy;
        }
        ++cases;
    }
    // The generalized form also admits zero/negative upstream pressure when
    // p+p_infinity remains positive, without dividing by p itself.
    for(double pressure:{0.,-1.}){
        double out[5];StiffenedShockConditions(2.,0.,pressure,1.,2.5,20.,.01,-4.,out);
        if(!std::isfinite(out[3]))throw std::runtime_error("shifted pressure branch");++cases;
    }
    std::cout<<"SHOCK_CHECKS {\"jump_cases\":"<<cases<<",\"legacy_cases\":"<<legacy<<"}\n";
    amrex::Finalize();
}
