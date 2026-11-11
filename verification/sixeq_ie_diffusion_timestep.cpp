// SPDX-License-Identifier: GPL-3.0-or-later WITH GPL-3.0-linking-exception
// Copyright (C) 2026 The Regents of the University of Michigan
// This file is part of COMPAS, free software under the GNU General Public License,
// version 3 or later, with an additional permission under section 7. It comes with
// ABSOLUTELY NO WARRANTY. See COPYRIGHT and LICENSE for details.

// Independent periodic heat-equation updates on anisotropic 2D/3D grids.
#include <SixEqIE_DiffusionTimeStep.H>
#include <AMReX_ParmParse.H>
#include <vector>
#include <iostream>
#include <stdexcept>
int main(int argc,char**argv){
    amrex::Initialize(argc,argv);
    {
        std::string invalid;amrex::ParmParse().query("invalid",invalid);
        const amrex::Real width[3]={.2,.7,1.3};
        if(!invalid.empty()){
            SixEqIE_DiffusionTimeStep(invalid=="coefficient"?-1.:2.3,
                                     invalid=="safety"?1.01:.8,width,3);
            throw std::runtime_error("invalid diffusion input accepted");
        }
        for(int dimensions:{2,3}){
            const int nx=6,ny=8,nz=dimensions==3?4:1;
            const auto index=[=](int i,int j,int k){return ((k+nz)%nz*ny+(j+ny)%ny)*nx+(i+nx)%nx;};
            std::vector<double> initial(nx*ny*nz),state(initial.size()),next(initial.size());
            for(int k=0;k<nz;++k)for(int j=0;j<ny;++j)for(int i=0;i<nx;++i)
                initial[index(i,j,k)]=1.+((i+j+k)%2?.4:-.4);
            const auto update=[&](double dt){
                for(int k=0;k<nz;++k)for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){
                    const int p=index(i,j,k);double rhs=(state[index(i-1,j,k)]-2.*state[p]+state[index(i+1,j,k)])/(width[0]*width[0]);
                    rhs+=(state[index(i,j-1,k)]-2.*state[p]+state[index(i,j+1,k)])/(width[1]*width[1]);
                    if(dimensions==3)rhs+=(state[index(i,j,k-1)]-2.*state[p]+state[index(i,j,k+1)])/(width[2]*width[2]);
                    next[p]=state[p]+2.3*dt*rhs;
                }
                state.swap(next);
            };
            state=initial;
            const double safe=SixEqIE_DiffusionTimeStep(2.3,.8,width,dimensions);
            for(int step=0;step<64;++step){
                update(safe);
                for(auto T:state)if(!std::isfinite(T) || T<.6-1e-14 || T>1.4+1e-14)
                    throw std::runtime_error("diffusion bound violated the discrete maximum principle");
            }
            state=initial;update(1.01*SixEqIE_DiffusionTimeStep(2.3,1.,width,dimensions));
            bool unstable=false;for(auto T:state)unstable=unstable || T<.6 || T>1.4;
            if(!unstable)throw std::runtime_error("above-bound checkerboard did not expose instability");
        }
        std::cout<<"DIFFUSION_TIMESTEP {\"anisotropic_grids\":2,\"safe_steps_per_grid\":64,\"above_bound_rejections\":2}\n";
    }
    amrex::Finalize();
}
