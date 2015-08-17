/***********************************************************************
* Filename : rtmlib-modC.c
* Create : Lin 2015-07-13
* Description:This is the C code to implement the same functions in ./rtmlib-modF.f90 
* Modified   : 
* Revision of last commit: $Rev$ 
* Author of last commit  : $Author$ $date$
* licence :
$licence$
* **********************************************************************/

#include <stdio.h>
#include <omp.h>
#include "gpu.h"


void rtm_op3_c_(int nt, int nz, int nx, int zrec, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c);


void rtm_op3_gpu(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c);



void rtm_op3_c_(int nt, int nz, int nx, int zrec, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c){ //(nz,nx)
 
 
   float c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0;
    
    int it, iz, ix;
	
	//call gpu function for data passthrough test
	//rtm_op3_gpu(nt, nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	rtm_gpu_init(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
      
    //time step backward
    fprintf(stderr,"nt=%d, nz=%d, nx=%d\n",nt, nz, nx);
    for(it = nt-3; it>=0; it--){
        //
//      !Vx
//      Vx0(:,:,it) =  Vx0(:,:,it) + Vx0(:,:,(it+2))
	if((it == 1000)){
		fprintf(stderr, "\nit = %d\n", it);
		rtm_gpu_func(it, nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
      	}
	else{
	  for(ix=0;ix<nx;ix++){
            for(iz=0;iz<nz;iz++){
               Vx0[index3d(iz, ix, it)] = Vx0[index3d(iz, ix, it)] + Vx0[index3d(iz, ix, it+2)];
               Vz0[index3d(iz, ix, it)] = Vz0[index3d(iz, ix, it)] + Vz0[index3d(iz, ix, it+2)];
               sigmaxx0[index3d(iz, ix, it)] = sigmaxx0[index3d(iz, ix, it)] + sigmaxx0[index3d(iz, ix, it+2)];
               sigmazz0[index3d(iz, ix, it)] = sigmazz0[index3d(iz, ix, it)] + sigmazz0[index3d(iz, ix, it+2)];
               sigmaxz0[index3d(iz, ix, it)] = sigmaxz0[index3d(iz, ix, it)] + sigmaxz0[index3d(iz, ix, it+2)];
            }
        }
}
	
//        
////      D_c(nz, nx, Vx0(:,:,it),sigmaxx0(:,:,(it+1)), aux_m2m3_c(:,:))    
        for(ix=4; ix<nx-5; ix++)
            for(iz=4; iz<nz-5; iz++){
              Vx0[index3d(iz,ix+5,it)] = Vx0[index3d(iz,ix+5,it)] + aux_m2m3_c[index(iz,ix)]*c1*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix+4,it)] = Vx0[index3d(iz,ix+4,it)] + aux_m2m3_c[index(iz,ix)]*c2*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix+3,it)] = Vx0[index3d(iz,ix+3,it)] + aux_m2m3_c[index(iz,ix)]*c3*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix+2,it)] = Vx0[index3d(iz,ix+2,it)] + aux_m2m3_c[index(iz,ix)]*c4*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix+1,it)] = Vx0[index3d(iz,ix+1,it)] + aux_m2m3_c[index(iz,ix)]*c5*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix  ,it)] = Vx0[index3d(iz,ix  ,it)] - aux_m2m3_c[index(iz,ix)]*c5*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix-1,it)] = Vx0[index3d(iz,ix-1,it)] - aux_m2m3_c[index(iz,ix)]*c4*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix-2,it)] = Vx0[index3d(iz,ix-2,it)] - aux_m2m3_c[index(iz,ix)]*c3*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix-3,it)] = Vx0[index3d(iz,ix-3,it)] - aux_m2m3_c[index(iz,ix)]*c2*sigmaxx0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix-4,it)] = Vx0[index3d(iz,ix-4,it)] - aux_m2m3_c[index(iz,ix)]*c1*sigmaxx0[index3d(iz,ix,it+1)];
	}
////      F_c(nz, nx, Vx0(:,:,it),sigmazz0(:,:,(it+1)), aux_m2_c(:,:))
 
        for(ix=4; ix<nx-5 ; ix++)
            for(iz=4; iz<nz-5 ; iz++){
              Vx0[index3d(iz,ix+5,it) ]= Vx0[index3d(iz,ix+5,it) ]+ aux_m2_c[index(iz,ix)]*c1*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix+4,it)] = Vx0[index3d(iz,ix+4,it)] + aux_m2_c[index(iz,ix)]*c2*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix+3,it)] = Vx0[index3d(iz,ix+3,it)] + aux_m2_c[index(iz,ix)]*c3*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix+2,it)] = Vx0[index3d(iz,ix+2,it)] + aux_m2_c[index(iz,ix)]*c4*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix+1,it)] = Vx0[index3d(iz,ix+1,it)] + aux_m2_c[index(iz,ix)]*c5*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix  ,it)] = Vx0[index3d(iz,ix  ,it)] - aux_m2_c[index(iz,ix)]*c5*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix-1,it)] = Vx0[index3d(iz,ix-1,it)] - aux_m2_c[index(iz,ix)]*c4*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix-2,it)] = Vx0[index3d(iz,ix-2,it)] - aux_m2_c[index(iz,ix)]*c3*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix-3,it)] = Vx0[index3d(iz,ix-3,it)] - aux_m2_c[index(iz,ix)]*c2*sigmazz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz,ix-4,it)] = Vx0[index3d(iz,ix-4,it)] - aux_m2_c[index(iz,ix)]*c1*sigmazz0[index3d(iz,ix,it+1)];
         }


//      H_c(nz, nx, Vx0(:,:,it),sigmaxz0(:,:,(it+1)), aux_m3_c(:,:))
 
        for(ix=5; ix<nx-5 ; ix++)
            for(iz=5; iz<nz-5 ; iz++){
              Vx0[index3d(iz+4,ix,it)] = Vx0[index3d(iz+4,ix,it)] + aux_m3_c[index(iz,ix)]*c1*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz+3,ix,it)] = Vx0[index3d(iz+3,ix,it)] + aux_m3_c[index(iz,ix)]*c2*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz+2,ix,it)] = Vx0[index3d(iz+2,ix,it)] + aux_m3_c[index(iz,ix)]*c3*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz+1,ix,it)] = Vx0[index3d(iz+1,ix,it)] + aux_m3_c[index(iz,ix)]*c4*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz  ,ix,it)] = Vx0[index3d(iz  ,ix,it)] + aux_m3_c[index(iz,ix)]*c5*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz-1,ix,it)] = Vx0[index3d(iz-1,ix,it)] - aux_m3_c[index(iz,ix)]*c5*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz-2,ix,it)] = Vx0[index3d(iz-2,ix,it)] - aux_m3_c[index(iz,ix)]*c4*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz-3,ix,it)] = Vx0[index3d(iz-3,ix,it)] - aux_m3_c[index(iz,ix)]*c3*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz-4,ix,it)] = Vx0[index3d(iz-4,ix,it)] - aux_m3_c[index(iz,ix)]*c2*sigmaxz0[index3d(iz,ix,it+1)];
              Vx0[index3d(iz-5,ix,it)] = Vx0[index3d(iz-5,ix,it)] - aux_m3_c[index(iz,ix)]*c1*sigmaxz0[index3d(iz,ix,it+1)];
     
        }

////      E_c(nz, nx, Vz0(:,:,it),sigmaxx0(:,:,(it+1)), aux_m2_c(:,:))    
        for(ix=4; ix<nx-5 ; ix++)
            for(iz=4; iz<nz-5 ; iz++){
              Vz0[index3d(iz+5,ix,it)] = Vz0[index3d(iz+5,ix,it)] + aux_m2_c[index(iz,ix)]*c1*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz+4,ix,it)] = Vz0[index3d(iz+4,ix,it)] + aux_m2_c[index(iz,ix)]*c2*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz+3,ix,it)] = Vz0[index3d(iz+3,ix,it)] + aux_m2_c[index(iz,ix)]*c3*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz+2,ix,it)] = Vz0[index3d(iz+2,ix,it)] + aux_m2_c[index(iz,ix)]*c4*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz+1,ix,it)] = Vz0[index3d(iz+1,ix,it)] + aux_m2_c[index(iz,ix)]*c5*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz  ,ix,it)] = Vz0[index3d(iz  ,ix,it)] - aux_m2_c[index(iz,ix)]*c5*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz-1,ix,it)] = Vz0[index3d(iz-1,ix,it)] - aux_m2_c[index(iz,ix)]*c4*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz-2,ix,it)] = Vz0[index3d(iz-2,ix,it)] - aux_m2_c[index(iz,ix)]*c3*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz-3,ix,it)] = Vz0[index3d(iz-3,ix,it)] - aux_m2_c[index(iz,ix)]*c2*sigmaxx0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz-4,ix,it)] = Vz0[index3d(iz-4,ix,it)] - aux_m2_c[index(iz,ix)]*c1*sigmaxx0[index3d(iz,ix,it+1)];
         }
    
//
////      G_c(nz, nx, Vz0(:,:,it),sigmazz0(:,:,(it+1)), aux_m2m3_c(:,:))
//
       for(ix=4; ix<nx-5 ; ix++)
           for(iz=4; iz<nz-5 ; iz++){
             Vz0[index3d(iz+5,ix, it)] = Vz0[index3d(iz+5,ix,it)] + aux_m2m3_c[index(iz,ix)]*c1*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz+4,ix, it)] = Vz0[index3d(iz+4,ix,it)] + aux_m2m3_c[index(iz,ix)]*c2*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz+3,ix, it)] = Vz0[index3d(iz+3,ix,it)] + aux_m2m3_c[index(iz,ix)]*c3*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz+2,ix, it)] = Vz0[index3d(iz+2,ix,it)] + aux_m2m3_c[index(iz,ix)]*c4*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz+1,ix, it)] = Vz0[index3d(iz+1,ix,it)] + aux_m2m3_c[index(iz,ix)]*c5*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz  ,ix, it)] = Vz0[index3d(iz  ,ix,it)] - aux_m2m3_c[index(iz,ix)]*c5*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz-1,ix, it)] = Vz0[index3d(iz-1,ix,it)] - aux_m2m3_c[index(iz,ix)]*c4*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz-2,ix, it)] = Vz0[index3d(iz-2,ix,it)] - aux_m2m3_c[index(iz,ix)]*c3*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz-3,ix, it)] = Vz0[index3d(iz-3,ix,it)] - aux_m2m3_c[index(iz,ix)]*c2*sigmazz0[index3d(iz,ix,it+1)];
             Vz0[index3d(iz-4,ix, it)] = Vz0[index3d(iz-4,ix,it)] - aux_m2m3_c[index(iz,ix)]*c1*sigmazz0[index3d(iz,ix,it+1)];
        
           }
//
////      J_c(nz, nx, Vz0(:,:,it),sigmaxz0(:,:,(it+1)), aux_m3_c(:,:))
        for(ix=5; ix<nx-5 ; ix++)
            for(iz=5; iz<nz-5 ; iz++){
              Vz0[index3d(iz,ix+4,it)] = Vz0[index3d(iz,ix+4,it)] + aux_m3_c[index(iz,ix)]*c1*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix+3,it)] = Vz0[index3d(iz,ix+3,it)] + aux_m3_c[index(iz,ix)]*c2*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix+2,it)] = Vz0[index3d(iz,ix+2,it)] + aux_m3_c[index(iz,ix)]*c3*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix+1,it)] = Vz0[index3d(iz,ix+1,it)] + aux_m3_c[index(iz,ix)]*c4*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix  ,it)] = Vz0[index3d(iz,ix  ,it)] + aux_m3_c[index(iz,ix)]*c5*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix-1,it)] = Vz0[index3d(iz,ix-1,it)] - aux_m3_c[index(iz,ix)]*c5*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix-2,it)] = Vz0[index3d(iz,ix-2,it)] - aux_m3_c[index(iz,ix)]*c4*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix-3,it)] = Vz0[index3d(iz,ix-3,it)] - aux_m3_c[index(iz,ix)]*c3*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix-4,it)] = Vz0[index3d(iz,ix-4,it)] - aux_m3_c[index(iz,ix)]*c2*sigmaxz0[index3d(iz,ix,it+1)];
              Vz0[index3d(iz,ix-5,it)] = Vz0[index3d(iz,ix-5,it)] - aux_m3_c[index(iz,ix)]*c1*sigmaxz0[index3d(iz,ix,it+1)];
        
            }

////      Ax_c(nz, nx, sigmaxx0(:,:,it),Vx0(:,:,(it+1)), m1_x(:,:))
        for(ix=5; ix<nx-5 ; ix++)
            for(iz=4; iz<nz-5 ; iz++){
              sigmaxx0[index3d(iz,ix+4,it)] = sigmaxx0[index3d(iz,ix+4,it)] + m1_x[index(iz,ix)]*c1*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix+3,it)] = sigmaxx0[index3d(iz,ix+3,it)] + m1_x[index(iz,ix)]*c2*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix+2,it)] = sigmaxx0[index3d(iz,ix+2,it)] + m1_x[index(iz,ix)]*c3*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix+1,it)] = sigmaxx0[index3d(iz,ix+1,it)] + m1_x[index(iz,ix)]*c4*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix ,it)]  = sigmaxx0[index3d(iz,ix  ,it)] + m1_x[index(iz,ix)]*c5*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix-1,it)] = sigmaxx0[index3d(iz,ix-1,it)] - m1_x[index(iz,ix)]*c5*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix-2,it)] = sigmaxx0[index3d(iz,ix-2,it)] - m1_x[index(iz,ix)]*c4*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix-3,it)] = sigmaxx0[index3d(iz,ix-3,it)] - m1_x[index(iz,ix)]*c3*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix-4,it)] = sigmaxx0[index3d(iz,ix-4,it)] - m1_x[index(iz,ix)]*c2*Vx0[index3d(iz,ix,it+1)];
              sigmaxx0[index3d(iz,ix-5,it)] = sigmaxx0[index3d(iz,ix-5,it)] - m1_x[index(iz,ix)]*c1*Vx0[index3d(iz,ix,it+1)];
            }


//      Bz_c(nz, nx, sigmazz0(:,:,it),Vz0(:,:,(it+1)), m1_z)
        for(ix=4; ix<nx-5 ; ix++)
            for(iz=5; iz<nz-5 ; iz++){
              sigmazz0[index3d(iz+4,ix,it)] = sigmazz0[index3d(iz+4,ix,it)] + m1_z[index(iz,ix)]*c1*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz+3,ix,it)] = sigmazz0[index3d(iz+3,ix,it)] + m1_z[index(iz,ix)]*c2*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz+2,ix,it)] = sigmazz0[index3d(iz+2,ix,it)] + m1_z[index(iz,ix)]*c3*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz+1,ix,it)] = sigmazz0[index3d(iz+1,ix,it)] + m1_z[index(iz,ix)]*c4*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz  ,ix,it)] = sigmazz0[index3d(iz  ,ix,it)] + m1_z[index(iz,ix)]*c5*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz-1,ix,it)] = sigmazz0[index3d(iz-1,ix,it)] - m1_z[index(iz,ix)]*c5*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz-2,ix,it)] = sigmazz0[index3d(iz-2,ix,it)] - m1_z[index(iz,ix)]*c4*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz-3,ix,it)] = sigmazz0[index3d(iz-3,ix,it)] - m1_z[index(iz,ix)]*c3*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz-4,ix,it)] = sigmazz0[index3d(iz-4,ix,it)] - m1_z[index(iz,ix)]*c2*Vz0[index3d(iz,ix,it+1)];
              sigmazz0[index3d(iz-5,ix,it)] = sigmazz0[index3d(iz-5,ix,it)] - m1_z[index(iz,ix)]*c1*Vz0[index3d(iz,ix,it+1)];
          }
    

//      Bx_c(nz, nx, sigmaxz0(:,:,it),Vx0(:,:,(it+1)), m1_x(:,:))
        for(ix=5; ix<nx-5 ; ix++)
            for(iz=4; iz<nz-5 ; iz++){
              sigmaxz0[index3d(iz+5,ix,it)] = sigmaxz0[index3d(iz+5,ix,it)] + m1_x[index(iz,ix)]*c1*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz+4,ix,it)] = sigmaxz0[index3d(iz+4,ix,it)] + m1_x[index(iz,ix)]*c2*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz+3,ix,it)] = sigmaxz0[index3d(iz+3,ix,it)] + m1_x[index(iz,ix)]*c3*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz+2,ix,it)] = sigmaxz0[index3d(iz+2,ix,it)] + m1_x[index(iz,ix)]*c4*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz+1,ix,it)] = sigmaxz0[index3d(iz+1,ix,it)] + m1_x[index(iz,ix)]*c5*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz  ,ix,it)] = sigmaxz0[index3d(iz  ,ix,it)] - m1_x[index(iz,ix)]*c5*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz-1,ix,it)] = sigmaxz0[index3d(iz-1,ix,it)] - m1_x[index(iz,ix)]*c4*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz-2,ix,it)] = sigmaxz0[index3d(iz-2,ix,it)] - m1_x[index(iz,ix)]*c3*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz-3,ix,it)] = sigmaxz0[index3d(iz-3,ix,it)] - m1_x[index(iz,ix)]*c2*Vx0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz-4,ix,it)] = sigmaxz0[index3d(iz-4,ix,it)] - m1_x[index(iz,ix)]*c1*Vx0[index3d(iz,ix,it+1)];
           }
    

//      Az_c(nz, nx, sigmaxz0(:,:,it),Vz0(:,:,(it+1)), m1_z(:,:))
//
        for(ix=4; ix<nx-5 ; ix++)
            for(iz=5; iz<nz-5 ; iz++){
              sigmaxz0[index3d(iz,ix+5,it)] = sigmaxz0[index3d(iz,ix+5,it)] + m1_z[index(iz,ix)]*c1*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix+4,it)] = sigmaxz0[index3d(iz,ix+4,it)] + m1_z[index(iz,ix)]*c2*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix+3,it)] = sigmaxz0[index3d(iz,ix+3,it)] + m1_z[index(iz,ix)]*c3*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix+2,it)] = sigmaxz0[index3d(iz,ix+2,it)] + m1_z[index(iz,ix)]*c4*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix+1,it)] = sigmaxz0[index3d(iz,ix+1,it)] + m1_z[index(iz,ix)]*c5*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix  ,it)] = sigmaxz0[index3d(iz,ix  ,it)] - m1_z[index(iz,ix)]*c5*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix-1,it)] = sigmaxz0[index3d(iz,ix-1,it)] - m1_z[index(iz,ix)]*c4*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix-2,it)] = sigmaxz0[index3d(iz,ix-2,it)] - m1_z[index(iz,ix)]*c3*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix-3,it)] = sigmaxz0[index3d(iz,ix-3,it)] - m1_z[index(iz,ix)]*c2*Vz0[index3d(iz,ix,it+1)];
              sigmaxz0[index3d(iz,ix-4,it)] = sigmaxz0[index3d(iz,ix-4,it)] - m1_z[index(iz,ix)]*c1*Vz0[index3d(iz,ix,it+1)];
          }
//    
 	if((it%(nt/59)) == 0){
		fprintf(stderr,"#");
	}
//
	}
	fprintf(stderr,"]\n");
	//rtm_gpu_final(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	rtm_gpu_final(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);

}

void rtm_op3_gpu(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c){
	
	rtm_gpu_init(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
	//rtm_gpu_func(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
	rtm_gpu_final(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
}

