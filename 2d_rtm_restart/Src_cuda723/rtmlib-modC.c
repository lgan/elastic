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

//void Ax_c_(int nx, int nz, float * model, float * data, float * m1_x);
//void Bx_c_(int nx, int nz, float * model, float * data, float * m1_x);
//void Az_c_(int nx, int nz, float * model, float * data, float * m1_z);
//void Bz_c_(int nx, int nz, float * model, float * data, float * m1_z);
//void D_c_(int nx, int nz, float * model, float * data, float * aux_m2m3_c);
//void E_c_(int nx, int nz, float * model, float * data, float * aux_m2_c);
//void F_c_(int nx, int nz, float * model, float * data, float * aux_m2_c);
//void G_c_(int nx, int nz, float * model, float * data, float * aux_m2m3_c);
//void H_c_(int nx, int nz, float * model, float * data, float * aux_m3_c);
//void J_c_(int nx, int nz, float * model, float * data, float * aux_m3_c);
#include <stdio.h>
#include <omp.h>
#include "gpu.h"
#define index(a,b) ((a)+(b)*nz)     //an index transpose happens here
#define indextx(a,b) ((a)+(b)*nt)     //an index transpose happens here
#define index3d(a,b,c) ( (c)*nx*nz + (b)*nz + (a))     //an index transpose happens here

void rtm_op3_c_(int nt, int nz, int nx, int zrec, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
//        float * Vx,  float * Vz,  float * sigmaxx,  float * sigmazz,  float * sigmaxz, //(nt, nx)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c);


void rtm_op3_gpu();



void rtm_op3_c_(int nt, int nz, int nx, int zrec, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
//        float * Vx,  float * Vz,  float * sigmaxx,  float * sigmazz,  float * sigmaxz, //(nt, nx)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c){ //(nz,nx)
 
 
   float c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0;
    
    int it, iz, ix;
   
    //time step backward
    fprintf(stderr,"nt=%d, nz=%d, nx=%d\n",nt, nz, nx);
    for(it = nt-3; it>=0; it--){
        //
//      !Vx
//      Vx0(:,:,it) =  Vx0(:,:,it) + Vx0(:,:,(it+2))
        for(ix=0;ix<nx;ix++){
            for(iz=0;iz<nz;iz++){
               Vx0[index3d(iz, ix, it)] = Vx0[index3d(iz, ix, it)] + Vx0[index3d(iz, ix, it+2)];
               Vz0[index3d(iz, ix, it)] = Vz0[index3d(iz, ix, it)] + Vz0[index3d(iz, ix, it+2)];
               sigmaxx0[index3d(iz, ix, it)] = sigmaxx0[index3d(iz, ix, it)] + sigmaxx0[index3d(iz, ix, it+2)];
               sigmazz0[index3d(iz, ix, it)] = sigmazz0[index3d(iz, ix, it)] + sigmazz0[index3d(iz, ix, it+2)];
               sigmaxz0[index3d(iz, ix, it)] = sigmaxz0[index3d(iz, ix, it)] + sigmaxz0[index3d(iz, ix, it+2)];
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
//              
// 	 Vx0[index3d(iz,ix  ,it)] = Vx0[index3d(iz,ix  ,it)] + aux_m2m3_c[index(iz,ix-5)]*c1*sigmaxx0[index3d(iz,ix-5,it+1)]
//							 	  + aux_m2m3_c[index(iz,ix-4)]*c2*sigmaxx0[index3d(iz,ix-4,it+1)];
//							 	  + aux_m2m3_c[index(iz,ix-3)]*c3*sigmaxx0[index3d(iz,ix-3,it+1)];
//							 	  + aux_m2m3_c[index(iz,ix-2)]*c4*sigmaxx0[index3d(iz,ix-2,it+1)];
//							 	  + aux_m2m3_c[index(iz,ix-1)]*c5*sigmaxx0[index3d(iz,ix-1,it+1)];
//							 	  - aux_m2m3_c[index(iz,ix)]*c5*sigmaxx0[index3d(iz,ix,it+1)];
//							 	  - aux_m2m3_c[index(iz,ix+1)]*c4*sigmaxx0[index3d(iz,ix+1,it+1)];
//							 	  - aux_m2m3_c[index(iz,ix+2)]*c3*sigmaxx0[index3d(iz,ix+2,it+1)];
//							 	  - aux_m2m3_c[index(iz,ix+3)]*c2*sigmaxx0[index3d(iz,ix+3,it+1)];
//							 	  - aux_m2m3_c[index(iz,ix+4)]*c1*sigmaxx0[index3d(iz,ix+4,it+1)];
//					         	
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
	rtm_op3_gpu();
//
}

void rtm_op3_gpu(){
	rtm_gpu_init();
}

