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
    
	int it, iz, ix;
	int iit, iiz, iix;
 
 
   	float c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0;
   
	//expaned arrays to store
	float *ex_aux_m2m3_c = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_aux_m2_c = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_aux_m3_c = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_sigmaxx0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_sigmazz0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_sigmaxz0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_Vx0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_Vz0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_m1_x = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_m1_z = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	
	//********Lin for debug********************************//
	float *debug = (float*)malloc(sizeof(float)*(nx)*(nz)*(nt));
	int fx, fz, ft;
	for(ft = 0;ft < nt; ft++)
	    for(fx=0; fx<nx; fx++)
	  	for(fz=0; fz<nz; fz++)
			debug[index3d(fz, fx, ft)] = 0;
	//********Lin for debug********************************//
		

	for(ix=0;ix<nx+10;ix++)
		for(iz=0;iz<nz+10;iz++){	
			ex_sigmaxx0[index_ex_ori(iz,ix)]=0.;
			ex_sigmazz0[index_ex_ori(iz,ix)]=0.;
			ex_sigmaxz0[index_ex_ori(iz,ix)]=0.;
			ex_aux_m2m3_c[index_ex_ori(iz,ix)]=0.;
			ex_aux_m2_c[index_ex_ori(iz,ix)]=0.;
			ex_aux_m3_c[index_ex_ori(iz,ix)]=0.;
			ex_Vx0[index_ex_ori(iz,ix)]=0.;
			ex_Vz0[index_ex_ori(iz,ix)]=0.;
			ex_m1_x[index_ex_ori(iz,ix)]=0.;
			ex_m1_z[index_ex_ori(iz,ix)]=0.;
	}

	for(ix=4;ix<nx-5;ix++)
		for(iz=4;iz<nz-5;iz++){	
			ex_aux_m2m3_c[index_ex(iz,ix)] = aux_m2m3_c[index(iz,ix)];
			ex_aux_m2_c[index_ex(iz,ix)] = aux_m2_c[index(iz,ix)];
			if(iz>4 && ix>4)
			ex_aux_m3_c[index_ex(iz,ix)] = aux_m3_c[index(iz,ix)];
			if(ix>4)
			ex_m1_x[index_ex(iz,ix)]=m1_x[index(iz,ix)];
			if(iz>4)
			ex_m1_z[index_ex(iz,ix)]=m1_z[index(iz,ix)];
	}
   

	//call gpu function for data passthrough test
	rtm_gpu_init(nt, nz, nx);
      
    //time step backward


       //
//      !Vx

for(it = nt-3; it>=0; it--){


if(it<2000 && it %500==0){

		fprintf(stderr, "\nit = %d\n", it);
		rtm_gpu_func(it, nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, ex_m1_x, ex_m1_z, ex_aux_m2_c, ex_aux_m3_c, ex_aux_m2m3_c, debug);
 
}else{
		//fprintf(stderr, "Using CPU Now\n");

 
      for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++){
              Vx0[index3d(iz,ix  ,it)] = Vx0[index3d(iz,ix  ,it)]	+ Vx0[index3d(iz, ix, it+2)]	
									+ ex_aux_m2m3_c[index_ex(iz,ix-5)]*c1*sigmaxx0[index3d(iz,ix-5,it+1)]							
							 		+ ex_aux_m2m3_c[index_ex(iz,ix-4)]*c2*sigmaxx0[index3d(iz,ix-4,it+1)]		
									+ ex_aux_m2m3_c[index_ex(iz,ix-3)]*c3*sigmaxx0[index3d(iz,ix-3,it+1)]	
									+ ex_aux_m2m3_c[index_ex(iz,ix-2)]*c4*sigmaxx0[index3d(iz,ix-2,it+1)]	
									+ ex_aux_m2m3_c[index_ex(iz,ix-1)]*c5*sigmaxx0[index3d(iz,ix-1,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix)]  *c5*sigmaxx0[index3d(iz,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix+1)]*c4*sigmaxx0[index3d(iz,ix+1,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix+2)]*c3*sigmaxx0[index3d(iz,ix+2,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix+3)]*c2*sigmaxx0[index3d(iz,ix+3,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix+4)]*c1*sigmaxx0[index3d(iz,ix+4,it+1)]
	

									+ ex_aux_m2_c[index_ex(iz,ix-5)]*c1*sigmazz0[index3d(iz,ix-5,it+1)]							
							 		+ ex_aux_m2_c[index_ex(iz,ix-4)]*c2*sigmazz0[index3d(iz,ix-4,it+1)]		
									+ ex_aux_m2_c[index_ex(iz,ix-3)]*c3*sigmazz0[index3d(iz,ix-3,it+1)]	
									+ ex_aux_m2_c[index_ex(iz,ix-2)]*c4*sigmazz0[index3d(iz,ix-2,it+1)]	
									+ ex_aux_m2_c[index_ex(iz,ix-1)]*c5*sigmazz0[index3d(iz,ix-1,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix)]  *c5*sigmazz0[index3d(iz,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix+1)]*c4*sigmazz0[index3d(iz,ix+1,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix+2)]*c3*sigmazz0[index3d(iz,ix+2,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix+3)]*c2*sigmazz0[index3d(iz,ix+3,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix+4)]*c1*sigmazz0[index3d(iz,ix+4,it+1)]	
	


									+ ex_aux_m3_c[index_ex(iz-4,ix)]*c1*sigmaxz0[index3d(iz-4,ix,it+1)]		
									+ ex_aux_m3_c[index_ex(iz-3,ix)]*c2*sigmaxz0[index3d(iz-3,ix,it+1)]	
									+ ex_aux_m3_c[index_ex(iz-2,ix)]*c3*sigmaxz0[index3d(iz-2,ix,it+1)]	
									+ ex_aux_m3_c[index_ex(iz-1,ix)]*c4*sigmaxz0[index3d(iz-1,ix,it+1)]	
									+ ex_aux_m3_c[index_ex(iz,ix)]  *c5*sigmaxz0[index3d(iz,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+1,ix)]*c5*sigmaxz0[index3d(iz+1,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+2,ix)]*c4*sigmaxz0[index3d(iz+2,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+3,ix)]*c3*sigmaxz0[index3d(iz+3,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+4,ix)]*c2*sigmaxz0[index3d(iz+4,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+5,ix)]*c1*sigmaxz0[index3d(iz+5,ix,it+1)]	;						
	
		//*********Lin for debug**************//
//		if(fabs(debug[index3d(iz,ix  ,it)] - Vz0[index3d(iz,ix  ,it)]) > 1e-7 ){
//			fprintf(stderr,"ix=%d, iz=%d, cpu = %.7f\n", ix,iz,Vz0[index3d(iz,ix  ,it)]);
//			fprintf(stderr,"ix=%d, iz=%d, gpu = %.7f\n\n", ix,iz,debug[index3d(iz,ix  ,it)]);}
//		//*********Lin for debug**************//
//	}

     	
             Vz0[index3d(iz,ix  ,it)] = Vz0[index3d(iz,ix  ,it)]  	+ Vz0[index3d(iz,ix  ,it+2)] 
									+ ex_aux_m2_c[index_ex(iz-5,ix)]*c1*sigmaxx0[index3d(iz-5,ix,it+1)]							
							 		+ ex_aux_m2_c[index_ex(iz-4,ix)]*c2*sigmaxx0[index3d(iz-4,ix,it+1)]		
									+ ex_aux_m2_c[index_ex(iz-3,ix)]*c3*sigmaxx0[index3d(iz-3,ix,it+1)]	
									+ ex_aux_m2_c[index_ex(iz-2,ix)]*c4*sigmaxx0[index3d(iz-2,ix,it+1)]	
									+ ex_aux_m2_c[index_ex(iz-1,ix)]*c5*sigmaxx0[index3d(iz-1,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix)]  *c5*sigmaxx0[index3d(iz,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz+1,ix)]*c4*sigmaxx0[index3d(iz+1,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz+2,ix)]*c3*sigmaxx0[index3d(iz+2,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz+3,ix)]*c2*sigmaxx0[index3d(iz+3,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz+4,ix)]*c1*sigmaxx0[index3d(iz+4,ix,it+1)]	//;
	

              //Vz0[index3d(iz,ix  ,it)] = Vz0[index3d(iz,ix  ,it)]	
              								+ ex_aux_m2m3_c[index_ex(iz-5,ix)]*c1*sigmazz0[index3d(iz-5,ix,it+1)]							
							 		+ ex_aux_m2m3_c[index_ex(iz-4,ix)]*c2*sigmazz0[index3d(iz-4,ix,it+1)]		
									+ ex_aux_m2m3_c[index_ex(iz-3,ix)]*c3*sigmazz0[index3d(iz-3,ix,it+1)]	
									+ ex_aux_m2m3_c[index_ex(iz-2,ix)]*c4*sigmazz0[index3d(iz-2,ix,it+1)]	
									+ ex_aux_m2m3_c[index_ex(iz-1,ix)]*c5*sigmazz0[index3d(iz-1,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix)]  *c5*sigmazz0[index3d(iz,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz+1,ix)]*c4*sigmazz0[index3d(iz+1,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz+2,ix)]*c3*sigmazz0[index3d(iz+2,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz+3,ix)]*c2*sigmazz0[index3d(iz+3,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz+4,ix)]*c1*sigmazz0[index3d(iz+4,ix,it+1)]	//;
	


//              Vz0[index3d(iz,ix  ,it)] = Vz0[index3d(iz,ix  ,it)]	
									+ ex_aux_m3_c[index_ex(iz,ix-4)]*c1*sigmaxz0[index3d(iz,ix-4,it+1)]		
									+ ex_aux_m3_c[index_ex(iz,ix-3)]*c2*sigmaxz0[index3d(iz,ix-3,it+1)]	
									+ ex_aux_m3_c[index_ex(iz,ix-2)]*c3*sigmaxz0[index3d(iz,ix-2,it+1)]	
									+ ex_aux_m3_c[index_ex(iz,ix-1)]*c4*sigmaxz0[index3d(iz,ix-1,it+1)]	
									+ ex_aux_m3_c[index_ex(iz,ix)]  *c5*sigmaxz0[index3d(iz,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+1)]*c5*sigmaxz0[index3d(iz,ix+1,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+2)]*c4*sigmaxz0[index3d(iz,ix+2,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+3)]*c3*sigmaxz0[index3d(iz,ix+3,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+4)]*c2*sigmaxz0[index3d(iz,ix+4,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+5)]*c1*sigmaxz0[index3d(iz,ix+5,it+1)]	;						

		

              sigmaxx0[index3d(iz,ix  ,it)] = sigmaxx0[index3d(iz,ix  ,it)]	+ sigmaxx0[index3d(iz,ix  ,it+2)] 
										+ ex_m1_x[index_ex(iz,ix-4)]*c1*Vx0[index3d(iz,ix-4,it+1)]		
										+ ex_m1_x[index_ex(iz,ix-3)]*c2*Vx0[index3d(iz,ix-3,it+1)]	
										+ ex_m1_x[index_ex(iz,ix-2)]*c3*Vx0[index3d(iz,ix-2,it+1)]	
										+ ex_m1_x[index_ex(iz,ix-1)]*c4*Vx0[index3d(iz,ix-1,it+1)]	
										+ ex_m1_x[index_ex(iz,ix)]  *c5*Vx0[index3d(iz,ix,it+1)]	
										- ex_m1_x[index_ex(iz,ix+1)]*c5*Vx0[index3d(iz,ix+1,it+1)]	
										- ex_m1_x[index_ex(iz,ix+2)]*c4*Vx0[index3d(iz,ix+2,it+1)]	
										- ex_m1_x[index_ex(iz,ix+3)]*c3*Vx0[index3d(iz,ix+3,it+1)]	
										- ex_m1_x[index_ex(iz,ix+4)]*c2*Vx0[index3d(iz,ix+4,it+1)]	
										- ex_m1_x[index_ex(iz,ix+5)]*c1*Vx0[index3d(iz,ix+5,it+1)]	;						
	

              sigmazz0[index3d(iz,ix  ,it)] = sigmazz0[index3d(iz,ix  ,it)]	+ sigmazz0[index3d(iz,ix  ,it+2)] 
										+ ex_m1_z[index_ex(iz-4,ix)]*c1*Vz0[index3d(iz-4,ix,it+1)]		
										+ ex_m1_z[index_ex(iz-3,ix)]*c2*Vz0[index3d(iz-3,ix,it+1)]	
										+ ex_m1_z[index_ex(iz-2,ix)]*c3*Vz0[index3d(iz-2,ix,it+1)]	
										+ ex_m1_z[index_ex(iz-1,ix)]*c4*Vz0[index3d(iz-1,ix,it+1)]	
										+ ex_m1_z[index_ex(iz,ix)]  *c5*Vz0[index3d(iz,ix,it+1)]	
										- ex_m1_z[index_ex(iz+1,ix)]*c5*Vz0[index3d(iz+1,ix,it+1)]	
										- ex_m1_z[index_ex(iz+2,ix)]*c4*Vz0[index3d(iz+2,ix,it+1)]	
										- ex_m1_z[index_ex(iz+3,ix)]*c3*Vz0[index3d(iz+3,ix,it+1)]	
										- ex_m1_z[index_ex(iz+4,ix)]*c2*Vz0[index3d(iz+4,ix,it+1)]	
										- ex_m1_z[index_ex(iz+5,ix)]*c1*Vz0[index3d(iz+5,ix,it+1)]	;						
	

             sigmaxz0[index3d(iz,ix  ,it)] = sigmaxz0[index3d(iz,ix  ,it)]	+ sigmaxz0[index3d(iz,ix  ,it+2)]	 
										+ ex_m1_x[index_ex(iz-5,ix)]*c1*Vx0[index3d(iz-5,ix,it+1)]							
							 			+ ex_m1_x[index_ex(iz-4,ix)]*c2*Vx0[index3d(iz-4,ix,it+1)]		
										+ ex_m1_x[index_ex(iz-3,ix)]*c3*Vx0[index3d(iz-3,ix,it+1)]	
										+ ex_m1_x[index_ex(iz-2,ix)]*c4*Vx0[index3d(iz-2,ix,it+1)]	
										+ ex_m1_x[index_ex(iz-1,ix)]*c5*Vx0[index3d(iz-1,ix,it+1)]	
										- ex_m1_x[index_ex(iz,ix)]  *c5*Vx0[index3d(iz,ix,it+1)]	
										- ex_m1_x[index_ex(iz+1,ix)]*c4*Vx0[index3d(iz+1,ix,it+1)]	
										- ex_m1_x[index_ex(iz+2,ix)]*c3*Vx0[index3d(iz+2,ix,it+1)]	
										- ex_m1_x[index_ex(iz+3,ix)]*c2*Vx0[index3d(iz+3,ix,it+1)]	
										- ex_m1_x[index_ex(iz+4,ix)]*c1*Vx0[index3d(iz+4,ix,it+1)]	//;
	
        
	      //sigmaxz0[index3d(iz,ix  ,it)] = sigmaxz0[index3d(iz,ix  ,it)]	
										+ ex_m1_z[index_ex(iz,ix-5)]*c1*Vz0[index3d(iz,ix-5,it+1)]							
							 			+ ex_m1_z[index_ex(iz,ix-4)]*c2*Vz0[index3d(iz,ix-4,it+1)]		
										+ ex_m1_z[index_ex(iz,ix-3)]*c3*Vz0[index3d(iz,ix-3,it+1)]	
										+ ex_m1_z[index_ex(iz,ix-2)]*c4*Vz0[index3d(iz,ix-2,it+1)]	
										+ ex_m1_z[index_ex(iz,ix-1)]*c5*Vz0[index3d(iz,ix-1,it+1)]	
										- ex_m1_z[index_ex(iz,ix)]  *c5*Vz0[index3d(iz,ix,it+1)]	
										- ex_m1_z[index_ex(iz,ix+1)]*c4*Vz0[index3d(iz,ix+1,it+1)]	
										- ex_m1_z[index_ex(iz,ix+2)]*c3*Vz0[index3d(iz,ix+2,it+1)]	
										- ex_m1_z[index_ex(iz,ix+3)]*c2*Vz0[index3d(iz,ix+3,it+1)]	
										- ex_m1_z[index_ex(iz,ix+4)]*c1*Vz0[index3d(iz,ix+4,it+1)]	;
	
	}
}
//    
 	if((it%(nt/59)) == 0){
		fprintf(stderr,"#");
	}
//
	}
	fprintf(stderr,"]\n");
	rtm_gpu_final();
}

void rtm_op3_gpu(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c){
	
	rtm_gpu_init(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
	//rtm_gpu_func(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
	rtm_gpu_final(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
}

