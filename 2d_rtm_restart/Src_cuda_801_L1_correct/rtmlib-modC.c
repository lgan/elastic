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
#include <sys/time.h>
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

	//time recording  1 for CPU, 2 for GPU
	double ctime, gtime, ctime_ps, gtime_ps, speedup; //cpu, gpu time; cpu, gpu time per step
	struct timeval start1, end1;
	struct timeval start2, end2;

	
 
 
   	float c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0;
   
	//expaned arrays to store
	float *ex_aux_m2m3_c = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_aux_m2_c = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_aux_m3_c = (float*)malloc(sizeof(float)*(nx+10)*(nz+10));
	float *ex_sigmaxx0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	float *ex_sigmazz0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	float *ex_sigmaxz0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	float *ex_Vx0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	float *ex_Vz0 = (float*)malloc(sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	float *ex_m1_x = (float*)malloc(sizeof(float)*(nx+10)*(nz+10)*(nt+10));
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
	
	for(it=0;it<nt;it++)
		for(ix=0;ix<nx;ix++)
			for(iz=0;iz<nz;iz++){
			
			ex_sigmaxx0[index3d_ex(iz,ix,it)]=sigmaxx0[index3d(iz,ix,it)];
			ex_sigmazz0[index3d_ex(iz,ix,it)]=sigmazz0[index3d(iz,ix,it)];
			ex_sigmaxz0[index3d_ex(iz,ix,it)]=sigmaxz0[index3d(iz,ix,it)];
			ex_Vx0[index3d_ex(iz,ix,it)]=Vx0[index3d(iz,ix,it)];
			ex_Vz0[index3d_ex(iz,ix,it)]=Vz0[index3d(iz,ix,it)];
			
	}

   

	//call gpu function for data passthrough test
	rtm_gpu_init(nt, nz, nx);
      
    //time step backward


       //
//      !Vx
gettimeofday(&start1, NULL);

//set GPU_start_step < 0 for CPU only computing
for(it = nt-3; (it>GPU_start_step && it>=0); it--){

  for(ix=0; ix<nx; ix++)
            for(iz=0; iz<nz; iz++){
  
         ex_Vx0[index3d_ex(iz,ix  ,it)] = ex_Vx0[index3d_ex(iz,ix  ,it)]	+ ex_Vx0[index3d_ex(iz, ix, it+2)]	
									+ ex_aux_m2m3_c[index_ex(iz,ix-5)]*c1*ex_sigmaxx0[index3d_ex(iz,ix-5,it+1)]							
							 		+ ex_aux_m2m3_c[index_ex(iz,ix-4)]*c2*ex_sigmaxx0[index3d_ex(iz,ix-4,it+1)]		
									+ ex_aux_m2m3_c[index_ex(iz,ix-3)]*c3*ex_sigmaxx0[index3d_ex(iz,ix-3,it+1)]	
									+ ex_aux_m2m3_c[index_ex(iz,ix-2)]*c4*ex_sigmaxx0[index3d_ex(iz,ix-2,it+1)]	
									+ ex_aux_m2m3_c[index_ex(iz,ix-1)]*c5*ex_sigmaxx0[index3d_ex(iz,ix-1,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix)]  *c5*ex_sigmaxx0[index3d_ex(iz,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix+1)]*c4*ex_sigmaxx0[index3d_ex(iz,ix+1,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix+2)]*c3*ex_sigmaxx0[index3d_ex(iz,ix+2,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix+3)]*c2*ex_sigmaxx0[index3d_ex(iz,ix+3,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix+4)]*c1*ex_sigmaxx0[index3d_ex(iz,ix+4,it+1)]
	

									+ ex_aux_m2_c[index_ex(iz,ix-5)]*c1*ex_sigmazz0[index3d_ex(iz,ix-5,it+1)]							
							 		+ ex_aux_m2_c[index_ex(iz,ix-4)]*c2*ex_sigmazz0[index3d_ex(iz,ix-4,it+1)]		
									+ ex_aux_m2_c[index_ex(iz,ix-3)]*c3*ex_sigmazz0[index3d_ex(iz,ix-3,it+1)]	
									+ ex_aux_m2_c[index_ex(iz,ix-2)]*c4*ex_sigmazz0[index3d_ex(iz,ix-2,it+1)]	
									+ ex_aux_m2_c[index_ex(iz,ix-1)]*c5*ex_sigmazz0[index3d_ex(iz,ix-1,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix)]  *c5*ex_sigmazz0[index3d_ex(iz,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix+1)]*c4*ex_sigmazz0[index3d_ex(iz,ix+1,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix+2)]*c3*ex_sigmazz0[index3d_ex(iz,ix+2,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix+3)]*c2*ex_sigmazz0[index3d_ex(iz,ix+3,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix+4)]*c1*ex_sigmazz0[index3d_ex(iz,ix+4,it+1)]	
	


									+ ex_aux_m3_c[index_ex(iz-4,ix)]*c1*ex_sigmaxz0[index3d_ex(iz-4,ix,it+1)]		
									+ ex_aux_m3_c[index_ex(iz-3,ix)]*c2*ex_sigmaxz0[index3d_ex(iz-3,ix,it+1)]	
									+ ex_aux_m3_c[index_ex(iz-2,ix)]*c3*ex_sigmaxz0[index3d_ex(iz-2,ix,it+1)]	
									+ ex_aux_m3_c[index_ex(iz-1,ix)]*c4*ex_sigmaxz0[index3d_ex(iz-1,ix,it+1)]	
									+ ex_aux_m3_c[index_ex(iz,ix)]  *c5*ex_sigmaxz0[index3d_ex(iz,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+1,ix)]*c5*ex_sigmaxz0[index3d_ex(iz+1,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+2,ix)]*c4*ex_sigmaxz0[index3d_ex(iz+2,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+3,ix)]*c3*ex_sigmaxz0[index3d_ex(iz+3,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+4,ix)]*c2*ex_sigmaxz0[index3d_ex(iz+4,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz+5,ix)]*c1*ex_sigmaxz0[index3d_ex(iz+5,ix,it+1)]	;						
	

             ex_Vz0[index3d_ex(iz,ix  ,it)] = ex_Vz0[index3d_ex(iz,ix  ,it)]  	+ ex_Vz0[index3d_ex(iz,ix  ,it+2)] 
									+ ex_aux_m2_c[index_ex(iz-5,ix)]*c1*ex_sigmaxx0[index3d_ex(iz-5,ix,it+1)]							
							 		+ ex_aux_m2_c[index_ex(iz-4,ix)]*c2*ex_sigmaxx0[index3d_ex(iz-4,ix,it+1)]		
									+ ex_aux_m2_c[index_ex(iz-3,ix)]*c3*ex_sigmaxx0[index3d_ex(iz-3,ix,it+1)]	
									+ ex_aux_m2_c[index_ex(iz-2,ix)]*c4*ex_sigmaxx0[index3d_ex(iz-2,ix,it+1)]	
									+ ex_aux_m2_c[index_ex(iz-1,ix)]*c5*ex_sigmaxx0[index3d_ex(iz-1,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz,ix)]  *c5*ex_sigmaxx0[index3d_ex(iz,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz+1,ix)]*c4*ex_sigmaxx0[index3d_ex(iz+1,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz+2,ix)]*c3*ex_sigmaxx0[index3d_ex(iz+2,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz+3,ix)]*c2*ex_sigmaxx0[index3d_ex(iz+3,ix,it+1)]	
									- ex_aux_m2_c[index_ex(iz+4,ix)]*c1*ex_sigmaxx0[index3d_ex(iz+4,ix,it+1)]	//;
	

              //Vz0[index3d_ex(iz,ix  ,it)] = Vz0[index3d_ex(iz,ix  ,it)]	
              								+ ex_aux_m2m3_c[index_ex(iz-5,ix)]*c1*ex_sigmazz0[index3d_ex(iz-5,ix,it+1)]							
							 		+ ex_aux_m2m3_c[index_ex(iz-4,ix)]*c2*ex_sigmazz0[index3d_ex(iz-4,ix,it+1)]		
									+ ex_aux_m2m3_c[index_ex(iz-3,ix)]*c3*ex_sigmazz0[index3d_ex(iz-3,ix,it+1)]	
									+ ex_aux_m2m3_c[index_ex(iz-2,ix)]*c4*ex_sigmazz0[index3d_ex(iz-2,ix,it+1)]	
									+ ex_aux_m2m3_c[index_ex(iz-1,ix)]*c5*ex_sigmazz0[index3d_ex(iz-1,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz,ix)]  *c5*ex_sigmazz0[index3d_ex(iz,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz+1,ix)]*c4*ex_sigmazz0[index3d_ex(iz+1,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz+2,ix)]*c3*ex_sigmazz0[index3d_ex(iz+2,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz+3,ix)]*c2*ex_sigmazz0[index3d_ex(iz+3,ix,it+1)]	
									- ex_aux_m2m3_c[index_ex(iz+4,ix)]*c1*ex_sigmazz0[index3d_ex(iz+4,ix,it+1)]	//;
	


//              Vz0[index3d_ex(iz,ix  ,it)] = Vz0[index3d_ex(iz,ix  ,it)]	
									+ ex_aux_m3_c[index_ex(iz,ix-4)]*c1*ex_sigmaxz0[index3d_ex(iz,ix-4,it+1)]		
									+ ex_aux_m3_c[index_ex(iz,ix-3)]*c2*ex_sigmaxz0[index3d_ex(iz,ix-3,it+1)]	
									+ ex_aux_m3_c[index_ex(iz,ix-2)]*c3*ex_sigmaxz0[index3d_ex(iz,ix-2,it+1)]	
									+ ex_aux_m3_c[index_ex(iz,ix-1)]*c4*ex_sigmaxz0[index3d_ex(iz,ix-1,it+1)]	
									+ ex_aux_m3_c[index_ex(iz,ix)]  *c5*ex_sigmaxz0[index3d_ex(iz,ix,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+1)]*c5*ex_sigmaxz0[index3d_ex(iz,ix+1,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+2)]*c4*ex_sigmaxz0[index3d_ex(iz,ix+2,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+3)]*c3*ex_sigmaxz0[index3d_ex(iz,ix+3,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+4)]*c2*ex_sigmaxz0[index3d_ex(iz,ix+4,it+1)]	
									- ex_aux_m3_c[index_ex(iz,ix+5)]*c1*ex_sigmaxz0[index3d_ex(iz,ix+5,it+1)]	;						

		

              ex_sigmaxx0[index3d_ex(iz,ix  ,it)] = ex_sigmaxx0[index3d_ex(iz,ix  ,it)]	+ ex_sigmaxx0[index3d_ex(iz,ix  ,it+2)] 
										+ ex_m1_x[index_ex(iz,ix-4)]*c1*ex_Vx0[index3d_ex(iz,ix-4,it+1)]		
										+ ex_m1_x[index_ex(iz,ix-3)]*c2*ex_Vx0[index3d_ex(iz,ix-3,it+1)]	
										+ ex_m1_x[index_ex(iz,ix-2)]*c3*ex_Vx0[index3d_ex(iz,ix-2,it+1)]	
										+ ex_m1_x[index_ex(iz,ix-1)]*c4*ex_Vx0[index3d_ex(iz,ix-1,it+1)]	
										+ ex_m1_x[index_ex(iz,ix)]  *c5*ex_Vx0[index3d_ex(iz,ix,it+1)]	
										- ex_m1_x[index_ex(iz,ix+1)]*c5*ex_Vx0[index3d_ex(iz,ix+1,it+1)]	
										- ex_m1_x[index_ex(iz,ix+2)]*c4*ex_Vx0[index3d_ex(iz,ix+2,it+1)]	
										- ex_m1_x[index_ex(iz,ix+3)]*c3*ex_Vx0[index3d_ex(iz,ix+3,it+1)]	
										- ex_m1_x[index_ex(iz,ix+4)]*c2*ex_Vx0[index3d_ex(iz,ix+4,it+1)]	
										- ex_m1_x[index_ex(iz,ix+5)]*c1*ex_Vx0[index3d_ex(iz,ix+5,it+1)]	;						
	

              ex_sigmazz0[index3d_ex(iz,ix  ,it)] = ex_sigmazz0[index3d_ex(iz,ix  ,it)]	+ ex_sigmazz0[index3d_ex(iz,ix  ,it+2)] 
										+ ex_m1_z[index_ex(iz-4,ix)]*c1*ex_Vz0[index3d_ex(iz-4,ix,it+1)]		
										+ ex_m1_z[index_ex(iz-3,ix)]*c2*ex_Vz0[index3d_ex(iz-3,ix,it+1)]	
										+ ex_m1_z[index_ex(iz-2,ix)]*c3*ex_Vz0[index3d_ex(iz-2,ix,it+1)]	
										+ ex_m1_z[index_ex(iz-1,ix)]*c4*ex_Vz0[index3d_ex(iz-1,ix,it+1)]	
										+ ex_m1_z[index_ex(iz,ix)]  *c5*ex_Vz0[index3d_ex(iz,ix,it+1)]	
										- ex_m1_z[index_ex(iz+1,ix)]*c5*ex_Vz0[index3d_ex(iz+1,ix,it+1)]	
										- ex_m1_z[index_ex(iz+2,ix)]*c4*ex_Vz0[index3d_ex(iz+2,ix,it+1)]	
										- ex_m1_z[index_ex(iz+3,ix)]*c3*ex_Vz0[index3d_ex(iz+3,ix,it+1)]	
										- ex_m1_z[index_ex(iz+4,ix)]*c2*ex_Vz0[index3d_ex(iz+4,ix,it+1)]	
										- ex_m1_z[index_ex(iz+5,ix)]*c1*ex_Vz0[index3d_ex(iz+5,ix,it+1)]	;						
	

             ex_sigmaxz0[index3d_ex(iz,ix  ,it)] = ex_sigmaxz0[index3d_ex(iz,ix  ,it)]	+ ex_sigmaxz0[index3d_ex(iz,ix  ,it+2)]	 
										+ ex_m1_x[index_ex(iz-5,ix)]*c1*ex_Vx0[index3d_ex(iz-5,ix,it+1)]							
							 			+ ex_m1_x[index_ex(iz-4,ix)]*c2*ex_Vx0[index3d_ex(iz-4,ix,it+1)]		
										+ ex_m1_x[index_ex(iz-3,ix)]*c3*ex_Vx0[index3d_ex(iz-3,ix,it+1)]	
										+ ex_m1_x[index_ex(iz-2,ix)]*c4*ex_Vx0[index3d_ex(iz-2,ix,it+1)]	
										+ ex_m1_x[index_ex(iz-1,ix)]*c5*ex_Vx0[index3d_ex(iz-1,ix,it+1)]	
										- ex_m1_x[index_ex(iz,ix)]  *c5*ex_Vx0[index3d_ex(iz,ix,it+1)]	
										- ex_m1_x[index_ex(iz+1,ix)]*c4*ex_Vx0[index3d_ex(iz+1,ix,it+1)]	
										- ex_m1_x[index_ex(iz+2,ix)]*c3*ex_Vx0[index3d_ex(iz+2,ix,it+1)]	
										- ex_m1_x[index_ex(iz+3,ix)]*c2*ex_Vx0[index3d_ex(iz+3,ix,it+1)]	
										- ex_m1_x[index_ex(iz+4,ix)]*c1*ex_Vx0[index3d_ex(iz+4,ix,it+1)]	//;
	
        
	      //sigmaxz0[index3d(iz,ix  ,it)] = sigmaxz0[index3d(iz,ix  ,it)]	
										+ ex_m1_z[index_ex(iz,ix-5)]*c1*ex_Vz0[index3d_ex(iz,ix-5,it+1)]							
							 			+ ex_m1_z[index_ex(iz,ix-4)]*c2*ex_Vz0[index3d_ex(iz,ix-4,it+1)]		
										+ ex_m1_z[index_ex(iz,ix-3)]*c3*ex_Vz0[index3d_ex(iz,ix-3,it+1)]	
										+ ex_m1_z[index_ex(iz,ix-2)]*c4*ex_Vz0[index3d_ex(iz,ix-2,it+1)]	
										+ ex_m1_z[index_ex(iz,ix-1)]*c5*ex_Vz0[index3d_ex(iz,ix-1,it+1)]	
										- ex_m1_z[index_ex(iz,ix)]  *c5*ex_Vz0[index3d_ex(iz,ix,it+1)]	
										- ex_m1_z[index_ex(iz,ix+1)]*c4*ex_Vz0[index3d_ex(iz,ix+1,it+1)]	
										- ex_m1_z[index_ex(iz,ix+2)]*c3*ex_Vz0[index3d_ex(iz,ix+2,it+1)]	
										- ex_m1_z[index_ex(iz,ix+3)]*c2*ex_Vz0[index3d_ex(iz,ix+3,it+1)]	
										- ex_m1_z[index_ex(iz,ix+4)]*c1*ex_Vz0[index3d_ex(iz,ix+4,it+1)]	;
	
	}

	if((it%(nt/59)) == 0){
		fprintf(stderr,"#");}

}

	
	fprintf(stderr,"]\n");

		gettimeofday(&end1, NULL);

if(GPU_start_step>=0){
		fprintf(stderr, "\nit = %d\n", it);
		gettimeofday(&start2, NULL);
		rtm_gpu_func(GPU_start_step, nt, nz, nx, ex_Vx0, ex_Vz0, ex_sigmaxx0, ex_sigmazz0, ex_sigmaxz0, ex_m1_x, ex_m1_z, ex_aux_m2_c, ex_aux_m3_c, ex_aux_m2m3_c, debug);
		gettimeofday(&end2, NULL);
		
}

	rtm_gpu_final();
	fprintf(stderr, "<<<<<<<<<<<<<<<<<PERFORMANCE PROFILING>>>>>>>>>>>>>>>>\n");
if((GPU_start_step < nt-3) && (GPU_start_step >=0)){
	//time computing
	gtime = (end2.tv_sec-start2.tv_sec)+(end2.tv_usec-start2.tv_usec)/1000000.0;
	ctime = (end1.tv_sec-start1.tv_sec)+(end1.tv_usec-start1.tv_usec)/1000000.0;
	//ctime = ctime - gtime;
	gtime_ps = gtime/(GPU_start_step+1);
	ctime_ps = ctime/(nt-3-GPU_start_step-1);
	speedup = ctime_ps/gtime_ps;
	fprintf(stderr, "CPU time per step %.8f\n",ctime_ps);
	fprintf(stderr, "GPU time per step %.8f\n",gtime_ps);
	fprintf(stderr, "SPEEDUP %.2f\n", speedup);
}else fprintf(stderr,"Set GPU_start_step to (0~nt-3) for speedup\n");
	fprintf(stderr, "<<<<<<<<<<<<<<<<<PERFORMANCE PROFILING>>>>>>>>>>>>>>>>\n");

	for(it=0;it<nt;it++)
		for(ix=0;ix<nx;ix++)
			for(iz=0;iz<nz;iz++){
			
			sigmaxx0[index3d(iz,ix,it)]=ex_sigmaxx0[index3d_ex(iz,ix,it)];
			sigmazz0[index3d(iz,ix,it)]=ex_sigmazz0[index3d_ex(iz,ix,it)];
			sigmaxz0[index3d(iz,ix,it)]=ex_sigmaxz0[index3d_ex(iz,ix,it)];
			Vx0[index3d(iz,ix,it)]=ex_Vx0[index3d_ex(iz,ix,it)];
			Vz0[index3d(iz,ix,it)]=ex_Vz0[index3d_ex(iz,ix,it)];
			
	}


}

void rtm_op3_gpu(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c){
	
	rtm_gpu_init(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
	//rtm_gpu_func(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
	rtm_gpu_final(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	
}

