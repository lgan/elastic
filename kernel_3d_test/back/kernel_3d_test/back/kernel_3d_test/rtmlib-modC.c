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
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include "gpu.h"

void func_print(int ny, int nx, int nz, float *y){
	int iy, iz, ix;
	for(iy = 0; iy < ny+10; iy++)
		for(ix = 0; ix < nx+10; ix++)
			for(iz = 0; iz < nz+10; iz++){
				fprintf(stderr, "%.6f\t", y[n3d_index_ex_ori(iz,ix, iy)]);
			}
	fprintf(stderr,"\n");
			
	}	


void func_check(int ny, int nx, int nz, float *a, float *b){
	int iy, iz, ix;
	int flag = 1;
	for(iy = 0; iy < ny; iy++)
		for(ix = 0; ix < nx; ix++)
			for(iz = 0; iz < nz; iz++){
				if(fabs(a[n3d_index_ex(iz,ix, iy)] - b[n3d_index_ex(iz,ix, iy)]) > 10e-14){
					fprintf(stderr, "%d,%d,%d, %.14f\n",iy,ix, iz,    a[n3d_index_ex(iz,ix, iy)] );
					fprintf(stderr, "%d,%d,%d, %.14f\n\n",iy, ix, iz, b[n3d_index_ex(iz,ix, iy)] );
					flag=0;}
			}
	if(flag) fprintf(stderr, "Check Okay\n");
			
}	



int main(){

	int it, iz, ix, iy;
	int nx = NX;
	int nz = NZ;
	int ny = NY;

	//time recording  1 for CPU, 2 for GPU
	float ctime, gtime, ctime_ps, gtime_ps, speedup; //cpu, gpu time; cpu, gpu time per step
	float gpu_kernel_time[3];	// GPU time, [0] for copy in, [1] for computation, [2] for copy out	
	struct timeval start1, end1;
	struct timeval start2, end2;

   	float c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0;


	fprintf(stderr, "***********************************************\n");
	fprintf(stderr, "Running 3D Elastic RTM Kernel for %d time steps\n", Steps_write_back);
	fprintf(stderr, "Data size: (Y = %d , X = %d , Z = %d)\n", ny, nx, nz);
	fprintf(stderr, "************************************************\n");

	//********Lin for debug********************************//
	float *debug = NULL;//(float*)malloc(sizeof(float)*(nx)*(nz)*(nt));
	
	//****LIN***** 
	//Temporary extended 3D arrays, later replaced by function inputs/outputs
	//Time step +1
	float * tmp = NULL;
	
	float * ex_Vx0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_Vz0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_Vy0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxx0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmazz0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmayy0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxy0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxz0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmayz0_in = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	

	//Time step +2
	float * ex_Vx0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_Vz0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_Vy0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxx0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmazz0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmayy0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxy0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxz0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmayz0_in1 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));


	//time step 0 and output
	float * ex_Vx0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_Vz0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_Vy0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxx0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmazz0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmayy0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxy0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmaxz0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * ex_sigmayz0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));

	
	float * gpu_ex_Vx0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * gpu_ex_Vz0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * gpu_ex_Vy0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * gpu_ex_sigmaxx0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * gpu_ex_sigmazz0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * gpu_ex_sigmayy0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * gpu_ex_sigmaxy0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * gpu_ex_sigmaxz0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float * gpu_ex_sigmayz0_out = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));

 
   
	//expaned arrays to store different Operators 
	float *ex_m2 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float *ex_m3 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float *ex_m2m3 = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float *ex_m1_x = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float *ex_m1_y = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));
	float *ex_m1_z = (float*)malloc(sizeof(float)*(ny+10)*(nx+10)*(nz+10));

	for(iy=0;iy<ny+10;iy++)
	    for(ix=0;ix<nx+10;ix++)
		for(iz=0;iz<nz+10;iz++){	
			ex_Vx0_in[n3d_index_ex_ori(iz,ix,iy)]=0;//((rand()%1000)/500. -1) ;
			ex_Vy0_in[n3d_index_ex_ori(iz,ix,iy)]=0;//((rand()%1000)/500. -1) ;
			ex_Vz0_in[n3d_index_ex_ori(iz,ix,iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxx0_in[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmazz0_in[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmayy0_in[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxy0_in[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxz0_in[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmayz0_in[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;

			ex_Vx0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_Vy0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_Vz0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxx0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmazz0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmayy0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxy0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxz0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmayz0_in1[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;

			ex_Vx0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_Vy0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_Vz0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxx0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmazz0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmayy0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxy0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmaxz0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_sigmayz0_out[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;


			gpu_ex_Vx0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_Vx0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			gpu_ex_Vy0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_Vy0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			gpu_ex_Vz0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_Vz0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			gpu_ex_sigmaxx0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_sigmaxx0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			gpu_ex_sigmayy0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_sigmayy0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			gpu_ex_sigmazz0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_sigmazz0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			gpu_ex_sigmaxy0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_sigmaxy0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			gpu_ex_sigmaxz0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_sigmaxz0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			gpu_ex_sigmayz0_out[n3d_index_ex_ori(iz,ix, iy)] =0;//ex_sigmayz0_out[n3d_index_ex_ori(iz,ix, iy)]  ;
			

			ex_m2[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_m3[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_m2m3[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_m1_x[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_m1_y[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;
			ex_m1_z[n3d_index_ex_ori(iz,ix, iy)]=0;//((rand()%1000)/500. -1) ;

	}
	

	for(iy=0;iy<ny;iy++)
	    for(ix=0;ix<nx;ix++)
		for(iz=0;iz<nz;iz++){	
			ex_Vx0_in[n3d_index_ex(iz,ix,iy)]=((rand()%1000)/500. -1) ;
			ex_Vy0_in[n3d_index_ex(iz,ix,iy)]=((rand()%1000)/500. -1) ;
			ex_Vz0_in[n3d_index_ex(iz,ix,iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxx0_in[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmazz0_in[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmayy0_in[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxy0_in[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxz0_in[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmayz0_in[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;

			ex_Vx0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_Vy0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_Vz0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxx0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmazz0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmayy0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxy0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxz0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmayz0_in1[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;

			ex_Vx0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_Vy0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_Vz0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxx0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmazz0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmayy0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxy0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmaxz0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_sigmayz0_out[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;


			gpu_ex_Vx0_out[n3d_index_ex(iz,ix, iy)] = ex_Vx0_out[n3d_index_ex(iz,ix, iy)]  ;
			gpu_ex_Vy0_out[n3d_index_ex(iz,ix, iy)] = ex_Vy0_out[n3d_index_ex(iz,ix, iy)]  ;
			gpu_ex_Vz0_out[n3d_index_ex(iz,ix, iy)] = ex_Vz0_out[n3d_index_ex(iz,ix, iy)]  ;
			gpu_ex_sigmaxx0_out[n3d_index_ex(iz,ix, iy)] = ex_sigmaxx0_out[n3d_index_ex(iz,ix, iy)]  ;
			gpu_ex_sigmayy0_out[n3d_index_ex(iz,ix, iy)] = ex_sigmayy0_out[n3d_index_ex(iz,ix, iy)]  ;
			gpu_ex_sigmazz0_out[n3d_index_ex(iz,ix, iy)] = ex_sigmazz0_out[n3d_index_ex(iz,ix, iy)]  ;
			gpu_ex_sigmaxy0_out[n3d_index_ex(iz,ix, iy)] = ex_sigmaxy0_out[n3d_index_ex(iz,ix, iy)]  ;
			gpu_ex_sigmaxz0_out[n3d_index_ex(iz,ix, iy)] = ex_sigmaxz0_out[n3d_index_ex(iz,ix, iy)]  ;
			gpu_ex_sigmayz0_out[n3d_index_ex(iz,ix, iy)] = ex_sigmayz0_out[n3d_index_ex(iz,ix, iy)]  ;
			

			ex_m2[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_m3[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_m2m3[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_m1_x[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_m1_y[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;
			ex_m1_z[n3d_index_ex(iz,ix, iy)]=((rand()%1000)/500. -1) ;

	}

fprintf(stderr, "Random Data Initallized ==========> OK\n");




//*********LIN(GPU Part)***************
	fprintf(stderr, "GPU Computation\n");


	rtm_gpu_init(ny,nz,nx);
	rtm_gpu_func(ny, nz, nx, 
		    ex_Vy0_in, ex_Vx0_in, ex_Vz0_in, ex_sigmayy0_in, ex_sigmaxx0_in, ex_sigmazz0_in, ex_sigmaxy0_in, ex_sigmaxz0_in, ex_sigmayz0_in,
		    ex_Vy0_in1, ex_Vx0_in1, ex_Vz0_in1, ex_sigmayy0_in1, ex_sigmaxx0_in1, ex_sigmazz0_in1, ex_sigmaxy0_in1, ex_sigmaxz0_in1, ex_sigmayz0_in1,
		    gpu_ex_Vy0_out, gpu_ex_Vx0_out, gpu_ex_Vz0_out, gpu_ex_sigmayy0_out, gpu_ex_sigmaxx0_out, gpu_ex_sigmazz0_out, gpu_ex_sigmaxy0_out, gpu_ex_sigmaxz0_out, gpu_ex_sigmayz0_out,
		    ex_m1_y, ex_m1_x, ex_m1_z, ex_m2, ex_m3, ex_m2m3,
	            debug, gpu_kernel_time);
	rtm_gpu_final();

	fprintf(stderr,"GPU Computing ==============> OK\n");

//*************************************

#ifdef DEBUG

gettimeofday(&start1, NULL);

//set GPU_start_step < 0 for CPU only computing
//
fprintf(stderr, "CPU Computation\n");


for(it = 0 ; it<Steps_write_back ; it++){

	
	gettimeofday(&start2, NULL);

  	for(iy=0;iy<ny;iy++)
            for(ix=0; ix<nx; ix++)
                for(iz=0; iz<nz; iz++){
  
         	ex_Vx0_out[n3d_index_ex(iz,ix  ,iy)] = ex_Vx0_out[n3d_index_ex(iz,ix  ,iy)]	+ ex_Vx0_in1[n3d_index_ex(iz, ix, iy)]	

									+ ex_m2m3[n3d_index_ex(iz,ix-5, iy)]*c1*ex_sigmaxx0_in[n3d_index_ex(iz,ix-5,iy)]							
							 		+ ex_m2m3[n3d_index_ex(iz,ix-4, iy)]*c2*ex_sigmaxx0_in[n3d_index_ex(iz,ix-4,iy)]		
									+ ex_m2m3[n3d_index_ex(iz,ix-3, iy)]*c3*ex_sigmaxx0_in[n3d_index_ex(iz,ix-3,iy)]	
									+ ex_m2m3[n3d_index_ex(iz,ix-2, iy)]*c4*ex_sigmaxx0_in[n3d_index_ex(iz,ix-2,iy)]	
									+ ex_m2m3[n3d_index_ex(iz,ix-1, iy)]*c5*ex_sigmaxx0_in[n3d_index_ex(iz,ix-1,iy)]	
									- ex_m2m3[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2m3[n3d_index_ex(iz,ix+1, iy)]*c4*ex_sigmaxx0_in[n3d_index_ex(iz,ix+1,iy)]	
									- ex_m2m3[n3d_index_ex(iz,ix+2, iy)]*c3*ex_sigmaxx0_in[n3d_index_ex(iz,ix+2,iy)]	
									- ex_m2m3[n3d_index_ex(iz,ix+3, iy)]*c2*ex_sigmaxx0_in[n3d_index_ex(iz,ix+3,iy)]	
									- ex_m2m3[n3d_index_ex(iz,ix+4, iy)]*c1*ex_sigmaxx0_in[n3d_index_ex(iz,ix+4,iy)]
	

									+ ex_m2[n3d_index_ex(iz,ix-5, iy)]*c1*ex_sigmayy0_in[n3d_index_ex(iz,ix-5,iy)]							
							 		+ ex_m2[n3d_index_ex(iz,ix-4, iy)]*c2*ex_sigmayy0_in[n3d_index_ex(iz,ix-4,iy)]		
									+ ex_m2[n3d_index_ex(iz,ix-3, iy)]*c3*ex_sigmayy0_in[n3d_index_ex(iz,ix-3,iy)]	
									+ ex_m2[n3d_index_ex(iz,ix-2, iy)]*c4*ex_sigmayy0_in[n3d_index_ex(iz,ix-2,iy)]	
									+ ex_m2[n3d_index_ex(iz,ix-1, iy)]*c5*ex_sigmayy0_in[n3d_index_ex(iz,ix-1,iy)]	
									- ex_m2[n3d_index_ex(iz,  ix, iy)]*c5*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2[n3d_index_ex(iz,ix+1, iy)]*c4*ex_sigmayy0_in[n3d_index_ex(iz,ix+1,iy)]	
									- ex_m2[n3d_index_ex(iz,ix+2, iy)]*c3*ex_sigmayy0_in[n3d_index_ex(iz,ix+2,iy)]	
									- ex_m2[n3d_index_ex(iz,ix+3, iy)]*c2*ex_sigmayy0_in[n3d_index_ex(iz,ix+3,iy)]	
									- ex_m2[n3d_index_ex(iz,ix+4, iy)]*c1*ex_sigmayy0_in[n3d_index_ex(iz,ix+4,iy)]	
	

									+ ex_m2[n3d_index_ex(iz,ix-5, iy)]*c1*ex_sigmazz0_in[n3d_index_ex(iz,ix-5,iy)]							
							 		+ ex_m2[n3d_index_ex(iz,ix-4, iy)]*c2*ex_sigmazz0_in[n3d_index_ex(iz,ix-4,iy)]		
									+ ex_m2[n3d_index_ex(iz,ix-3, iy)]*c3*ex_sigmazz0_in[n3d_index_ex(iz,ix-3,iy)]	
									+ ex_m2[n3d_index_ex(iz,ix-2, iy)]*c4*ex_sigmazz0_in[n3d_index_ex(iz,ix-2,iy)]	
									+ ex_m2[n3d_index_ex(iz,ix-1, iy)]*c5*ex_sigmazz0_in[n3d_index_ex(iz,ix-1,iy)]	
									- ex_m2[n3d_index_ex(iz,  ix, iy)]*c5*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2[n3d_index_ex(iz,ix+1, iy)]*c4*ex_sigmazz0_in[n3d_index_ex(iz,ix+1,iy)]	
									- ex_m2[n3d_index_ex(iz,ix+2, iy)]*c3*ex_sigmazz0_in[n3d_index_ex(iz,ix+2,iy)]	
									- ex_m2[n3d_index_ex(iz,ix+3, iy)]*c2*ex_sigmazz0_in[n3d_index_ex(iz,ix+3,iy)]	
									- ex_m2[n3d_index_ex(iz,ix+4, iy)]*c1*ex_sigmazz0_in[n3d_index_ex(iz,ix+4,iy)]	
	

									+ ex_m3[n3d_index_ex(iz,ix, iy-4)]*c1*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy-4)]		
									+ ex_m3[n3d_index_ex(iz,ix, iy-3)]*c2*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy-3)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy-2)]*c3*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy-2)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy-1)]*c4*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy-1)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+1)]*c5*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy+1)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+2)]*c4*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy+2)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+3)]*c3*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy+3)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+4)]*c2*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy+4)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+5)]*c1*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy+5)]							
	

									+ ex_m3[n3d_index_ex(iz-4,ix, iy)]*c1*ex_sigmaxz0_in[n3d_index_ex(iz-4,ix,iy)]		
									+ ex_m3[n3d_index_ex(iz-3,ix, iy)]*c2*ex_sigmaxz0_in[n3d_index_ex(iz-3,ix,iy)]	
									+ ex_m3[n3d_index_ex(iz-2,ix, iy)]*c3*ex_sigmaxz0_in[n3d_index_ex(iz-2,ix,iy)]	
									+ ex_m3[n3d_index_ex(iz-1,ix, iy)]*c4*ex_sigmaxz0_in[n3d_index_ex(iz-1,ix,iy)]	
									+ ex_m3[n3d_index_ex(iz,  ix, iy)]*c5*ex_sigmaxz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+1,ix, iy)]*c5*ex_sigmaxz0_in[n3d_index_ex(iz+1,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+2,ix, iy)]*c4*ex_sigmaxz0_in[n3d_index_ex(iz+2,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+3,ix, iy)]*c3*ex_sigmaxz0_in[n3d_index_ex(iz+3,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+4,ix, iy)]*c2*ex_sigmaxz0_in[n3d_index_ex(iz+4,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+5,ix, iy)]*c1*ex_sigmaxz0_in[n3d_index_ex(iz+5,ix,iy)]	;						
	

         	ex_Vy0_out[n3d_index_ex(iz,ix  ,iy)] = ex_Vy0_out[n3d_index_ex(iz,ix  ,iy)]	+ ex_Vy0_in1[n3d_index_ex(iz, ix, iy)]	

									+ ex_m2m3[n3d_index_ex(iz,ix, iy-5)]*c1*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy-5)]							
							 		+ ex_m2m3[n3d_index_ex(iz,ix, iy-4)]*c2*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy-4)]		
									+ ex_m2m3[n3d_index_ex(iz,ix, iy-3)]*c3*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy-3)]	
									+ ex_m2m3[n3d_index_ex(iz,ix, iy-2)]*c4*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy-2)]	
									+ ex_m2m3[n3d_index_ex(iz,ix, iy-1)]*c5*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy-1)]	
									- ex_m2m3[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2m3[n3d_index_ex(iz,ix, iy+1)]*c4*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy+1)]	
									- ex_m2m3[n3d_index_ex(iz,ix, iy+2)]*c3*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy+2)]	
									- ex_m2m3[n3d_index_ex(iz,ix, iy+3)]*c2*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy+3)]	
									- ex_m2m3[n3d_index_ex(iz,ix, iy+4)]*c1*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy+4)]
	

									+ ex_m2[n3d_index_ex(iz,ix, iy-5)]*c1*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy-5)]							
							 		+ ex_m2[n3d_index_ex(iz,ix, iy-4)]*c2*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy-4)]		
									+ ex_m2[n3d_index_ex(iz,ix, iy-3)]*c3*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy-3)]	
									+ ex_m2[n3d_index_ex(iz,ix, iy-2)]*c4*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy-2)]	
									+ ex_m2[n3d_index_ex(iz,ix, iy-1)]*c5*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy-1)]	
									- ex_m2[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2[n3d_index_ex(iz,ix, iy+1)]*c4*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy+1)]	
									- ex_m2[n3d_index_ex(iz,ix, iy+2)]*c3*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy+2)]	
									- ex_m2[n3d_index_ex(iz,ix, iy+3)]*c2*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy+3)]	
									- ex_m2[n3d_index_ex(iz,ix, iy+4)]*c1*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy+4)]	
	

									+ ex_m2[n3d_index_ex(iz,ix, iy-5)]*c1*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy-5)]							
							 		+ ex_m2[n3d_index_ex(iz,ix, iy-4)]*c2*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy-4)]		
									+ ex_m2[n3d_index_ex(iz,ix, iy-3)]*c3*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy-3)]	
									+ ex_m2[n3d_index_ex(iz,ix, iy-2)]*c4*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy-2)]	
									+ ex_m2[n3d_index_ex(iz,ix, iy-1)]*c5*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy-1)]	
									- ex_m2[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2[n3d_index_ex(iz,ix, iy+1)]*c4*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy+1)]	
									- ex_m2[n3d_index_ex(iz,ix, iy+2)]*c3*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy+2)]	
									- ex_m2[n3d_index_ex(iz,ix, iy+3)]*c2*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy+3)]	
									- ex_m2[n3d_index_ex(iz,ix, iy+4)]*c1*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy+4)]	
	

									+ ex_m3[n3d_index_ex(iz-4,ix, iy)]*c1*ex_sigmayz0_in[n3d_index_ex(iz-4,ix,iy)]		
									+ ex_m3[n3d_index_ex(iz-3,ix, iy)]*c2*ex_sigmayz0_in[n3d_index_ex(iz-3,ix,iy)]	
									+ ex_m3[n3d_index_ex(iz-2,ix, iy)]*c3*ex_sigmayz0_in[n3d_index_ex(iz-2,ix,iy)]	
									+ ex_m3[n3d_index_ex(iz-1,ix, iy)]*c4*ex_sigmayz0_in[n3d_index_ex(iz-1,ix,iy)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+1,ix, iy)]*c5*ex_sigmayz0_in[n3d_index_ex(iz+1,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+2,ix, iy)]*c4*ex_sigmayz0_in[n3d_index_ex(iz+2,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+3,ix, iy)]*c3*ex_sigmayz0_in[n3d_index_ex(iz+3,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+4,ix, iy)]*c2*ex_sigmayz0_in[n3d_index_ex(iz+4,ix,iy)]	
									- ex_m3[n3d_index_ex(iz+5,ix, iy)]*c1*ex_sigmayz0_in[n3d_index_ex(iz+5,ix,iy)]							
	

									+ ex_m3[n3d_index_ex(iz,ix-4, iy)]*c1*ex_sigmaxy0_in[n3d_index_ex(iz,ix-4,iy)]		
									+ ex_m3[n3d_index_ex(iz,ix-3, iy)]*c2*ex_sigmaxy0_in[n3d_index_ex(iz,ix-3,iy)]	
									+ ex_m3[n3d_index_ex(iz,ix-2, iy)]*c3*ex_sigmaxy0_in[n3d_index_ex(iz,ix-2,iy)]	
									+ ex_m3[n3d_index_ex(iz,ix-1, iy)]*c4*ex_sigmaxy0_in[n3d_index_ex(iz,ix-1,iy)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmaxy0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+1, iy)]*c5*ex_sigmaxy0_in[n3d_index_ex(iz,ix+1,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+2, iy)]*c4*ex_sigmaxy0_in[n3d_index_ex(iz,ix+2,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+3, iy)]*c3*ex_sigmaxy0_in[n3d_index_ex(iz,ix+3,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+4, iy)]*c2*ex_sigmaxy0_in[n3d_index_ex(iz,ix+4,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+5, iy)]*c1*ex_sigmaxy0_in[n3d_index_ex(iz,ix+5,iy)]	;						




         	ex_Vz0_out[n3d_index_ex(iz,ix  ,iy)] = ex_Vz0_out[n3d_index_ex(iz,ix  ,iy)]	+ ex_Vz0_in1[n3d_index_ex(iz, ix, iy)]	

									+ ex_m2m3[n3d_index_ex(iz-5,ix, iy)]*c1*ex_sigmazz0_in[n3d_index_ex(iz-5,ix,iy)]							
							 		+ ex_m2m3[n3d_index_ex(iz-4,ix, iy)]*c2*ex_sigmazz0_in[n3d_index_ex(iz-4,ix,iy)]		
									+ ex_m2m3[n3d_index_ex(iz-3,ix, iy)]*c3*ex_sigmazz0_in[n3d_index_ex(iz-3,ix,iy)]	
									+ ex_m2m3[n3d_index_ex(iz-2,ix, iy)]*c4*ex_sigmazz0_in[n3d_index_ex(iz-2,ix,iy)]	
									+ ex_m2m3[n3d_index_ex(iz-1,ix, iy)]*c5*ex_sigmazz0_in[n3d_index_ex(iz-1,ix,iy)]	
									- ex_m2m3[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmazz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2m3[n3d_index_ex(iz+1,ix, iy)]*c4*ex_sigmazz0_in[n3d_index_ex(iz+1,ix,iy)]	
									- ex_m2m3[n3d_index_ex(iz+2,ix, iy)]*c3*ex_sigmazz0_in[n3d_index_ex(iz+2,ix,iy)]	
									- ex_m2m3[n3d_index_ex(iz+3,ix, iy)]*c2*ex_sigmazz0_in[n3d_index_ex(iz+3,ix,iy)]	
									- ex_m2m3[n3d_index_ex(iz+4,ix, iy)]*c1*ex_sigmazz0_in[n3d_index_ex(iz+4,ix,iy)]
	

									+ ex_m2[n3d_index_ex(iz-5,ix, iy)]*c1*ex_sigmaxx0_in[n3d_index_ex(iz-5,ix,iy)]							
							 		+ ex_m2[n3d_index_ex(iz-4,ix, iy)]*c2*ex_sigmaxx0_in[n3d_index_ex(iz-4,ix,iy)]		
									+ ex_m2[n3d_index_ex(iz-3,ix, iy)]*c3*ex_sigmaxx0_in[n3d_index_ex(iz-3,ix,iy)]	
									+ ex_m2[n3d_index_ex(iz-2,ix, iy)]*c4*ex_sigmaxx0_in[n3d_index_ex(iz-2,ix,iy)]	
									+ ex_m2[n3d_index_ex(iz-1,ix, iy)]*c5*ex_sigmaxx0_in[n3d_index_ex(iz-1,ix,iy)]	
									- ex_m2[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmaxx0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2[n3d_index_ex(iz+1,ix, iy)]*c4*ex_sigmaxx0_in[n3d_index_ex(iz+1,ix,iy)]	
									- ex_m2[n3d_index_ex(iz+2,ix, iy)]*c3*ex_sigmaxx0_in[n3d_index_ex(iz+2,ix,iy)]	
									- ex_m2[n3d_index_ex(iz+3,ix, iy)]*c2*ex_sigmaxx0_in[n3d_index_ex(iz+3,ix,iy)]	
									- ex_m2[n3d_index_ex(iz+4,ix, iy)]*c1*ex_sigmaxx0_in[n3d_index_ex(iz+4,ix,iy)]
	

									+ ex_m2[n3d_index_ex(iz-5,ix, iy)]*c1*ex_sigmayy0_in[n3d_index_ex(iz-5,ix,iy)]							
							 		+ ex_m2[n3d_index_ex(iz-4,ix, iy)]*c2*ex_sigmayy0_in[n3d_index_ex(iz-4,ix,iy)]		
									+ ex_m2[n3d_index_ex(iz-3,ix, iy)]*c3*ex_sigmayy0_in[n3d_index_ex(iz-3,ix,iy)]	
									+ ex_m2[n3d_index_ex(iz-2,ix, iy)]*c4*ex_sigmayy0_in[n3d_index_ex(iz-2,ix,iy)]	
									+ ex_m2[n3d_index_ex(iz-1,ix, iy)]*c5*ex_sigmayy0_in[n3d_index_ex(iz-1,ix,iy)]	
									- ex_m2[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmayy0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m2[n3d_index_ex(iz+1,ix, iy)]*c4*ex_sigmayy0_in[n3d_index_ex(iz+1,ix,iy)]	
									- ex_m2[n3d_index_ex(iz+2,ix, iy)]*c3*ex_sigmayy0_in[n3d_index_ex(iz+2,ix,iy)]	
									- ex_m2[n3d_index_ex(iz+3,ix, iy)]*c2*ex_sigmayy0_in[n3d_index_ex(iz+3,ix,iy)]	
									- ex_m2[n3d_index_ex(iz+4,ix, iy)]*c1*ex_sigmayy0_in[n3d_index_ex(iz+4,ix,iy)]
	
									+ ex_m3[n3d_index_ex(iz,ix, iy-4)]*c1*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy-4)]		
									+ ex_m3[n3d_index_ex(iz,ix, iy-3)]*c2*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy-3)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy-2)]*c3*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy-2)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy-1)]*c4*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy-1)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+1)]*c5*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy+1)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+2)]*c4*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy+2)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+3)]*c3*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy+3)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+4)]*c2*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy+4)]	
									- ex_m3[n3d_index_ex(iz,ix, iy+5)]*c1*ex_sigmayz0_in[n3d_index_ex(iz,ix,iy+5)]							
	

									+ ex_m3[n3d_index_ex(iz,ix-4, iy)]*c1*ex_sigmaxz0_in[n3d_index_ex(iz,ix-4,iy)]		
									+ ex_m3[n3d_index_ex(iz,ix-3, iy)]*c2*ex_sigmaxz0_in[n3d_index_ex(iz,ix-3,iy)]	
									+ ex_m3[n3d_index_ex(iz,ix-2, iy)]*c3*ex_sigmaxz0_in[n3d_index_ex(iz,ix-2,iy)]	
									+ ex_m3[n3d_index_ex(iz,ix-1, iy)]*c4*ex_sigmaxz0_in[n3d_index_ex(iz,ix-1,iy)]	
									+ ex_m3[n3d_index_ex(iz,ix, iy)]  *c5*ex_sigmaxz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+1, iy)]*c5*ex_sigmaxz0_in[n3d_index_ex(iz,ix+1,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+2, iy)]*c4*ex_sigmaxz0_in[n3d_index_ex(iz,ix+2,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+3, iy)]*c3*ex_sigmaxz0_in[n3d_index_ex(iz,ix+3,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+4, iy)]*c2*ex_sigmaxz0_in[n3d_index_ex(iz,ix+4,iy)]	
									- ex_m3[n3d_index_ex(iz,ix+5, iy)]*c1*ex_sigmaxz0_in[n3d_index_ex(iz,ix+5,iy)]	;						


		

              ex_sigmaxx0_out[n3d_index_ex(iz,ix  ,iy)] = ex_sigmaxx0_out[n3d_index_ex(iz,ix  , iy)]	+ ex_sigmaxx0_in1[n3d_index_ex(iz,ix  , iy)] 
									+ ex_m1_x[n3d_index_ex(iz,ix-4, iy)]*c1*ex_Vx0_in[n3d_index_ex(iz,ix-4,iy)]		
									+ ex_m1_x[n3d_index_ex(iz,ix-3, iy)]*c2*ex_Vx0_in[n3d_index_ex(iz,ix-3,iy)]	
									+ ex_m1_x[n3d_index_ex(iz,ix-2, iy)]*c3*ex_Vx0_in[n3d_index_ex(iz,ix-2,iy)]	
									+ ex_m1_x[n3d_index_ex(iz,ix-1, iy)]*c4*ex_Vx0_in[n3d_index_ex(iz,ix-1,iy)]	
									+ ex_m1_x[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vx0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_x[n3d_index_ex(iz,ix+1, iy)]*c5*ex_Vx0_in[n3d_index_ex(iz,ix+1,iy)]	
									- ex_m1_x[n3d_index_ex(iz,ix+2, iy)]*c4*ex_Vx0_in[n3d_index_ex(iz,ix+2,iy)]	
									- ex_m1_x[n3d_index_ex(iz,ix+3, iy)]*c3*ex_Vx0_in[n3d_index_ex(iz,ix+3,iy)]	
									- ex_m1_x[n3d_index_ex(iz,ix+4, iy)]*c2*ex_Vx0_in[n3d_index_ex(iz,ix+4,iy)]	
									- ex_m1_x[n3d_index_ex(iz,ix+5, iy)]*c1*ex_Vx0_in[n3d_index_ex(iz,ix+5,iy)]	;						

	    
              ex_sigmayy0_out[n3d_index_ex(iz,ix  ,iy)] = ex_sigmayy0_out[n3d_index_ex(iz,ix  , iy)]	+ ex_sigmayy0_in1[n3d_index_ex(iz,ix  , iy)] 
									+ ex_m1_y[n3d_index_ex(iz,ix, iy-4)]*c1*ex_Vy0_in[n3d_index_ex(iz,ix,iy-4)]		
									+ ex_m1_y[n3d_index_ex(iz,ix, iy-3)]*c2*ex_Vy0_in[n3d_index_ex(iz,ix,iy-3)]	
									+ ex_m1_y[n3d_index_ex(iz,ix, iy-2)]*c3*ex_Vy0_in[n3d_index_ex(iz,ix,iy-2)]	
									+ ex_m1_y[n3d_index_ex(iz,ix, iy-1)]*c4*ex_Vy0_in[n3d_index_ex(iz,ix,iy-1)]	
									+ ex_m1_y[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vy0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_y[n3d_index_ex(iz,ix, iy+1)]*c5*ex_Vy0_in[n3d_index_ex(iz,ix,iy+1)]	
									- ex_m1_y[n3d_index_ex(iz,ix, iy+2)]*c4*ex_Vy0_in[n3d_index_ex(iz,ix,iy+2)]	
									- ex_m1_y[n3d_index_ex(iz,ix, iy+3)]*c3*ex_Vy0_in[n3d_index_ex(iz,ix,iy+3)]	
									- ex_m1_y[n3d_index_ex(iz,ix, iy+4)]*c2*ex_Vy0_in[n3d_index_ex(iz,ix,iy+4)]	
									- ex_m1_y[n3d_index_ex(iz,ix, iy+5)]*c1*ex_Vy0_in[n3d_index_ex(iz,ix,iy+5)]	;		


              ex_sigmazz0_out[n3d_index_ex(iz,ix  ,iy)] = ex_sigmazz0_out[n3d_index_ex(iz,ix  , iy)]	+ ex_sigmazz0_in1[n3d_index_ex(iz,ix  , iy)] 
									+ ex_m1_z[n3d_index_ex(iz-4,ix, iy)]*c1*ex_Vz0_in[n3d_index_ex(iz-4,ix,iy)]		
									+ ex_m1_z[n3d_index_ex(iz-3,ix, iy)]*c2*ex_Vz0_in[n3d_index_ex(iz-3,ix,iy)]	
									+ ex_m1_z[n3d_index_ex(iz-2,ix, iy)]*c3*ex_Vz0_in[n3d_index_ex(iz-2,ix,iy)]	
									+ ex_m1_z[n3d_index_ex(iz-1,ix, iy)]*c4*ex_Vz0_in[n3d_index_ex(iz-1,ix,iy)]	
									+ ex_m1_z[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_z[n3d_index_ex(iz+1,ix, iy)]*c5*ex_Vz0_in[n3d_index_ex(iz+1,ix,iy)]	
									- ex_m1_z[n3d_index_ex(iz+2,ix, iy)]*c4*ex_Vz0_in[n3d_index_ex(iz+2,ix,iy)]	
									- ex_m1_z[n3d_index_ex(iz+3,ix, iy)]*c3*ex_Vz0_in[n3d_index_ex(iz+3,ix,iy)]	
									- ex_m1_z[n3d_index_ex(iz+4,ix, iy)]*c2*ex_Vz0_in[n3d_index_ex(iz+4,ix,iy)]	
									- ex_m1_z[n3d_index_ex(iz+5,ix, iy)]*c1*ex_Vz0_in[n3d_index_ex(iz+5,ix,iy)]	;		 
	
	


              ex_sigmaxy0_out[n3d_index_ex(iz,ix  ,iy)] = ex_sigmaxy0_out[n3d_index_ex(iz,ix  , iy)]	+ ex_sigmaxy0_in1[n3d_index_ex(iz,ix  , iy)] 
									+ ex_m1_y[n3d_index_ex(iz,ix-4, iy)]*c1*ex_Vy0_in[n3d_index_ex(iz,ix-4,iy)]		
									+ ex_m1_y[n3d_index_ex(iz,ix-3, iy)]*c2*ex_Vy0_in[n3d_index_ex(iz,ix-3,iy)]	
									+ ex_m1_y[n3d_index_ex(iz,ix-2, iy)]*c3*ex_Vy0_in[n3d_index_ex(iz,ix-2,iy)]	
									+ ex_m1_y[n3d_index_ex(iz,ix-1, iy)]*c4*ex_Vy0_in[n3d_index_ex(iz,ix-1,iy)]	
									+ ex_m1_y[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vy0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_y[n3d_index_ex(iz,ix+1, iy)]*c5*ex_Vy0_in[n3d_index_ex(iz,ix+1,iy)]	
									- ex_m1_y[n3d_index_ex(iz,ix+2, iy)]*c4*ex_Vy0_in[n3d_index_ex(iz,ix+2,iy)]	
									- ex_m1_y[n3d_index_ex(iz,ix+3, iy)]*c3*ex_Vy0_in[n3d_index_ex(iz,ix+3,iy)]	
									- ex_m1_y[n3d_index_ex(iz,ix+4, iy)]*c2*ex_Vy0_in[n3d_index_ex(iz,ix+4,iy)]	
									- ex_m1_y[n3d_index_ex(iz,ix+5, iy)]*c1*ex_Vy0_in[n3d_index_ex(iz,ix+5,iy)]	

	    
									+ ex_m1_x[n3d_index_ex(iz,ix, iy-4)]*c1*ex_Vx0_in[n3d_index_ex(iz,ix,iy-4)]		
									+ ex_m1_x[n3d_index_ex(iz,ix, iy-3)]*c2*ex_Vx0_in[n3d_index_ex(iz,ix,iy-3)]	
									+ ex_m1_x[n3d_index_ex(iz,ix, iy-2)]*c3*ex_Vx0_in[n3d_index_ex(iz,ix,iy-2)]	
									+ ex_m1_x[n3d_index_ex(iz,ix, iy-1)]*c4*ex_Vx0_in[n3d_index_ex(iz,ix,iy-1)]	
									+ ex_m1_x[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vx0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_x[n3d_index_ex(iz,ix, iy+1)]*c5*ex_Vx0_in[n3d_index_ex(iz,ix,iy+1)]	
									- ex_m1_x[n3d_index_ex(iz,ix, iy+2)]*c4*ex_Vx0_in[n3d_index_ex(iz,ix,iy+2)]	
									- ex_m1_x[n3d_index_ex(iz,ix, iy+3)]*c3*ex_Vx0_in[n3d_index_ex(iz,ix,iy+3)]	
									- ex_m1_x[n3d_index_ex(iz,ix, iy+4)]*c2*ex_Vx0_in[n3d_index_ex(iz,ix,iy+4)]	
									- ex_m1_x[n3d_index_ex(iz,ix, iy+5)]*c1*ex_Vx0_in[n3d_index_ex(iz,ix,iy+5)]	;		


              ex_sigmaxz0_out[n3d_index_ex(iz,ix  ,iy)] = ex_sigmaxz0_out[n3d_index_ex(iz,ix  , iy)]	+ ex_sigmaxz0_in1[n3d_index_ex(iz,ix  , iy)] 
									+ ex_m1_x[n3d_index_ex(iz-4,ix, iy)]*c1*ex_Vx0_in[n3d_index_ex(iz-4,ix,iy)]		
									+ ex_m1_x[n3d_index_ex(iz-3,ix, iy)]*c2*ex_Vx0_in[n3d_index_ex(iz-3,ix,iy)]	
									+ ex_m1_x[n3d_index_ex(iz-2,ix, iy)]*c3*ex_Vx0_in[n3d_index_ex(iz-2,ix,iy)]	
									+ ex_m1_x[n3d_index_ex(iz-1,ix, iy)]*c4*ex_Vx0_in[n3d_index_ex(iz-1,ix,iy)]	
									+ ex_m1_x[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vx0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_x[n3d_index_ex(iz+1,ix, iy)]*c5*ex_Vx0_in[n3d_index_ex(iz+1,ix,iy)]	
									- ex_m1_x[n3d_index_ex(iz+2,ix, iy)]*c4*ex_Vx0_in[n3d_index_ex(iz+2,ix,iy)]	
									- ex_m1_x[n3d_index_ex(iz+3,ix, iy)]*c3*ex_Vx0_in[n3d_index_ex(iz+3,ix,iy)]	
									- ex_m1_x[n3d_index_ex(iz+4,ix, iy)]*c2*ex_Vx0_in[n3d_index_ex(iz+4,ix,iy)]	
									- ex_m1_x[n3d_index_ex(iz+5,ix, iy)]*c1*ex_Vx0_in[n3d_index_ex(iz+5,ix,iy)]	
							
									+ ex_m1_z[n3d_index_ex(iz,ix-4, iy)]*c1*ex_Vz0_in[n3d_index_ex(iz,ix-4,iy)]		
									+ ex_m1_z[n3d_index_ex(iz,ix-3, iy)]*c2*ex_Vz0_in[n3d_index_ex(iz,ix-3,iy)]	
									+ ex_m1_z[n3d_index_ex(iz,ix-2, iy)]*c3*ex_Vz0_in[n3d_index_ex(iz,ix-2,iy)]	
									+ ex_m1_z[n3d_index_ex(iz,ix-1, iy)]*c4*ex_Vz0_in[n3d_index_ex(iz,ix-1,iy)]	
									+ ex_m1_z[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_z[n3d_index_ex(iz,ix+1, iy)]*c5*ex_Vz0_in[n3d_index_ex(iz,ix+1,iy)]	
									- ex_m1_z[n3d_index_ex(iz,ix+2, iy)]*c4*ex_Vz0_in[n3d_index_ex(iz,ix+2,iy)]	
									- ex_m1_z[n3d_index_ex(iz,ix+3, iy)]*c3*ex_Vz0_in[n3d_index_ex(iz,ix+3,iy)]	
									- ex_m1_z[n3d_index_ex(iz,ix+4, iy)]*c2*ex_Vz0_in[n3d_index_ex(iz,ix+4,iy)]	
									- ex_m1_z[n3d_index_ex(iz,ix+5, iy)]*c1*ex_Vz0_in[n3d_index_ex(iz,ix+5,iy)]	;						

              ex_sigmayz0_out[n3d_index_ex(iz,ix  ,iy)] = ex_sigmayz0_out[n3d_index_ex(iz,ix  , iy)]	+ ex_sigmayz0_in1[n3d_index_ex(iz,ix  , iy)] 
									+ ex_m1_y[n3d_index_ex(iz-4,ix, iy)]*c1*ex_Vy0_in[n3d_index_ex(iz-4,ix,iy)]		
									+ ex_m1_y[n3d_index_ex(iz-3,ix, iy)]*c2*ex_Vy0_in[n3d_index_ex(iz-3,ix,iy)]	
									+ ex_m1_y[n3d_index_ex(iz-2,ix, iy)]*c3*ex_Vy0_in[n3d_index_ex(iz-2,ix,iy)]	
									+ ex_m1_y[n3d_index_ex(iz-1,ix, iy)]*c4*ex_Vy0_in[n3d_index_ex(iz-1,ix,iy)]	
									+ ex_m1_y[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vy0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_y[n3d_index_ex(iz+1,ix, iy)]*c5*ex_Vy0_in[n3d_index_ex(iz+1,ix,iy)]	
									- ex_m1_y[n3d_index_ex(iz+2,ix, iy)]*c4*ex_Vy0_in[n3d_index_ex(iz+2,ix,iy)]	
									- ex_m1_y[n3d_index_ex(iz+3,ix, iy)]*c3*ex_Vy0_in[n3d_index_ex(iz+3,ix,iy)]	
									- ex_m1_y[n3d_index_ex(iz+4,ix, iy)]*c2*ex_Vy0_in[n3d_index_ex(iz+4,ix,iy)]	
									- ex_m1_y[n3d_index_ex(iz+5,ix, iy)]*c1*ex_Vy0_in[n3d_index_ex(iz+5,ix,iy)]	
	
									+ ex_m1_z[n3d_index_ex(iz,ix, iy-4)]*c1*ex_Vz0_in[n3d_index_ex(iz,ix,iy-4)]		
									+ ex_m1_z[n3d_index_ex(iz,ix, iy-3)]*c2*ex_Vz0_in[n3d_index_ex(iz,ix,iy-3)]	
									+ ex_m1_z[n3d_index_ex(iz,ix, iy-2)]*c3*ex_Vz0_in[n3d_index_ex(iz,ix,iy-2)]	
									+ ex_m1_z[n3d_index_ex(iz,ix, iy-1)]*c4*ex_Vz0_in[n3d_index_ex(iz,ix,iy-1)]	
									+ ex_m1_z[n3d_index_ex(iz,ix, iy)]  *c5*ex_Vz0_in[n3d_index_ex(iz,ix,iy)]	
									- ex_m1_z[n3d_index_ex(iz,ix, iy+1)]*c5*ex_Vz0_in[n3d_index_ex(iz,ix,iy+1)]	
									- ex_m1_z[n3d_index_ex(iz,ix, iy+2)]*c4*ex_Vz0_in[n3d_index_ex(iz,ix,iy+2)]	
									- ex_m1_z[n3d_index_ex(iz,ix, iy+3)]*c3*ex_Vz0_in[n3d_index_ex(iz,ix,iy+3)]	
									- ex_m1_z[n3d_index_ex(iz,ix, iy+4)]*c2*ex_Vz0_in[n3d_index_ex(iz,ix,iy+4)]	
									- ex_m1_z[n3d_index_ex(iz,ix, iy+5)]*c1*ex_Vz0_in[n3d_index_ex(iz,ix,iy+5)]	;		

		}
		gettimeofday(&end2, NULL);

		ctime = (end2.tv_sec-start2.tv_sec)+(end2.tv_usec-start2.tv_usec)/1000000.0;
		fprintf(stderr, "CPU time at %d step:  %.8f\n",it+1, ctime);
	
		
		if(0 < Steps_write_back-1){
		fprintf(stderr, "CPU pointer chanced\n");
		tmp = ex_Vx0_out;
		ex_Vx0_out = ex_Vx0_in;
		ex_Vx0_in = tmp;

		tmp = ex_Vx0_out;
		ex_Vx0_out = ex_Vx0_in1;
		ex_Vx0_in1 = tmp; 

	

		tmp = ex_Vz0_out;
		ex_Vz0_out = ex_Vz0_in;
		ex_Vz0_in = tmp;

		tmp = ex_Vz0_out;
		ex_Vz0_out = ex_Vz0_in1;
		ex_Vz0_in1 = tmp; 


	
		tmp = ex_Vy0_out;
		ex_Vy0_out = ex_Vy0_in;
		ex_Vy0_in = tmp;

		tmp = ex_Vy0_out;
		ex_Vy0_out = ex_Vy0_in1;
		ex_Vy0_in1 = tmp; 

	

		tmp = ex_sigmaxx0_out;
		ex_sigmaxx0_out = ex_sigmaxx0_in;
		ex_sigmaxx0_in = tmp;

		tmp = ex_sigmaxx0_out;
		ex_sigmaxx0_out = ex_sigmaxx0_in1;
		ex_sigmaxx0_in1 = tmp; 



	
		tmp = ex_sigmazz0_out;
		ex_sigmazz0_out = ex_sigmazz0_in;
		ex_sigmazz0_in = tmp;

		tmp = ex_sigmazz0_out;
		ex_sigmazz0_out = ex_sigmazz0_in1;
		ex_sigmazz0_in1 = tmp; 



	
		tmp = ex_sigmayy0_out;
		ex_sigmayy0_out = ex_sigmayy0_in;
		ex_sigmayy0_in = tmp;

		tmp = ex_sigmayy0_out;
		ex_sigmayy0_out = ex_sigmayy0_in1;
		ex_sigmayy0_in1 = tmp; 


	
		tmp = ex_sigmaxy0_out;
		ex_sigmaxy0_out = ex_sigmaxy0_in;
		ex_sigmaxy0_in = tmp;

		tmp = ex_sigmaxy0_out;
		ex_sigmaxy0_out = ex_sigmaxy0_in1;
		ex_sigmaxy0_in1 = tmp; 



	
		tmp = ex_sigmaxz0_out;
		ex_sigmaxz0_out = ex_sigmaxz0_in;
		ex_sigmaxz0_in = tmp;

		tmp = ex_sigmaxz0_out;
		ex_sigmaxz0_out = ex_sigmaxz0_in1;
		ex_sigmaxz0_in1 = tmp; 


	
		tmp = ex_sigmayz0_out;
		ex_sigmayz0_out = ex_sigmayz0_in;
		ex_sigmayz0_in = tmp;
		
		tmp = ex_sigmayz0_out;
		ex_sigmayz0_out = ex_sigmayz0_in1;
		ex_sigmayz0_in1 = tmp; 


		}

	}
	fprintf(stderr,"CPU Computing ==============> OK\n");

	
	func_check(ny, nx, nz, ex_sigmaxx0_out, gpu_ex_sigmaxx0_out);
	func_check(ny, nx, nz, ex_sigmayy0_out, gpu_ex_sigmayy0_out);
	func_check(ny, nx, nz, ex_sigmazz0_out, gpu_ex_sigmazz0_out);
	func_check(ny, nx, nz, ex_Vx0_out, gpu_ex_Vx0_out);
	func_check(ny, nx, nz, ex_Vy0_out, gpu_ex_Vy0_out);
	func_check(ny, nx, nz, ex_Vz0_out, gpu_ex_Vz0_out);
	func_check(ny, nx, nz, ex_sigmaxy0_out, gpu_ex_sigmaxy0_out);
	func_check(ny, nx, nz, ex_sigmaxz0_out, gpu_ex_sigmaxz0_out);
	func_check(ny, nx, nz, ex_sigmayz0_out, gpu_ex_sigmayz0_out);


	gettimeofday(&end1, NULL);

	fprintf(stderr, "<<<<<<<<<<<<<<<<<PERFORMANCE PROFILING>>>>>>>>>>>>>>>>\n");

	ctime = (end1.tv_sec-start1.tv_sec)+(end1.tv_usec-start1.tv_usec)/1000000.0;
	fprintf(stderr, "CPU time for %d steps:  %.8f\n",Steps_write_back, ctime);
	
	fprintf(stderr, "Computing   Speedup: %.2f\n",(ctime/gpu_kernel_time[1])); 
	fprintf(stderr, "Application Speedup: %.2f\n",(ctime/(gpu_kernel_time[1]+gpu_kernel_time[0]+gpu_kernel_time[2]))); 
	
	fprintf(stderr, "<<<<<<<<<<<<<<<<<PERFORMANCE PROFILING>>>>>>>>>>>>>>>>\n");

#endif

	free(ex_Vx0_in);
	free(ex_Vz0_in);
	free(ex_Vy0_in);
	free(ex_sigmaxx0_in);
	free(ex_sigmazz0_in);
	free(ex_sigmayy0_in);
	free(ex_sigmaxy0_in);
	free(ex_sigmaxz0_in);
	free(ex_sigmayz0_in);

	free(ex_Vx0_in1);
	free(ex_Vz0_in1);
	free(ex_Vy0_in1);
	free(ex_sigmaxx0_in1);
	free(ex_sigmazz0_in1);
	free(ex_sigmayy0_in1);
	free(ex_sigmaxy0_in1);
	free(ex_sigmaxz0_in1);
	free(ex_sigmayz0_in1);


	free(ex_Vx0_out);
	free(ex_Vz0_out);
	free(ex_Vy0_out);
	free(ex_sigmaxx0_out);
	free(ex_sigmazz0_out);
	free(ex_sigmayy0_out);
	free(ex_sigmaxy0_out);
	free(ex_sigmaxz0_out);
	free(ex_sigmayz0_out);


	free(ex_m2);
	free(ex_m3);
	free(ex_m2m3);
	free(ex_m1_x);
	free(ex_m1_z);
	free(ex_m1_y);

	return 1;
}


