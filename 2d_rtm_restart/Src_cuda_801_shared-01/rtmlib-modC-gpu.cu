#include <stdio.h>
#include "gpu.h"


//propagation data
 float *g_Vx0;
 float *g_Vz0;
 float *g_sigmaxx0; 
 float *g_sigmazz0;
 float *g_sigmaxz0;


//extended propagation data residence in GPU device
__device__ float *g_ex_sigmaxx0; 
__device__ float *g_ex_sigmazz0; 
__device__ float *g_ex_sigmaxz0; 
__device__ float *g_ex_Vx0;
__device__ float *g_ex_Vz0;



//constant data, extended, with 10 more layers over the CPU version 
 float *g_ex_m1_x;
 float *g_ex_m1_z;
 float *g_ex_aux_m2_c; 
 float *g_ex_aux_m3_c; 
 float *g_ex_aux_m2m3_c; 

__global__ void rtm_gpu_kernel(int it,int nt, int nz, int nx,
        float * g_ex_Vx0, float * g_ex_Vz0, float * g_ex_sigmaxx0, float * g_ex_sigmazz0, float * g_ex_sigmaxz0, //(nz, nx, nt)
        float * g_ex_m1_x,float * g_ex_m1_z,float * g_ex_aux_m2_c, float * g_ex_aux_m3_c, float * g_ex_aux_m2m3_c);//(nz+10,	nx+10)


extern "C" void rtm_gpu_init(int nt, int nz, int nx) 
{
	//set cuda devices and put all data onto gpu memory
	
	cudaError_t cuda_ret;
     	cudaError_t err;

	//Set Device 
    	cuda_ret = cudaSetDevice(0);
	if(cuda_ret != cudaSuccess){
		fprintf(stderr, "Failed to Set The cuda Device !\n");
	}
	else{
		fprintf(stderr, "GPU Device Set ====> OK\n");
	}

	// data init
	cudaMalloc(&g_ex_Vx0,sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	cudaMalloc(&g_ex_Vz0,sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	cudaMalloc(&g_ex_sigmaxx0,sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	cudaMalloc(&g_ex_sigmaxx0,sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	cudaMalloc(&g_ex_sigmazz0,sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	cudaMalloc(&g_ex_sigmazz0,sizeof(float)*(nx+10)*(nz+10)*(nt+10));
	cudaMalloc(&g_ex_sigmaxz0,sizeof(float)*(nx+10)*(nz+10)*(nt+10));

//	cudaMalloc(&g_Vx0,sizeof(float)*nx*nz*nt);
//	cudaMalloc(&g_Vz0,sizeof(float)*nx*nz*nt);
//	cudaMalloc(&g_sigmaxx0,sizeof(float)*nx*nz*nt);
//	cudaMalloc(&g_sigmazz0,sizeof(float)*nx*nz*nt);
//	cudaMalloc(&g_sigmaxz0,sizeof(float)*nx*nz*nt);

	cudaMalloc(&g_ex_m1_x,sizeof(float)*(nx+10)*(nz+10));
	cudaMalloc(&g_ex_m1_z,sizeof(float)*(nx+10)*(nz+10));
	cudaMalloc(&g_ex_aux_m2_c,sizeof(float)*(nx+10)*(nz+10));
	cudaMalloc(&g_ex_aux_m3_c,sizeof(float)*(nx+10)*(nz+10));
	cudaMalloc(&g_ex_aux_m2m3_c,sizeof(float)*(nx+10)*(nz+10));


	cudaFuncSetCacheConfig(rtm_gpu_kernel,cudaFuncCachePreferShared);

	fprintf(stderr,"GPU Data Init ====> OK\n");

	// data copy

//	cudaMemcpy(g_Vx0, Vx0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_Vz0, Vz0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_sigmaxx0, sigmaxx0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_sigmaxz0, sigmaxz0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_sigmazz0, sigmazz0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_m1_x, m1_x, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_m1_z, m1_z, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_aux_m2_c, aux_m2_c, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_aux_m3_c, aux_m3_c, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
//	cudaMemcpy(g_aux_m2m3_c, aux_m2m3_c, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
//	fprintf(stderr,"Data Copy To GPU OK\n");
}



extern "C" void rtm_gpu_copy_in(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * ex_m1_x,float * ex_m1_z,float * ex_aux_m2_c, float * ex_aux_m3_c, float * ex_aux_m2m3_c)//(nz,	nx)
{
	// data copy

	cudaMemcpy(g_ex_Vx0, Vx0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_Vz0, Vz0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_sigmaxx0, sigmaxx0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_sigmaxz0, sigmaxz0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_sigmazz0, sigmazz0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_m1_x, ex_m1_x, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_m1_z, ex_m1_z, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_aux_m2_c, ex_aux_m2_c, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_aux_m3_c, ex_aux_m3_c, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyHostToDevice);
	cudaMemcpy(g_ex_aux_m2m3_c, ex_aux_m2m3_c, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyHostToDevice);
	
//	fprintf(stderr,"Copy out for debuging \n");
//	cudaMemcpy(ex_m1_x, g_ex_m1_x, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyDeviceToHost);
//	cudaMemcpy(ex_m1_z, g_ex_m1_z, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyDeviceToHost);
//	cudaMemcpy(ex_aux_m2_c, g_ex_aux_m2_c, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyDeviceToHost);
//	cudaMemcpy(ex_aux_m3_c, g_ex_aux_m3_c, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyDeviceToHost);
//	cudaMemcpy(ex_aux_m2m3_c, g_ex_aux_m2m3_c, sizeof(float)*(nx+10)*(nz+10), cudaMemcpyDeviceToHost);
	
	fprintf(stderr,"Data Copy To GPU  ====> OK\n");
}



extern "C" void rtm_gpu_copy_out(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0)//, //(nz, nx, nt)
{
	// data copy back from GPU mem
	cudaMemcpy(Vx0, g_ex_Vx0, sizeof(float)*(nx+10)*(nz+10)*(nt+10),  		cudaMemcpyDeviceToHost);
	cudaMemcpy(Vz0, g_ex_Vz0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), 			cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmaxx0, g_ex_sigmaxx0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmaxz0, g_ex_sigmaxz0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), 	cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmazz0, g_ex_sigmazz0, sizeof(float)*(nx+10)*(nz+10)*(nt+10), 	cudaMemcpyDeviceToHost);
	//cudaMemcpy(sigmazz0, g_sigmazz0,  sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
	fprintf(stderr,"Data Copy To CPU ====> OK\n");

}


extern "C" void rtm_gpu_final()
{

	//release GPU memory space
	cudaFree(&g_ex_Vx0);
	cudaFree(&g_ex_Vz0);
	cudaFree(&g_ex_sigmaxx0);
	cudaFree(&g_ex_sigmazz0);
	cudaFree(&g_ex_sigmaxz0);
	cudaFree(&g_ex_m1_x);
	cudaFree(&g_ex_m1_z);
	cudaFree(&g_ex_aux_m2_c);
	cudaFree(&g_ex_aux_m3_c);
	cudaFree(&g_ex_aux_m2m3_c);
	fprintf(stderr,"GPU Mem Released ====> OK\n");
}


__global__ void rtm_gpu_kernel(int it,int nt, int nz, int nx,
        float * g_ex_Vx0, float * g_ex_Vz0, float * g_ex_sigmaxx0, float * g_ex_sigmazz0, float * g_ex_sigmaxz0, //(nz, nx, nt)
        float * g_ex_m1_x,float * g_ex_m1_z,float * g_ex_aux_m2_c, float * g_ex_aux_m3_c, float * g_ex_aux_m2m3_c)//(nz+10,	nx+10)
{

	float c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0;

	//GPU thread index
	int iz, ix;
	iz = blockIdx.x*blockDim.x + threadIdx.x;
	ix = blockIdx.y*blockDim.y + threadIdx.y;
	//gt = it;
 	

	__shared__ float sh_ex_aux_m2m3_c[(TZ+10)*(TX+10)];
	__shared__ float sh_ex_aux_m2_c[(TZ+10)*(TX+10)];
	__shared__ float sh_ex_aux_m3_c[(TZ+10)*(TX+10)];
	__shared__ float sh_ex_m1_x[(TZ+10)*(TX+10)];
	__shared__ float sh_ex_m1_z[(TZ+10)*(TX+10)];

	//sh_ex_aux_m2m3_c[threadIdx][];

	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y)] = g_ex_aux_m2m3_c[index_ex(iz,ix)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y)] = g_ex_aux_m2_c[index_ex(iz,ix)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y)] = g_ex_aux_m3_c[index_ex(iz,ix)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y)] = g_ex_m1_x[index_ex(iz,ix)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y)] = g_ex_m1_z[index_ex(iz,ix)];

	if(threadIdx.x<5){
	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x-5,threadIdx.y)] = g_ex_aux_m2m3_c[index_ex(iz-5,ix)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x-5,threadIdx.y)] = g_ex_aux_m2_c[index_ex(iz-5,ix)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x-5,threadIdx.y)] = g_ex_aux_m3_c[index_ex(iz-5,ix)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x-5,threadIdx.y)] = g_ex_m1_x[index_ex(iz-5,ix)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x-5,threadIdx.y)] = g_ex_m1_z[index_ex(iz-5,ix)];
	}

	if(threadIdx.x>=TZ-5){
	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x+5,threadIdx.y)] = g_ex_aux_m2m3_c[index_ex(iz+5,ix)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x+5,threadIdx.y)] = g_ex_aux_m2_c[index_ex(iz+5,ix)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x+5,threadIdx.y)] = g_ex_aux_m3_c[index_ex(iz+5,ix)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x+5,threadIdx.y)] = g_ex_m1_x[index_ex(iz+5,ix)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x+5,threadIdx.y)] = g_ex_m1_z[index_ex(iz+5,ix)];
	}
	

	if(threadIdx.y<5){
	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y-5)] = g_ex_aux_m2m3_c[index_ex(iz,ix-5)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y-5)] = g_ex_aux_m2_c[index_ex(iz,ix-5)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y-5)] = g_ex_aux_m3_c[index_ex(iz,ix-5)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y-5)] = g_ex_m1_x[index_ex(iz,ix-5)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y-5)] = g_ex_m1_z[index_ex(iz,ix-5)];
	}


	if(threadIdx.y>=TX-5){
	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y+5)] = g_ex_aux_m2m3_c[index_ex(iz,ix+5)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y+5)] = g_ex_aux_m2_c[index_ex(iz,ix+5)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y+5)] = g_ex_aux_m3_c[index_ex(iz,ix+5)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y+5)] = g_ex_m1_x[index_ex(iz,ix+5)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y+5)] = g_ex_m1_z[index_ex(iz,ix+5)];
	}



	if(threadIdx.x <5 && threadIdx.y <5){
	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x-5,threadIdx.y-5)] = g_ex_aux_m2m3_c[index_ex(iz-5,ix-5)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x-5,threadIdx.y-5)] = g_ex_aux_m2_c[index_ex(iz-5,ix-5)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x-5,threadIdx.y-5)] = g_ex_aux_m3_c[index_ex(iz-5,ix-5)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x-5,threadIdx.y-5)] = g_ex_m1_x[index_ex(iz-5,ix-5)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x-5,threadIdx.y-5)] = g_ex_m1_z[index_ex(iz-5,ix-5)];
	}

	if(threadIdx.x >= 5+TZ && threadIdx.y >= 5+TX){
	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x+5,threadIdx.y+5)] = g_ex_aux_m2m3_c[index_ex(iz+5,ix+5)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x+5,threadIdx.y+5)] = g_ex_aux_m2_c[index_ex(iz+5,ix+5)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x+5,threadIdx.y+5)] = g_ex_aux_m3_c[index_ex(iz+5,ix+5)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x+5,threadIdx.y+5)] = g_ex_m1_x[index_ex(iz+5,ix+5)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x+5,threadIdx.y+5)] = g_ex_m1_z[index_ex(iz+5,ix+5)];
	}


	if(threadIdx.x >= TZ+5 && threadIdx.y <5){
	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x+5,threadIdx.y-5)] = g_ex_aux_m2m3_c[index_ex(iz+5,ix-5)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x+5,threadIdx.y-5)] = g_ex_aux_m2_c[index_ex(iz+5,ix-5)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x+5,threadIdx.y-5)] = g_ex_aux_m3_c[index_ex(iz+5,ix-5)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x+5,threadIdx.y-5)] = g_ex_m1_x[index_ex(iz+5,ix-5)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x+5,threadIdx.y-5)] = g_ex_m1_z[index_ex(iz+5,ix-5)];
	}


	if(threadIdx.x <5 && threadIdx.y >= TX-5){
	sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x-5,threadIdx.y+5)] = g_ex_aux_m2m3_c[index_ex(iz-5,ix+5)];
	sh_ex_aux_m2_c[index_blk_ex(threadIdx.x-5,threadIdx.y+5)] = g_ex_aux_m2_c[index_ex(iz-5,ix+5)];
	sh_ex_aux_m3_c[index_blk_ex(threadIdx.x-5,threadIdx.y+5)] = g_ex_aux_m3_c[index_ex(iz-5,ix+5)];
	sh_ex_m1_x[index_blk_ex(threadIdx.x-5,threadIdx.y+5)] = g_ex_m1_x[index_ex(iz-5,ix+5)];
	sh_ex_m1_z[index_blk_ex(threadIdx.x-5,threadIdx.y+5)] = g_ex_m1_z[index_ex(iz-5,ix+5)];
	}



	__syncthreads();

              g_ex_Vx0[index3d_ex(iz,ix  ,it)] = g_ex_Vx0[index3d_ex(iz,ix  ,it)]	+ g_ex_Vx0[index3d_ex(iz, ix, it+2)]
									+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y-5)]*c1*g_ex_sigmaxx0[index3d_ex(iz,ix-5,it+1)]							
							 		+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y-4)]*c2*g_ex_sigmaxx0[index3d_ex(iz,ix-4,it+1)]		
									+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y-3)]*c3*g_ex_sigmaxx0[index3d_ex(iz,ix-3,it+1)]	
									+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y-2)]*c4*g_ex_sigmaxx0[index3d_ex(iz,ix-2,it+1)]	
									+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y-1)]*c5*g_ex_sigmaxx0[index3d_ex(iz,ix-1,it+1)]	
									- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y)]  *c5*g_ex_sigmaxx0[index3d_ex(iz,ix,it+1)]	
									- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y+1)]*c4*g_ex_sigmaxx0[index3d_ex(iz,ix+1,it+1)]	
									- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y+2)]*c3*g_ex_sigmaxx0[index3d_ex(iz,ix+2,it+1)]	
									- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y+3)]*c2*g_ex_sigmaxx0[index3d_ex(iz,ix+3,it+1)]	
									- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,threadIdx.y+4)]*c1*g_ex_sigmaxx0[index3d_ex(iz,ix+4,it+1)]


									+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y-5)]*c1*g_ex_sigmazz0[index3d_ex(iz,ix-5,it+1)]							
							 		+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y-4)]*c2*g_ex_sigmazz0[index3d_ex(iz,ix-4,it+1)]		
									+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y-3)]*c3*g_ex_sigmazz0[index3d_ex(iz,ix-3,it+1)]	
									+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y-2)]*c4*g_ex_sigmazz0[index3d_ex(iz,ix-2,it+1)]	
									+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y-1)]*c5*g_ex_sigmazz0[index3d_ex(iz,ix-1,it+1)]	
									- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y)]  *c5*g_ex_sigmazz0[index3d_ex(iz,ix,it+1)]	
									- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y+1)]*c4*g_ex_sigmazz0[index3d_ex(iz,ix+1,it+1)]	
									- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y+2)]*c3*g_ex_sigmazz0[index3d_ex(iz,ix+2,it+1)]	
									- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y+3)]*c2*g_ex_sigmazz0[index3d_ex(iz,ix+3,it+1)]	
									- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,threadIdx.y+4)]*c1*g_ex_sigmazz0[index3d_ex(iz,ix+4,it+1)]	
	


									+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x-4,threadIdx.y)]*c1*g_ex_sigmaxz0[index3d_ex(iz-4,ix,it+1)]		
									+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x-3,threadIdx.y)]*c2*g_ex_sigmaxz0[index3d_ex(iz-3,ix,it+1)]	
									+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x-2,threadIdx.y)]*c3*g_ex_sigmaxz0[index3d_ex(iz-2,ix,it+1)]	
									+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x-1,threadIdx.y)]*c4*g_ex_sigmaxz0[index3d_ex(iz-1,ix,it+1)]	
									+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,  threadIdx.y)]  *c5*g_ex_sigmaxz0[index3d_ex(iz,ix,it+1)]	
									- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x+1,threadIdx.y)]*c5*g_ex_sigmaxz0[index3d_ex(iz+1,ix,it+1)]	
									- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x+2,threadIdx.y)]*c4*g_ex_sigmaxz0[index3d_ex(iz+2,ix,it+1)]	
									- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x+3,threadIdx.y)]*c3*g_ex_sigmaxz0[index3d_ex(iz+3,ix,it+1)]	
									- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x+4,threadIdx.y)]*c2*g_ex_sigmaxz0[index3d_ex(iz+4,ix,it+1)]	
									- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x+5,threadIdx.y)]*c1*g_ex_sigmaxz0[index3d_ex(iz+5,ix,it+1)]	;						

 
     __syncthreads();

            g_ex_Vz0[index3d_ex(iz,ix  ,it)] = g_ex_Vz0[index3d_ex(iz,ix,  it)]  	+ g_ex_Vz0[index3d_ex(iz,ix  ,it+2)] 
	     								+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x-5,threadIdx.y)]*c1*g_ex_sigmaxx0[index3d_ex(iz-5,ix,it+1)]							
	     						 		+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x-4,threadIdx.y)]*c2*g_ex_sigmaxx0[index3d_ex(iz-4,ix,it+1)]		
	     								+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x-3,threadIdx.y)]*c3*g_ex_sigmaxx0[index3d_ex(iz-3,ix,it+1)]	
	     								+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x-2,threadIdx.y)]*c4*g_ex_sigmaxx0[index3d_ex(iz-2,ix,it+1)]	
	     								+ sh_ex_aux_m2_c[index_blk_ex(threadIdx.x-1,threadIdx.y)]*c5*g_ex_sigmaxx0[index3d_ex(iz-1,ix,it+1)]	
	     								- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x,  threadIdx.y)]  *c5*g_ex_sigmaxx0[index3d_ex(iz,ix,it+1)]	
	     								- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x+1,threadIdx.y)]*c4*g_ex_sigmaxx0[index3d_ex(iz+1,ix,it+1)]	
	     								- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x+2,threadIdx.y)]*c3*g_ex_sigmaxx0[index3d_ex(iz+2,ix,it+1)]	
	     								- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x+3,threadIdx.y)]*c2*g_ex_sigmaxx0[index3d_ex(iz+3,ix,it+1)]	
	     								- sh_ex_aux_m2_c[index_blk_ex(threadIdx.x+4,threadIdx.y)]*c1*g_ex_sigmaxx0[index3d_ex(iz+4,ix,it+1)]	
	     
	
	             							+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x-5,threadIdx.y)]*c1*g_ex_sigmazz0[index3d_ex(iz-5,ix,it+1)]							
	     						 		+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x-4,threadIdx.y)]*c2*g_ex_sigmazz0[index3d_ex(iz-4,ix,it+1)]		
	     								+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x-3,threadIdx.y)]*c3*g_ex_sigmazz0[index3d_ex(iz-3,ix,it+1)]	
	     								+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x-2,threadIdx.y)]*c4*g_ex_sigmazz0[index3d_ex(iz-2,ix,it+1)]	
	     								+ sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x-1,threadIdx.y)]*c5*g_ex_sigmazz0[index3d_ex(iz-1,ix,it+1)]	
	     								- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x,  threadIdx.y)]  *c5*g_ex_sigmazz0[index3d_ex(iz,ix,it+1)]	
	     								- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x+1,threadIdx.y)]*c4*g_ex_sigmazz0[index3d_ex(iz+1,ix,it+1)]	
	     								- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x+2,threadIdx.y)]*c3*g_ex_sigmazz0[index3d_ex(iz+2,ix,it+1)]	
	     								- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x+3,threadIdx.y)]*c2*g_ex_sigmazz0[index3d_ex(iz+3,ix,it+1)]	
	     								- sh_ex_aux_m2m3_c[index_blk_ex(threadIdx.x+4,threadIdx.y)]*c1*g_ex_sigmazz0[index3d_ex(iz+4,ix,it+1)]	
	     
	     								+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y-4)]*c1*g_ex_sigmaxz0[index3d_ex(iz,ix-4,it+1)]		
	     								+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y-3)]*c2*g_ex_sigmaxz0[index3d_ex(iz,ix-3,it+1)]	
	     								+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y-2)]*c3*g_ex_sigmaxz0[index3d_ex(iz,ix-2,it+1)]	
	     								+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y-1)]*c4*g_ex_sigmaxz0[index3d_ex(iz,ix-1,it+1)]	
	     								+ sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y)]  *c5*g_ex_sigmaxz0[index3d_ex(iz,ix,it+1)]	
	     								- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y+1)]*c5*g_ex_sigmaxz0[index3d_ex(iz,ix+1,it+1)]	
	     								- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y+2)]*c4*g_ex_sigmaxz0[index3d_ex(iz,ix+2,it+1)]	
	     								- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y+3)]*c3*g_ex_sigmaxz0[index3d_ex(iz,ix+3,it+1)]	
	     								- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y+4)]*c2*g_ex_sigmaxz0[index3d_ex(iz,ix+4,it+1)]	
	     								- sh_ex_aux_m3_c[index_blk_ex(threadIdx.x,threadIdx.y+5)]*c1*g_ex_sigmaxz0[index3d_ex(iz,ix+5,it+1)]	;							
	


              g_ex_sigmaxx0[index3d_ex(iz,ix  ,it)] = g_ex_sigmaxx0[index3d_ex(iz,ix  ,it)]	+ g_ex_sigmaxx0[index3d_ex(iz,ix  ,it+2)] 
        									+ sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y-4)]*c1*g_ex_Vx0[index3d_ex(iz,ix-4,it+1)]		
        									+ sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y-3)]*c2*g_ex_Vx0[index3d_ex(iz,ix-3,it+1)]	
        									+ sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y-2)]*c3*g_ex_Vx0[index3d_ex(iz,ix-2,it+1)]	
        									+ sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y-1)]*c4*g_ex_Vx0[index3d_ex(iz,ix-1,it+1)]	
        									+ sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y)]  *c5*g_ex_Vx0[index3d_ex(iz,ix,it+1)]	
        									- sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y+1)]*c5*g_ex_Vx0[index3d_ex(iz,ix+1,it+1)]	
        									- sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y+2)]*c4*g_ex_Vx0[index3d_ex(iz,ix+2,it+1)]	
        									- sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y+3)]*c3*g_ex_Vx0[index3d_ex(iz,ix+3,it+1)]	
        									- sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y+4)]*c2*g_ex_Vx0[index3d_ex(iz,ix+4,it+1)]	
        									- sh_ex_m1_x[index_blk_ex(threadIdx.x,threadIdx.y+5)]*c1*g_ex_Vx0[index3d_ex(iz,ix+5,it+1)]	;						
 
    __syncthreads();
             g_ex_sigmazz0[index3d_ex(iz,ix  ,it)] = g_ex_sigmazz0[index3d_ex(iz,ix  ,it)]	+ g_ex_sigmazz0[index3d_ex(iz,ix  ,it+2)] 
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x-4,threadIdx.y)]*c1*g_ex_Vz0[index3d_ex(iz-4,ix,it+1)]		
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x-3,threadIdx.y)]*c2*g_ex_Vz0[index3d_ex(iz-3,ix,it+1)]	
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x-2,threadIdx.y)]*c3*g_ex_Vz0[index3d_ex(iz-2,ix,it+1)]	
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x-1,threadIdx.y)]*c4*g_ex_Vz0[index3d_ex(iz-1,ix,it+1)]	
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x,  threadIdx.y)]  *c5*g_ex_Vz0[index3d_ex(iz,ix,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x+1,threadIdx.y)]*c5*g_ex_Vz0[index3d_ex(iz+1,ix,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x+2,threadIdx.y)]*c4*g_ex_Vz0[index3d_ex(iz+2,ix,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x+3,threadIdx.y)]*c3*g_ex_Vz0[index3d_ex(iz+3,ix,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x+4,threadIdx.y)]*c2*g_ex_Vz0[index3d_ex(iz+4,ix,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x+5,threadIdx.y)]*c1*g_ex_Vz0[index3d_ex(iz+5,ix,it+1)]	;						
     __syncthreads();
     g_ex_sigmaxz0[index3d_ex(iz,ix  ,it)] = g_ex_sigmaxz0[index3d_ex(iz,ix  ,it)]	+ g_ex_sigmaxz0[index3d_ex(iz,ix  ,it+2)]	 
										+ sh_ex_m1_x[index_blk_ex(threadIdx.x-5,threadIdx.y)]*c1*g_ex_Vx0[index3d_ex(iz-5,ix,it+1)]							
							 			+ sh_ex_m1_x[index_blk_ex(threadIdx.x-4,threadIdx.y)]*c2*g_ex_Vx0[index3d_ex(iz-4,ix,it+1)]		
										+ sh_ex_m1_x[index_blk_ex(threadIdx.x-3,threadIdx.y)]*c3*g_ex_Vx0[index3d_ex(iz-3,ix,it+1)]	
										+ sh_ex_m1_x[index_blk_ex(threadIdx.x-2,threadIdx.y)]*c4*g_ex_Vx0[index3d_ex(iz-2,ix,it+1)]	
										+ sh_ex_m1_x[index_blk_ex(threadIdx.x-1,threadIdx.y)]*c5*g_ex_Vx0[index3d_ex(iz-1,ix,it+1)]	
										- sh_ex_m1_x[index_blk_ex(threadIdx.x,  threadIdx.y)]  *c5*g_ex_Vx0[index3d_ex(iz,ix,it+1)]	
										- sh_ex_m1_x[index_blk_ex(threadIdx.x+1,threadIdx.y)]*c4*g_ex_Vx0[index3d_ex(iz+1,ix,it+1)]	
										- sh_ex_m1_x[index_blk_ex(threadIdx.x+2,threadIdx.y)]*c3*g_ex_Vx0[index3d_ex(iz+2,ix,it+1)]	
										- sh_ex_m1_x[index_blk_ex(threadIdx.x+3,threadIdx.y)]*c2*g_ex_Vx0[index3d_ex(iz+3,ix,it+1)]	
										- sh_ex_m1_x[index_blk_ex(threadIdx.x+4,threadIdx.y)]*c1*g_ex_Vx0[index3d_ex(iz+4,ix,it+1)]	//;
	
        
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y-5)]*c1*g_ex_Vz0[index3d_ex(iz,ix-5,it+1)]							
							 			+ sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y-4)]*c2*g_ex_Vz0[index3d_ex(iz,ix-4,it+1)]		
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y-3)]*c3*g_ex_Vz0[index3d_ex(iz,ix-3,it+1)]	
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y-2)]*c4*g_ex_Vz0[index3d_ex(iz,ix-2,it+1)]	
										+ sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y-1)]*c5*g_ex_Vz0[index3d_ex(iz,ix-1,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y)]  *c5*g_ex_Vz0[index3d_ex(iz,ix,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y+1)]*c4*g_ex_Vz0[index3d_ex(iz,ix+1,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y+2)]*c3*g_ex_Vz0[index3d_ex(iz,ix+2,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y+3)]*c2*g_ex_Vz0[index3d_ex(iz,ix+3,it+1)]	
										- sh_ex_m1_z[index_blk_ex(threadIdx.x,threadIdx.y+4)]*c1*g_ex_Vz0[index3d_ex(iz,ix+4,it+1)]	;
		
	__syncthreads();


	}



extern "C" void rtm_gpu_func(int it_max, int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * ex_m1_x,float * ex_m1_z,float * ex_aux_m2_c, float * ex_aux_m3_c, float * ex_aux_m2m3_c,//)//(nz+10,nx+10)
	float * debug)
{	
     	cudaError_t err;
	cudaEvent_t start, stop;
	float elapsedTime = 0.0f;
	int g_it;

	//time record

	//data copy in 
     	rtm_gpu_copy_in(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, ex_m1_x, ex_m1_z, ex_aux_m2_c, ex_aux_m3_c, ex_aux_m2m3_c);
	
	err = cudaGetLastError();
	if(cudaSuccess != err){
		fprintf(stderr, "Cuda error1: %s.\n", cudaGetErrorString(err));
	}	
	
	//RTM computing
	fprintf(stderr,"GPU Computing from TS=%d ... ...(NZ=%d, NX=%d, TZ=%d, TX=%d)\n", it_max, nz, nx, TZ, TX);
	for(g_it = it_max; g_it>=0; g_it--){

	dim3 dimGrid(nz/TZ, nx/TX);
	dim3 dimBlock(TZ, TX);

	rtm_gpu_kernel<<<dimGrid, dimBlock>>>(g_it,nt, nz, nx, g_ex_Vx0, g_ex_Vz0, g_ex_sigmaxx0, g_ex_sigmazz0, g_ex_sigmaxz0, g_ex_m1_x, g_ex_m1_z, g_ex_aux_m2_c, g_ex_aux_m3_c, g_ex_aux_m2m3_c);
	cudaThreadSynchronize();

	err = cudaGetLastError();
	if(cudaSuccess != err){
		fprintf(stderr, "Cuda error2: %s.\n", cudaGetErrorString(err));
		}
	}

	//data copy out
	rtm_gpu_copy_out(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0);	

	err = cudaGetLastError();
	if(cudaSuccess != err){
		fprintf(stderr, "Cuda error3: %s.\n", cudaGetErrorString(err));
	}	


	//cudaEventRecord(stop, 0);
	//cudaEventSynchronize(stop);
	//cudaEventElapsedTime(&elapsedTime, start, stop);
	//fprintf(stderr,"GPU Computational Elapsed Time: %.4f\n",elapsedTime);
}

