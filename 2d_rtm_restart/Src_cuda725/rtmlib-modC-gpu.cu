#include <stdio.h>
#include "gpu.h"
__device__ float *g_Vx0;
__device__ float *g_Vz0;
__device__ float *g_sigmaxx0; 
__device__ float *g_sigmazz0;
__device__ float *g_sigmaxz0;
__device__ float *g_m1_x;
__device__ float *g_m1_z;
__device__ float *g_aux_m2_c; 
__device__ float *g_aux_m3_c; 
__device__ float *g_aux_m2m3_c; 


//void setup_cuda(int ngpus, int argc, char **argv){
	//insert from Bob' Born
//	;
//}
//void process_error( const cudaError_t &error, char *string=0, bool verbose=false ){
	//insert from Bob's Born
//	;
//}
extern "C" void rtm_gpu_init(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)//(nz,	nx)
{
	//set cuda devices and put all data onto gpu memory
	
	cudaError_t cuda_ret;
     	cudaError_t err;

	//Set Device 
	fprintf(stderr,"GPU init. \n");
    	cuda_ret = cudaSetDevice(0);
	if(cuda_ret != cudaSuccess){
		fprintf(stderr, "Failed to Set The cuda Device !\n");
	}
	else{
		fprintf(stderr, "GPU Device Set OK\n");
	}

	// data init
	cudaMalloc(&g_Vx0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_Vz0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_sigmaxx0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_sigmazz0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_sigmaxz0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_m1_x,sizeof(float)*nx*nz);
	cudaMalloc(&g_m1_z,sizeof(float)*nx*nz);
	cudaMalloc(&g_aux_m2_c,sizeof(float)*nx*nz);
	cudaMalloc(&g_aux_m3_c,sizeof(float)*nx*nz);
	cudaMalloc(&g_aux_m2m3_c,sizeof(float)*nx*nz);
	fprintf(stderr,"GPU Data Init OK\n");

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
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)//(nz,	nx)
{
	// data copy

	cudaMemcpy(g_Vx0, Vx0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
	cudaMemcpy(g_Vz0, Vz0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
	cudaMemcpy(g_sigmaxx0, sigmaxx0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
	cudaMemcpy(g_sigmaxz0, sigmaxz0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
	cudaMemcpy(g_sigmazz0, sigmazz0, sizeof(float)*nx*nz*nt, cudaMemcpyHostToDevice);
	cudaMemcpy(g_m1_x, m1_x, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
	cudaMemcpy(g_m1_z, m1_z, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
	cudaMemcpy(g_aux_m2_c, aux_m2_c, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
	cudaMemcpy(g_aux_m3_c, aux_m3_c, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
	cudaMemcpy(g_aux_m2m3_c, aux_m2m3_c, sizeof(float)*nx*nz, cudaMemcpyHostToDevice);
	fprintf(stderr,"Data Copy To GPU OK\n");
}


extern "C" void rtm_gpu_copy_out(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)//(nz,	nx)
{
	// data copy back from GPU mem
	cudaMemcpy(Vx0, g_Vx0, sizeof(float)*nx*nz*nt,  		cudaMemcpyDeviceToHost);
	cudaMemcpy( Vz0, g_Vz0,sizeof(float)*nx*nz*nt, 			cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmaxx0, g_sigmaxx0, sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmaxz0, g_sigmaxz0, sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmazz0, g_sigmazz0,  sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
	cudaMemcpy(m1_x, g_m1_x,  sizeof(float)*nx*nz, 			cudaMemcpyDeviceToHost);
	cudaMemcpy(m1_z, g_m1_z,  sizeof(float)*nx*nz, 			cudaMemcpyDeviceToHost);
	cudaMemcpy(aux_m2_c, g_aux_m2_c,  sizeof(float)*nx*nz, 		cudaMemcpyDeviceToHost);
	cudaMemcpy(aux_m3_c, g_aux_m3_c,  sizeof(float)*nx*nz, 		cudaMemcpyDeviceToHost);
	cudaMemcpy(aux_m2m3_c, g_aux_m2m3_c,  sizeof(float)*nx*nz, 	cudaMemcpyDeviceToHost);
	fprintf(stderr,"Data Copy To CPU OK\n");

}


extern "C" void rtm_gpu_final(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)//(nz,	nx)
{
	// data copy back from GPU mem
//	cudaMemcpy(Vx0, g_Vx0, sizeof(float)*nx*nz*nt,  		cudaMemcpyDeviceToHost);
//	cudaMemcpy( Vz0, g_Vz0,sizeof(float)*nx*nz*nt, 			cudaMemcpyDeviceToHost);
//	cudaMemcpy(sigmaxx0, g_sigmaxx0, sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
//	cudaMemcpy(sigmaxz0, g_sigmaxz0, sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
//	cudaMemcpy(sigmazz0, g_sigmazz0,  sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
//	cudaMemcpy(m1_x, g_m1_x,  sizeof(float)*nx*nz, 			cudaMemcpyDeviceToHost);
//	cudaMemcpy(m1_z, g_m1_z,  sizeof(float)*nx*nz, 			cudaMemcpyDeviceToHost);
//	cudaMemcpy(aux_m2_c, g_aux_m2_c,  sizeof(float)*nx*nz, 		cudaMemcpyDeviceToHost);
//	cudaMemcpy(aux_m3_c, g_aux_m3_c,  sizeof(float)*nx*nz, 		cudaMemcpyDeviceToHost);
//	cudaMemcpy(aux_m2m3_c, g_aux_m2m3_c,  sizeof(float)*nx*nz, 	cudaMemcpyDeviceToHost);
//	fprintf(stderr,"Data Copy To CPU OK\n");


	cudaFree(&g_Vx0);
	cudaFree(&g_Vz0);
	cudaFree(&g_sigmaxx0);
	cudaFree(&g_sigmazz0);
	cudaFree(&g_sigmaxz0);
	cudaFree(&g_m1_x);
	cudaFree(&g_m1_z);
	cudaFree(&g_aux_m2_c);
	cudaFree(&g_aux_m3_c);
	cudaFree(&g_aux_m2m3_c);
	fprintf(stderr,"GPU Mem Released OK\n");
}


__global__ void rtm_gpu_kernel(int it,int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0) //(nz, nx, nt)
        //float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)//(nz,	nx)
{
	//GPU thread index
	int gz, gx, gt;
	gz = blockIdx.x*blockDim.x + threadIdx.x;
	gx = blockIdx.y*blockDim.y + threadIdx.y;
	gt = it;
//	gt = blockIdx.z*blockDim.y + threadIdx.z;
 
       Vx0[index3d(gz, gx, gt)] = Vx0[index3d(gz, gx, gt)] + Vx0[index3d(gz, gx, gt+2)];
       Vz0[index3d(gz, gx, gt)] = Vz0[index3d(gz, gx, gt)] + Vz0[index3d(gz, gx, gt+2)];
       sigmaxx0[index3d(gz, gx, gt)] = sigmaxx0[index3d(gz, gx, gt)] + sigmaxx0[index3d(gz, gx, gt+2)];
       sigmazz0[index3d(gz, gx, gt)] = sigmazz0[index3d(gz, gx, gt)] + sigmazz0[index3d(gz, gx, gt+2)];
       sigmaxz0[index3d(gz, gx, gt)] = sigmaxz0[index3d(gz, gx, gt)] + sigmaxz0[index3d(gz, gx, gt+2)];
}


extern "C" void rtm_gpu_func(int it, int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)//(nz,	nx)
{	
     	cudaError_t err;
	cudaEvent_t start, stop;
	float elapsedTime = 0.0f;

	//time record
	
	dim3 dimGrid(nz/TZ, nx/TX);
	dim3 dimBlock(TZ, TX);

	//RTM kernel 
	fprintf(stderr,"GPU Computing...(NZ=%d, NX=%d, TZ=%d, TX=%d)\n", nz, nx, TZ, TX);

	//cudaEventRecord(start, 0);
	rtm_gpu_kernel<<<dimGrid, dimBlock>>>(it,nt, nz, nx, g_Vx0, g_Vz0, g_sigmaxx0, g_sigmazz0, g_sigmaxz0);
	cudaThreadSynchronize();

	err = cudaGetLastError();
	if(cudaSuccess != err)
		fprintf(stderr, "Cuda error: %s.\n", cudaGetErrorString(err));


	//cudaEventRecord(stop, 0);
	//cudaEventSynchronize(stop);
	//cudaEventElapsedTime(&elapsedTime, start, stop);
	//fprintf(stderr,"GPU Computational Elapsed Time: %.4f\n",elapsedTime);



	
}

