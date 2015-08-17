#include <stdio.h>
#include "gpu.h"
 float *g_Vx0;
// float *g_Vx0_out;
 float *g_Vz0;
 float *g_sigmaxx0; 
 float *g_sigmazz0;
 float *g_sigmaxz0;
 float *g_m1_x;
 float *g_m1_z;
 float *g_aux_m2_c; 
 float *g_aux_m3_c; 
 float *g_aux_m2m3_c; 


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
    	cuda_ret = cudaSetDevice(0);
	if(cuda_ret != cudaSuccess){
		fprintf(stderr, "Failed to Set The cuda Device !\n");
	}
	else{
		fprintf(stderr, "GPU Device Set ====> OK\n");
	}

	// data init
	cudaMalloc(&g_Vx0,sizeof(float)*nx*nz*nt);
//	cudaMalloc(&g_Vx0_out,sizeof(float)*nx*nz*nt);
//	cudaMemset(g_Vx0_out, 0, sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_Vz0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_sigmaxx0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_sigmazz0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_sigmaxz0,sizeof(float)*nx*nz*nt);
	cudaMalloc(&g_m1_x,sizeof(float)*nx*nz);
	cudaMalloc(&g_m1_z,sizeof(float)*nx*nz);
	cudaMalloc(&g_aux_m2_c,sizeof(float)*nx*nz);
	cudaMalloc(&g_aux_m3_c,sizeof(float)*nx*nz);
	cudaMalloc(&g_aux_m2m3_c,sizeof(float)*nx*nz);
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
	fprintf(stderr,"Data Copy To GPU  ====> OK\n");
}


extern "C" void rtm_gpu_copy_out(int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0)//, //(nz, nx, nt)
        //float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)//(nz,	nx)
{
	// data copy back from GPU mem
	cudaMemcpy(Vx0, g_Vx0, sizeof(float)*nx*nz*nt,  		cudaMemcpyDeviceToHost);
	cudaMemcpy(Vz0, g_Vz0,sizeof(float)*nx*nz*nt, 			cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmaxx0, g_sigmaxx0, sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmaxz0, g_sigmaxz0, sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
	cudaMemcpy(sigmazz0, g_sigmazz0,  sizeof(float)*nx*nz*nt, 	cudaMemcpyDeviceToHost);
//	cudaMemcpy(m1_x, g_m1_x,  sizeof(float)*nx*nz, 			cudaMemcpyDeviceToHost);
//	cudaMemcpy(m1_z, g_m1_z,  sizeof(float)*nx*nz, 			cudaMemcpyDeviceToHost);
//	cudaMemcpy(aux_m2_c, g_aux_m2_c,  sizeof(float)*nx*nz, 		cudaMemcpyDeviceToHost);
//	cudaMemcpy(aux_m3_c, g_aux_m3_c,  sizeof(float)*nx*nz, 		cudaMemcpyDeviceToHost);
//	cudaMemcpy(aux_m2m3_c, g_aux_m2m3_c,  sizeof(float)*nx*nz, 	cudaMemcpyDeviceToHost);
	fprintf(stderr,"Data Copy To CPU ====> OK\n");

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
///	cudaFree(&g_Vx0_out);
	cudaFree(&g_Vz0);
	cudaFree(&g_sigmaxx0);
	cudaFree(&g_sigmazz0);
	cudaFree(&g_sigmaxz0);
	cudaFree(&g_m1_x);
	cudaFree(&g_m1_z);
	cudaFree(&g_aux_m2_c);
	cudaFree(&g_aux_m3_c);
	cudaFree(&g_aux_m2m3_c);
	fprintf(stderr,"GPU Mem Released ====> OK\n");
}


__global__ void rtm_gpu_kernel(int it,int nt, int nz, int nx,
        float * g_Vx0, float * g_Vz0, float * g_sigmaxx0, float * g_sigmazz0, float * g_sigmaxz0, //(nz, nx, nt)
        float * g_m1_x,float * g_m1_z,float * g_aux_m2_c, float * g_aux_m3_c, float * g_aux_m2m3_c)//(nz,	nx)
{

	float c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0;

	//GPU thread index
	int iz, ix;
	iz = blockIdx.x*blockDim.x + threadIdx.x;
	ix = blockIdx.y*blockDim.y + threadIdx.y;
	//gt = it;
//	gt = blockIdx.z*blockDim.y + threadIdx.z;
 	
//       g_Vx0[index3d(gz, gx, gt)] = g_Vx0[index3d(gz, gx, gt)] + g_Vx0[index3d(gz, gx, gt+2)];
//       g_Vz0[index3d(gz, gx, gt)] = g_Vz0[index3d(gz, gx, gt)] + g_Vz0[index3d(gz, gx, gt+2)];
//       g_sigmaxx0[index3d(gz, gx, gt)] = g_sigmaxx0[index3d(gz, gx, gt)] + g_sigmaxx0[index3d(gz, gx, gt+2)];
//       g_sigmazz0[index3d(gz, gx, gt)] = g_sigmazz0[index3d(gz, gx, gt)] + g_sigmazz0[index3d(gz, gx, gt+2)];
//       g_sigmaxz0[index3d(gz, gx, gt)] = g_sigmaxz0[index3d(gz, gx, gt)] + g_sigmaxz0[index3d(gz, gx, gt+2)];
	if(ix>=9 && ix<(nx-9)  &&  iz>=4 && iz<(nz-5)){

              g_Vx0[index3d(iz,ix  ,it)] = g_Vx0[index3d(iz,ix  ,it)];//	+ g_aux_m2m3_c[index(iz,ix-5)]*c1*g_sigmaxx0[index3d(iz,ix-5,it+1)];
							 	      //  + g_aux_m2m3_c[index(iz,ix-4)]*c2*g_sigmaxx0[index3d(iz,ix-4,it+1)];		
								      // + g_aux_m2m3_c[index(iz,ix-3)]*c3*g_sigmaxx0[index3d(iz,ix-3,it+1)]	
								      // + g_aux_m2m3_c[index(iz,ix-2)]*c4*g_sigmaxx0[index3d(iz,ix-2,it+1)]	
								      // + g_aux_m2m3_c[index(iz,ix-1)]*c5*g_sigmaxx0[index3d(iz,ix-1,it+1)]	
								      // - g_aux_m2m3_c[index(iz,ix)]  *c5*g_sigmaxx0[index3d(iz,ix,it+1)]	
								      // - g_aux_m2m3_c[index(iz,ix+1)]*c4*g_sigmaxx0[index3d(iz,ix+1,it+1)]	
								      // - g_aux_m2m3_c[index(iz,ix+2)]*c3*g_sigmaxx0[index3d(iz,ix+2,it+1)]	
								      // - g_aux_m2m3_c[index(iz,ix+3)]*c2*g_sigmaxx0[index3d(iz,ix+3,it+1)]	
								      // - g_aux_m2m3_c[index(iz,ix+4)]*c1*g_sigmaxx0[index3d(iz,ix+4,it+1)]	;
	}


	__syncthreads();

//	       g_Vx0[index3d(iz,ix+5,it)] = g_Vx0[index3d(iz,ix+5,it)] + g_aux_m2m3_c[index(iz,ix)]*c1*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix+4,it)] = g_Vx0[index3d(iz,ix+4,it)] + g_aux_m2m3_c[index(iz,ix)]*c2*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix+3,it)] = g_Vx0[index3d(iz,ix+3,it)] + g_aux_m2m3_c[index(iz,ix)]*c3*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix+2,it)] = g_Vx0[index3d(iz,ix+2,it)] + g_aux_m2m3_c[index(iz,ix)]*c4*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix+1,it)] = g_Vx0[index3d(iz,ix+1,it)] + g_aux_m2m3_c[index(iz,ix)]*c5*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix  ,it)] = g_Vx0[index3d(iz,ix  ,it)] - g_aux_m2m3_c[index(iz,ix)]*c5*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix-1,it)] = g_Vx0[index3d(iz,ix-1,it)] - g_aux_m2m3_c[index(iz,ix)]*c4*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix-2,it)] = g_Vx0[index3d(iz,ix-2,it)] - g_aux_m2m3_c[index(iz,ix)]*c3*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix-3,it)] = g_Vx0[index3d(iz,ix-3,it)] - g_aux_m2m3_c[index(iz,ix)]*c2*g_sigmaxx0[index3d(iz,ix,it+1)];
//	       g_Vx0[index3d(iz,ix-4,it)] = g_Vx0[index3d(iz,ix-4,it)] - g_aux_m2m3_c[index(iz,ix)]*c1*g_sigmaxx0[index3d(iz,ix,it+1)];
//	}
}


extern "C" void rtm_gpu_func(int it, int nt, int nz, int nx, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)//(nz,	nx)
{	
     	cudaError_t err;
	cudaEvent_t start, stop;
	float elapsedTime = 0.0f;

	//time record
	//rtm_gpu_init(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	err = cudaGetLastError();
	if(cudaSuccess != err){
		fprintf(stderr, "Cuda error0: %s.\n", cudaGetErrorString(err));
	}	

     	rtm_gpu_copy_in(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);

	err = cudaGetLastError();
	if(cudaSuccess != err){
		fprintf(stderr, "Cuda error1: %s.\n", cudaGetErrorString(err));
	}	
	
	dim3 dimGrid(nz/TZ, nx/TX);
	dim3 dimBlock(TZ, TX);
	//RTM kernel 
	//cudaEventRecord(start, 0);
	fprintf(stderr,"GPU Computing...(NZ=%d, NX=%d, TZ=%d, TX=%d)\n", nz, nx, TZ, TX);
	rtm_gpu_kernel<<<dimGrid, dimBlock>>>(it,nt, nz, nx, g_Vx0, g_Vz0, g_sigmaxx0, g_sigmazz0, g_sigmaxz0, g_m1_x, g_m1_z, g_aux_m2_c, g_aux_m3_c, g_aux_m2m3_c);
	cudaThreadSynchronize();

	err = cudaGetLastError();
	if(cudaSuccess != err){
		fprintf(stderr, "Cuda error2: %s.\n", cudaGetErrorString(err));
	}	

	rtm_gpu_copy_out(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0);//, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	cudaThreadSynchronize();
	err = cudaGetLastError();
	if(cudaSuccess != err){
		fprintf(stderr, "Cuda error3: %s.\n", cudaGetErrorString(err));
	}	

//	rtm_gpu_final(nt, nz, nx, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c);
	err = cudaGetLastError();
	if(cudaSuccess != err){
		fprintf(stderr, "Cuda error4: %s.\n", cudaGetErrorString(err));
	}	
	

	//cudaEventRecord(stop, 0);
	//cudaEventSynchronize(stop);
	//cudaEventElapsedTime(&elapsedTime, start, stop);
	//fprintf(stderr,"GPU Computational Elapsed Time: %.4f\n",elapsedTime);
}

