#include <stdio.h>


//void setup_cuda(int ngpus, int argc, char **argv){
	//insert from Bob' Born
//	;
//}
//void process_error( const cudaError_t &error, char *string=0, bool verbose=false ){
	//insert from Bob's Born
//	;
//}


extern "C" void rtm_gpu_init(int nt, int nz, int nx, int zrec, 
        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
//        float * Vx,  float * Vz,  float * sigmaxx,  float * sigmazz,  float * sigmaxz, //(nt, nx)
        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c)
{
	//set cuda devices and put all data onto gpu memory
	
	cudaError_t cuda_ret;
     	cudaError_t err;

	//Set Device 
	fprintf(stderr,"GPU init. \n");
    	cuda_ret = cudaSetDevice(0);
	if(cuda_ret != cudaSuccess){
		fprintf(stderr, "Failed to set the cuda device !\n");
	}
	else{
		fprintf(stderr, "cuda device set OK\n");
	}

	// init data
	//cudaMalloc(&g_,sizeof()*nx*nz*nt);

}

