//parameter settings 

#define DEBUG  // if defined, CPU will compute and compare its results with GPU

#define TZ (10)	//z direction >>fast 
#define TX (10)	//x direction >>middle
#define TY (5)	//t direction >>slow
#define NY (200)	//t direction >>slow
#define NX (200)	//t direction >>slow
#define NZ (200)	//t direction >>slow
#define GPU_start_step (-1)//(nt-10)
#define Steps_write_back (1) // Every () steps, data write back from GPU to CPU for correlation


//#define TT (10)	//t direction >>slow
	
//Index for 2D RTM
#define index(a,b) ((a)+(b)*nz)     //an index transpose happens here
#define index_ex(a,b) ((a+5)+(b+5)*(nz+10))     //an index transpose happens here
#define index_ex_ori(a,b) ((a)+(b)*(nz+10))     //an index transpose happens here

#define index3d(a,b,c) ( (c)*nx*nz + (b)*nz + (a))     //an index transpose happens here
#define index3d_ex(a,b,c) ( (c+5)*(nx+10)*(nz+10) + (b+5)*(nz+10) + (a+5))     //an index transpose happens here, the expanded array
#define index3d_ex_ori(a,b,c) ( (c)*(nx+10)*(nz+10) + (b)*(nz+10) + (a))  //an index transpose happens here, the expanded array of offset

#define index_blk_ex(a,b) ((b+5)*(TZ+10) + (a+5))
#define index_blk(a,b) ((b)*(TZ+10) + (a))

//Index for 3D RTM
#define n3d_index_ex_ori(a,b,c) ((c)*(nx+10)*(nz+10) + (b)*(nz+10) +(a))
#define n3d_index_ex(a,b,c) ((c+5)*(nx+10)*(nz+10) + (b+5)*(nz+10) +(a+5))



#define gz iz
#define gx ix
#define gt it


extern "C" void rtm_gpu_init(int ny, int nz, int nx);


extern "C" void rtm_gpu_final();


extern "C" void rtm_gpu_func(int ny, int nz, int nx, 
        float *ex_Vy0_in,  float * ex_Vx0_in, float * ex_Vz0_in, float * ex_sigmayy0_in, float *ex_sigmaxx0_in, float * ex_sigmazz0_in, float * ex_sigmaxy0_in, float * ex_sigmaxz0_in, float * ex_sigmayz0_in,//(nz, nx, nt)
        float *ex_Vy0_in1,  float * ex_Vx0_in1, float * ex_Vz0_in1, float * ex_sigmayy0_in1, float *ex_sigmaxx0_in1, float * ex_sigmazz0_in1, float * ex_sigmaxy0_in1, float * ex_sigmaxz0_in1, float * ex_sigmayz0_in1,//(nz, nx, nt)
        float *ex_Vy0_out,  float * ex_Vx0_out, float * ex_Vz0_out, float * ex_sigmayy0_out, float *ex_sigmaxx0_out, float * ex_sigmazz0_out, float * ex_sigmaxy0_out, float * ex_sigmaxz0_out, float * ex_sigmayz0_out,//(nz, nx, nt)
        float * ex_m1_y, float * ex_m1_x,float * ex_m1_z,float * ex_m2, float * ex_m3, float * ex_m2m3,//)//(nz+10,nx+10)
	float * debug, float* gpu_kernel_time_ps);

