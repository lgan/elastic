//parameter settings 

#define TZ (20)	//z direction >>fast 
#define TX (20)	//x direction >>middle
#define TT (10)	//t direction >>slow
	

#define index(a,b) ((a)+(b)*nz)     //an index transpose happens here
#define index_ex(a,b) ((a+5)+(b+5)*(nz+10))     //an index transpose happens here
#define index_ex_ori(a,b) ((a)+(b)*(nz+10))     //an index transpose happens here

#define index3d(a,b,c) ( (c)*nx*nz + (b)*nz + (a))     //an index transpose happens here
#define index3d_ex(a,b,c) ( (c+5)*(nx+10)*(nz+10) + (b+5)*(nz+10) + (a+5))     //an index transpose happens here, the expanded array
#define index3d_ex_ori(a,b,c) ( (c)*(nx+10)*(nz+10) + (b)*(nz+10) + (a))  //an index transpose happens here, the expanded array of offset


#define GPU_start_step (nt-3)

#define gz iz
#define gx ix
#define gt it

//extern void rtm_gpu_init(int nt, int nz, int nx, 
//        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
//        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c);//int nt, int nz, int nx, int zrec, 

//extern void rtm_gpu_final(int nt, int nz, int nx, 
//        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
//        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c);//int nt, int nz, int nx, int zrec, 

//extern void rtm_gpu_func(int it, int nt, int nz, int nx, 
//        float * Vx0, float * Vz0, float * sigmaxx0, float * sigmazz0, float * sigmaxz0, //(nz, nx, nt)
//        float * m1_x,float * m1_z,float * aux_m2_c, float * aux_m3_c, float * aux_m2m3_c);//(nz,	nx)

