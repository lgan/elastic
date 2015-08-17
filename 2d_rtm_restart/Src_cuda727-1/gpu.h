//parameter settings 

#define TZ (30)	//z direction >>fast 
#define TX (30)	//x direction >>middle
#define TT (10)	//t direction >>slow
	

#define index(a,b) ((a)+(b)*nz)     //an index transpose happens here
#define indextx(a,b) ((a)+(b)*nt)     //an index transpose happens here
#define index3d(a,b,c) ( (c)*nx*nz + (b)*nz + (a))     //an index transpose happens here

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

