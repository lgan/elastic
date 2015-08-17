! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Jan, 22th, 2015

! Forward Born modeling program, based on the velocity-stress
! elastic modeling program

! FD elastic operator module
! This is the main subroutine that calls
! the FD derivatives, sources, etc.
module rtm_mod_dp

  use rtmlib_mod_double
  use sep
  use source_mod_double

  implicit none
  integer, private :: nt, nx, nz
  real, private    :: dt, dx, dz

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine elastic_init(nz_in,dz_in,nx_in,dx_in,dt_in)

    integer             :: nz_in,nx_in
    real                :: dz_in,dx_in,dt_in

    nz = nz_in
    dz = dz_in
    dt = dt_in
    nx = nx_in
    dx = dx_in

  end subroutine

  ! Function that is called by Main
  function elastic_op(vp0, vs0, rho0) result(stat)
    real,dimension(:,:)         :: vp0, vs0, rho0
    integer                     :: stat

    stat=1
    call elastic_op2(vp0, vs0, rho0)
    stat=0
  end function

  ! Subroutine used by function
  subroutine elastic_op2(vp0, vs0, rho0)
    integer                                 :: it, ix, iz
    real, dimension(nz,nx)                  :: vp0,vs0,rho0
    real                                    :: aux1
    integer                                 :: stat

    !Background properties
    real*8, allocatable, dimension(:,:)     :: mi0,l0,l2mi0
    real*8, allocatable, dimension(:,:)     :: rho0_Vx,rho0_Vz

    !model 1
    real*8, allocatable, dimension(:,:,:)   :: Vx_x,Vx_z,Vz_x,Vz_z
    real*8, allocatable, dimension(:,:,:)   :: sigmaxx_x,sigmaxx_z,sigmazz_x,sigmazz_z
    real*8, allocatable, dimension(:,:,:)   :: sigmaxz_x,sigmaxz_z

    !data 1
    real*8, allocatable, dimension(:,:,:)   :: Vx_x_f,Vx_z_f,Vz_x_f,Vz_z_f
    real*8, allocatable, dimension(:,:,:)   :: sigmaxx_x_f,sigmaxx_z_f,sigmazz_x_f,sigmazz_z_f
    real*8, allocatable, dimension(:,:,:)   :: sigmaxz_x_f,sigmaxz_z_f

    !data 2
    real*8, allocatable, dimension(:,:,:)   :: Vx_x2,Vx_z2,Vz_x2,Vz_z2
    real*8, allocatable, dimension(:,:,:)   :: sigmaxx_x2,sigmaxx_z2,sigmazz_x2,sigmazz_z2
    real*8, allocatable, dimension(:,:,:)   :: sigmaxz_x2,sigmaxz_z2

    !model 2
    real*8, allocatable, dimension(:,:,:)   :: Vx_x2_adj,Vx_z2_adj,Vz_x2_adj,Vz_z2_adj
    real*8, allocatable, dimension(:,:,:)   :: sigmaxx_x2_adj,sigmaxx_z2_adj,sigmazz_x2_adj,sigmazz_z2_adj
    real*8, allocatable, dimension(:,:,:)   :: sigmaxz_x2_adj,sigmaxz_z2_adj

    !Test sums
    real*8                                   :: sumvx,sumvz
    real*8                                   :: sumsigmaxx,sumsigmazz
    real*8                                   :: sumsigmaxz

    !Dot product test result
    real*8                                   :: dp1, dp2, error

    !Source terms
    real*8                                   :: source1,source2

    !Output
    real, allocatable, dimension(:,:,:)      :: data1
    real, allocatable, dimension(:)          :: out1,out2

    nt=1000
    write(0,*)nx,nz

    !Allocate lame parameters
    allocate(mi0(nz,nx))
    allocate(l0(nz,nx))
    allocate(l2mi0(nz,nx))
    allocate(rho0_Vx(nz,nx))
    allocate(rho0_Vz(nz,nx))

    !Allocate Velocities and stresses
    !model 1
    allocate(Vx_x(nz,nx,nt))
    allocate(Vx_z(nz,nx,nt))
    allocate(Vz_x(nz,nx,nt))
    allocate(Vz_z(nz,nx,nt))
    allocate(sigmaxx_x(nz,nx,nt))
    allocate(sigmaxx_z(nz,nx,nt))
    allocate(sigmazz_x(nz,nx,nt))
    allocate(sigmazz_z(nz,nx,nt))
    allocate(sigmaxz_x(nz,nx,nt))
    allocate(sigmaxz_z(nz,nx,nt))

    !data 1
    allocate(Vx_x_f(nz,nx,nt))
    allocate(Vx_z_f(nz,nx,nt))
    allocate(Vz_x_f(nz,nx,nt))
    allocate(Vz_z_f(nz,nx,nt))
    allocate(sigmaxx_x_f(nz,nx,nt))
    allocate(sigmaxx_z_f(nz,nx,nt))
    allocate(sigmazz_x_f(nz,nx,nt))
    allocate(sigmazz_z_f(nz,nx,nt))
    allocate(sigmaxz_x_f(nz,nx,nt))
    allocate(sigmaxz_z_f(nz,nx,nt))

    !data 2
    allocate(Vx_x2(nz,nx,nt))
    allocate(Vx_z2(nz,nx,nt))
    allocate(Vz_x2(nz,nx,nt))
    allocate(Vz_z2(nz,nx,nt))
    allocate(sigmaxx_x2(nz,nx,nt))
    allocate(sigmaxx_z2(nz,nx,nt))
    allocate(sigmazz_x2(nz,nx,nt))
    allocate(sigmazz_z2(nz,nx,nt))
    allocate(sigmaxz_x2(nz,nx,nt))
    allocate(sigmaxz_z2(nz,nx,nt))

    !model 2
    allocate(Vx_x2_adj(nz,nx,nt))
    allocate(Vx_z2_adj(nz,nx,nt))
    allocate(Vz_x2_adj(nz,nx,nt))
    allocate(Vz_z2_adj(nz,nx,nt))
    allocate(sigmaxx_x2_adj(nz,nx,nt))
    allocate(sigmaxx_z2_adj(nz,nx,nt))
    allocate(sigmazz_x2_adj(nz,nx,nt))
    allocate(sigmazz_z2_adj(nz,nx,nt))
    allocate(sigmaxz_x2_adj(nz,nx,nt))
    allocate(sigmaxz_z2_adj(nz,nx,nt))

    !Output
    allocate(data1(nz,nx,nt))
    allocate(out1(nt))
    allocate(out2(nt))

    !Initialize values
    !Background properties
    mi0         =0.0
    l0          =0.0
    l2mi0       =0.0
    rho0_Vx     =0.0
    rho0_Vz     =0.0

    !model 1
    Vx_x       =0.0
    Vx_z       =0.0
    Vz_x       =0.0
    Vz_z       =0.0
    sigmaxx_x  =0.0
    sigmaxx_z  =0.0
    sigmazz_x  =0.0
    sigmazz_z  =0.0
    sigmaxz_x  =0.0
    sigmaxz_z  =0.0

    !data 1
    Vx_x_f       =0.0
    Vx_z_f       =0.0
    Vz_x_f       =0.0
    Vz_z_f       =0.0
    sigmaxx_x_f  =0.0
    sigmaxx_z_f  =0.0
    sigmazz_x_f  =0.0
    sigmazz_z_f  =0.0
    sigmaxz_x_f  =0.0
    sigmaxz_z_f  =0.0

    !data 2
    Vx_x2       =0.0
    Vx_z2       =0.0
    Vz_x2       =0.0
    Vz_z2       =0.0
    sigmaxx_x2  =0.0
    sigmaxx_z2  =0.0
    sigmazz_x2  =0.0
    sigmazz_z2  =0.0
    sigmaxz_x2  =0.0
    sigmaxz_z2  =0.0

    !model 2
    Vx_x2_adj       =0.0
    Vx_z2_adj       =0.0
    Vz_x2_adj       =0.0
    Vz_z2_adj       =0.0
    sigmaxx_x2_adj  =0.0
    sigmaxx_z2_adj  =0.0
    sigmazz_x2_adj  =0.0
    sigmazz_z2_adj  =0.0
    sigmaxz_x2_adj  =0.0
    sigmaxz_z2_adj  =0.0

    !Output
    data1 = 0.
    out1 = 0.
    out2 = 0.

    !Define Lame parameters
    !These model properties are still in the original grid, centered at i,j
    do ix=1,nx
      do iz=1,nz
        l0(iz,ix)=vp0(iz,ix)
        mi0(iz,ix)=vs0(iz,ix)
      enddo
    enddo

    !Calculate inverse of density
    !These also calculate the values of density at the velocity staggered grid,
    !using the arithmetic average proposed in Moczo, 2002
    !Rho_Vx

    aux1=(dt/dx)

    !Rho_vx
    do ix=1,nx
      do iz=1,nz-1
        rho0_Vx(iz,ix)=aux1*1.0/(0.5*(rho0(iz+1,ix) + rho0(iz,ix)))
      enddo
    enddo
    !Rho_vz
    do ix=1,nx-1
      do iz=1,nz
        rho0_Vz(iz,ix)=aux1*1.0/(0.5*(rho0(iz,ix+1) + rho0(iz,ix)))
      enddo
    enddo
    !Right side
    do iz=1,nz
      rho0_Vz(iz,nx)=aux1*(1.0/rho0(iz,nx))
    enddo
    !Bottom
    do ix=1,nx
      rho0_Vx(nz,ix)=aux1*(1.0/rho0(nz,ix))
    enddo
    
    !Calculate weigthed average for the lambda variable
    !These put lambda and (lambda+2mi) in their respective staggered grids,
    !using the harmonic average proposed in Moczo, 2002
    l0 = aux1*l0
    mi0=aux1*mi0
    do ix=2,nx-1
      do iz=2,nz-1
        l0(iz,ix)=4.0/((1.0/l0(iz,ix))+(1.0/l0(iz+1,ix))+(1.0/l0(iz,ix+1))+(1.0/l0(iz+1,ix+1)))
        if (mi0(iz,ix)==0.or.mi0(iz+1,ix)==0.or.mi0(iz,ix+1)==0.or.mi0(iz+1,ix+1)==0) then
          mi0(iz,ix)=0.0
        else
          mi0(iz,ix)=4.0/((1.0/mi0(iz,ix))+(1.0/mi0(iz+1,ix))+(1.0/mi0(iz,ix+1))+(1.0/mi0(iz+1,ix+1)))
        endif
      enddo
    enddo
    l2mi0=l0+2*mi0

    !Initialize subroutines
    call rtmmod_init(nz,nx,dx,rho0_Vx,rho0_Vz,l0,mi0,l2mi0)
!    call source_init(20.0,dt,dx,250,nz,nx)

!###################################################
! Create source for model 1 and data 2
!###################################################
!    do it=1,nt
!      call source_op2(it,source1)
!      call explosive_sourcex(source1,100,100,sigmaxx_x(:,:,it))
!      call explosive_sourcez(source1,100,100,sigmazz_z(:,:,it))
!      call source_op3(it,source2)
!      call explosive_sourcex(source2,100,100,sigmaxx_x2_adj(:,:,it))
!      call explosive_sourcez(source2,100,100,sigmazz_z2_adj(:,:,it))
!    enddo

!    out1(:) = sigmazz_z(100,100,:)
!    out2(:) = sigmazz_z2_adj(100,100,:)

!    call to_history("n1", nt, "Dat/source1.H")
!    call to_history("o1", 0., "Dat/source1.H")
!    call to_history("d1", dt, "Dat/source1.H")
         
!    call sep_write(out1,"Dat/source1.H")

!    call to_history("n1", nt, "Dat/source2.H")
!    call to_history("o1", 0., "Dat/source2.H")
!    call to_history("d1", dt, "Dat/source2.H")
         
!    call sep_write(out2,"Dat/source2.H")

!###################################################
! Create random model 1 and data 2
!###################################################

    !model 1
    call random_number(Vx_x)
    call random_number(Vx_z)
    call random_number(Vz_x)
    call random_number(Vz_z)
    call random_number(sigmaxx_x)
    call random_number(sigmaxx_z)
    call random_number(sigmazz_x)
    call random_number(sigmazz_z)
    call random_number(sigmaxz_x)
    call random_number(sigmaxz_z)
    !data 2
    call random_number(Vx_x2)
    call random_number(Vx_z2)
    call random_number(Vz_x2)
    call random_number(Vz_z2)
    call random_number(sigmaxx_x2)
    call random_number(sigmaxx_z2)
    call random_number(sigmazz_x2)
    call random_number(sigmazz_z2)
    call random_number(sigmaxz_x2)
    call random_number(sigmaxz_z2)

    Vx_x = 2.0*Vx_x - 1.0
    Vx_z = 2.0*Vx_z - 1.0
    Vz_x = 2.0*Vz_x - 1.0
    Vz_z = 2.0*Vz_z - 1.0
    sigmaxx_x = 2.0*sigmaxx_x - 1.0
    sigmaxx_z = 2.0*sigmaxx_z - 1.0
    sigmazz_x = 2.0*sigmazz_x - 1.0
    sigmazz_z = 2.0*sigmazz_z - 1.0
    sigmaxz_x = 2.0*sigmaxz_x - 1.0
    sigmaxz_z = 2.0*sigmaxz_z - 1.0
    Vx_x2 = 2.0*Vx_x2 - 1.0
    Vx_z2 = 2.0*Vx_z2 - 1.0
    Vz_x2 = 2.0*Vz_x2 - 1.0
    Vz_z2 = 2.0*Vz_z2 - 1.0
    sigmaxx_x2 = 2.0*sigmaxx_x2 - 1.0
    sigmaxx_z2 = 2.0*sigmaxx_z2 - 1.0
    sigmazz_x2 = 2.0*sigmazz_x2 - 1.0
    sigmazz_z2 = 2.0*sigmazz_z2 - 1.0
    sigmaxz_x2 = 2.0*sigmaxz_x2 - 1.0
    sigmaxz_z2 = 2.0*sigmaxz_z2 - 1.0

    sumvx=0.
    sumvz=0.
    sumsigmaxx=0.
    sumsigmazz=0.
    sumsigmaxz=0.
    
    write(0,*)"============================================================="
    write(0,*)"TEST RANDOM INPUT MODEL 1"
    write(0,*)"============================================================="
    do it=1,nt
      do ix=1,nx
        do iz=1,nz
          sumvx = sumvx + Vx_x(iz,ix,it) + Vx_z(iz,ix,it)
          sumvz = sumvz + Vz_x(iz,ix,it) + Vz_z(iz,ix,it)
          sumsigmaxx = sumsigmaxx + sigmaxx_x(iz,ix,it) + sigmaxx_z(iz,ix,it)
          sumsigmazz = sumsigmazz + sigmazz_x(iz,ix,it) + sigmazz_z(iz,ix,it)
          sumsigmaxz = sumsigmaxz + sigmaxz_x(iz,ix,it) + sigmaxz_z(iz,ix,it)
        enddo
      enddo
    enddo

    write(0,*)sumvx,sumvz,sumsigmaxx,sumsigmazz,sumsigmaxz

    !Time loop for 2nd order time
    write(0,*)"============================================================="
    write(0,*)"BEGIN OF TIME MODELING"
    write(0,*)"============================================================="
!    write(0,*)"Processing:"
!    write(0,*)"0%           25%            50%           75%            100%"
!    write(0,'(a)',advance="no")"["

!###################################################
! Forward
!###################################################
! B d1 = m1
! First time steps
    !it=1
    Vx_x_f(:,:,1) = Vx_x(:,:,1)
    Vx_z_f(:,:,1) = Vx_z(:,:,1)
    Vz_x_f(:,:,1) = Vz_x(:,:,1)
    Vz_z_f(:,:,1) = Vz_z(:,:,1)
    sigmaxx_x_f(:,:,1) = sigmaxx_x(:,:,1)
    sigmaxx_z_f(:,:,1) = sigmaxx_z(:,:,1)
    sigmazz_x_f(:,:,1) = sigmazz_x(:,:,1)
    sigmazz_z_f(:,:,1) = sigmazz_z(:,:,1)
    sigmaxz_x_f(:,:,1) = sigmaxz_x(:,:,1)
    sigmaxz_z_f(:,:,1) = sigmaxz_z(:,:,1)

    !it=2
    !Vx_x
    Vx_x_f(:,:,2) = Vx_x(:,:,2)
    call Ax(.false.,.true.,sigmaxx_x_f(:,:,1),Vx_x_f(:,:,2))
    call Ax(.false.,.true.,sigmaxx_z_f(:,:,1),Vx_x_f(:,:,2))
    !Vx_z
    Vx_z_f(:,:,2) = Vx_z(:,:,2)
    call Bx(.false.,.true.,sigmaxz_x_f(:,:,1),Vx_z_f(:,:,2))
    call Bx(.false.,.true.,sigmaxz_z_f(:,:,1),Vx_z_f(:,:,2))
    !Vz_x
    Vz_x_f(:,:,2) = Vz_x(:,:,2)
    call Az(.false.,.true.,sigmaxz_x_f(:,:,1),Vz_x_f(:,:,2))
    call Az(.false.,.true.,sigmaxz_z_f(:,:,1),Vz_x_f(:,:,2))
    !Vz_z
    Vz_z_f(:,:,2) = Vz_z(:,:,2)
    call Bz(.false.,.true.,sigmazz_x_f(:,:,1),Vz_z_f(:,:,2))
    call Bz(.false.,.true.,sigmazz_z_f(:,:,1),Vz_z_f(:,:,2))
    !sigmaxx_x
    sigmaxx_x_f(:,:,2) = sigmaxx_x(:,:,2)
    call D(.false.,.true.,Vx_x_f(:,:,1),sigmaxx_x_f(:,:,2))
    call D(.false.,.true.,Vx_z_f(:,:,1),sigmaxx_x_f(:,:,2))
    !sigmaxx_z
    sigmaxx_z_f(:,:,2) = sigmaxx_z(:,:,2)
    call E(.false.,.true.,Vz_x_f(:,:,1),sigmaxx_z_f(:,:,2))
    call E(.false.,.true.,Vz_z_f(:,:,1),sigmaxx_z_f(:,:,2))
    !sigmazz_x
    sigmazz_x_f(:,:,2) = sigmazz_x(:,:,2)
    call F(.false.,.true.,Vx_x_f(:,:,1),sigmazz_x_f(:,:,2))
    call F(.false.,.true.,Vx_z_f(:,:,1),sigmazz_x_f(:,:,2))
    !sigmazz_z
    sigmazz_z_f(:,:,2) = sigmazz_z(:,:,2)
    call G(.false.,.true.,Vz_x_f(:,:,1),sigmazz_z_f(:,:,2))
    call G(.false.,.true.,Vz_z_f(:,:,1),sigmazz_z_f(:,:,2))
    !sigmaxz_x
    sigmaxz_x_f(:,:,2) = sigmaxz_x(:,:,2)
    call H(.false.,.true.,Vx_x_f(:,:,1),sigmaxz_x_f(:,:,2))
    call H(.false.,.true.,Vx_z_f(:,:,1),sigmaxz_x_f(:,:,2))
    !sigmaxz_z
    sigmaxz_z_f(:,:,2) = sigmaxz_z(:,:,2)
    call J(.false.,.true.,Vz_x_f(:,:,1),sigmaxz_z_f(:,:,2))
    call J(.false.,.true.,Vz_z_f(:,:,1),sigmaxz_z_f(:,:,2))

    !it=3 until it=nt
    do it=3,nt
      !Vx_x
      Vx_x_f(:,:,it) = Vx_x(:,:,it) + Vx_x_f(:,:,(it-2))
      call Ax(.false.,.true.,sigmaxx_x_f(:,:,(it-1)),Vx_x_f(:,:,it))
      call Ax(.false.,.true.,sigmaxx_z_f(:,:,(it-1)),Vx_x_f(:,:,it))
      !Vx_z
      Vx_z_f(:,:,it) = Vx_z(:,:,it) + Vx_z_f(:,:,(it-2))
      call Bx(.false.,.true.,sigmaxz_x_f(:,:,(it-1)),Vx_z_f(:,:,it))
      call Bx(.false.,.true.,sigmaxz_z_f(:,:,(it-1)),Vx_z_f(:,:,it))
      !Vz_x
      Vz_x_f(:,:,it) = Vz_x(:,:,it) + Vz_x_f(:,:,(it-2))
      call Az(.false.,.true.,sigmaxz_x_f(:,:,(it-1)),Vz_x_f(:,:,it))
      call Az(.false.,.true.,sigmaxz_z_f(:,:,(it-1)),Vz_x_f(:,:,it))
      !Vz_z
      Vz_z_f(:,:,it) = Vz_z(:,:,it) + Vz_z_f(:,:,(it-2))
      call Bz(.false.,.true.,sigmazz_x_f(:,:,(it-1)),Vz_z_f(:,:,it))
      call Bz(.false.,.true.,sigmazz_z_f(:,:,(it-1)),Vz_z_f(:,:,it))
      !sigmaxx_x
      sigmaxx_x_f(:,:,it) = sigmaxx_x(:,:,it) + sigmaxx_x_f(:,:,(it-2))
      call D(.false.,.true.,Vx_x_f(:,:,(it-1)),sigmaxx_x_f(:,:,it))
      call D(.false.,.true.,Vx_z_f(:,:,(it-1)),sigmaxx_x_f(:,:,it))
      !sigmaxx_z
      sigmaxx_z_f(:,:,it) = sigmaxx_z(:,:,it) + sigmaxx_z_f(:,:,(it-2))
      call E(.false.,.true.,Vz_x_f(:,:,(it-1)),sigmaxx_z_f(:,:,it))
      call E(.false.,.true.,Vz_z_f(:,:,(it-1)),sigmaxx_z_f(:,:,it))
      !sigmazz_x
      sigmazz_x_f(:,:,it) = sigmazz_x(:,:,it) + sigmazz_x_f(:,:,(it-2))
      call F(.false.,.true.,Vx_x_f(:,:,(it-1)),sigmazz_x_f(:,:,it))
      call F(.false.,.true.,Vx_z_f(:,:,(it-1)),sigmazz_x_f(:,:,it))
      !sigmazz_z
      sigmazz_z_f(:,:,it) = sigmazz_z(:,:,it) + sigmazz_z_f(:,:,(it-2))
      call G(.false.,.true.,Vz_x_f(:,:,(it-1)),sigmazz_z_f(:,:,it))
      call G(.false.,.true.,Vz_z_f(:,:,(it-1)),sigmazz_z_f(:,:,it))
      !sigmaxz_x
      sigmaxz_x_f(:,:,it) = sigmaxz_x(:,:,it) + sigmaxz_x_f(:,:,(it-2))
      call H(.false.,.true.,Vx_x_f(:,:,(it-1)),sigmaxz_x_f(:,:,it))
      call H(.false.,.true.,Vx_z_f(:,:,(it-1)),sigmaxz_x_f(:,:,it))
      !sigmaxz_z
      sigmaxz_z_f(:,:,it) = sigmaxz_z(:,:,it) + sigmaxz_z_f(:,:,(it-2))
      call J(.false.,.true.,Vz_x_f(:,:,(it-1)),sigmaxz_z_f(:,:,it))
      call J(.false.,.true.,Vz_z_f(:,:,(it-1)),sigmaxz_z_f(:,:,it))
    enddo

    sumvx=0.
    sumvz=0.
    sumsigmaxx=0.
    sumsigmazz=0.
    sumsigmaxz=0.

    write(0,*)"============================================================="
    write(0,*)"TEST RANDOM OUTPUT DATA 1"
    write(0,*)"============================================================="
    do it=1,nt
      do ix=1,nx
        do iz=1,nz
          sumvx = sumvx + Vx_x_f(iz,ix,it) + Vx_z_f(iz,ix,it)
          sumvz = sumvz + Vz_x_f(iz,ix,it) + Vz_z_f(iz,ix,it)
          sumsigmaxx = sumsigmaxx + sigmaxx_x_f(iz,ix,it) + sigmaxx_z_f(iz,ix,it)
          sumsigmazz = sumsigmazz + sigmazz_x_f(iz,ix,it) + sigmazz_z_f(iz,ix,it)
          sumsigmaxz = sumsigmaxz + sigmaxz_x_f(iz,ix,it) + sigmaxz_z_f(iz,ix,it)
        enddo
      enddo
    enddo

    write(0,*)sumvx,sumvz,sumsigmaxx,sumsigmazz,sumsigmaxz

!###################################################
! Forward to create data 2
!###################################################
! B d2 = m2
! First time steps
    !it=1
!    Vx_x2(:,:,1) = Vx_x2_adj(:,:,1)
!    Vx_z2(:,:,1) = Vx_z2_adj(:,:,1)
!    Vz_x2(:,:,1) = Vz_x2_adj(:,:,1)
!    Vz_z2(:,:,1) = Vz_z2_adj(:,:,1)
!    sigmaxx_x2(:,:,1) = sigmaxx_x2_adj(:,:,1)
!    sigmaxx_z2(:,:,1) = sigmaxx_z2_adj(:,:,1)
!    sigmazz_x2(:,:,1) = sigmazz_x2_adj(:,:,1)
!    sigmazz_z2(:,:,1) = sigmazz_z2_adj(:,:,1)
!    sigmaxz_x2(:,:,1) = sigmaxz_x2_adj(:,:,1)
!    sigmaxz_z2(:,:,1) = sigmaxz_z2_adj(:,:,1)

    !it=2
    !Vx_x
!    Vx_x2(:,:,2) = Vx_x2_adj(:,:,2)
!    call Ax(.false.,.true.,sigmaxx_x2(:,:,1),Vx_x2(:,:,2))
!    call Ax(.false.,.true.,sigmaxx_z2(:,:,1),Vx_x2(:,:,2))
    !Vx_z
!    Vx_z2(:,:,2) = Vx_z2_adj(:,:,2)
!    call Bx(.false.,.true.,sigmaxz_x2(:,:,1),Vx_z2(:,:,2))
!    call Bx(.false.,.true.,sigmaxz_z2(:,:,1),Vx_z2(:,:,2))
    !Vz_x
!    Vz_x2(:,:,2) = Vz_x2_adj(:,:,2)
!    call Az(.false.,.true.,sigmaxz_x2(:,:,1),Vz_x2(:,:,2))
!    call Az(.false.,.true.,sigmaxz_z2(:,:,1),Vz_x2(:,:,2))
    !Vz_z
!    Vz_z2(:,:,2) = Vz_z2_adj(:,:,2)
!    call Bz(.false.,.true.,sigmazz_x2(:,:,1),Vz_z2(:,:,2))
!    call Bz(.false.,.true.,sigmazz_z2(:,:,1),Vz_z2(:,:,2))
    !sigmaxx_x
!    sigmaxx_x2(:,:,2) = sigmaxx_x2_adj(:,:,2)
!    call D(.false.,.true.,Vx_x2(:,:,1),sigmaxx_x2(:,:,2))
!    call D(.false.,.true.,Vx_z2(:,:,1),sigmaxx_x2(:,:,2))
    !sigmaxx_z
!    sigmaxx_z2(:,:,2) = sigmaxx_z2_adj(:,:,2)
!    call E(.false.,.true.,Vz_x2(:,:,1),sigmaxx_z2(:,:,2))
!    call E(.false.,.true.,Vz_z2(:,:,1),sigmaxx_z2(:,:,2))
    !sigmazz_x
!    sigmazz_x2(:,:,2) = sigmazz_x2_adj(:,:,2)
!    call F(.false.,.true.,Vx_x2(:,:,1),sigmazz_x2(:,:,2))
!    call F(.false.,.true.,Vx_z2(:,:,1),sigmazz_x2(:,:,2))
    !sigmazz_z
!    sigmazz_z2(:,:,2) = sigmazz_z2_adj(:,:,2)
!    call G(.false.,.true.,Vz_x2(:,:,1),sigmazz_z2(:,:,2))
!    call G(.false.,.true.,Vz_z2(:,:,1),sigmazz_z2(:,:,2))
    !sigmaxz_x
!    sigmaxz_x2(:,:,2) = sigmaxz_x2_adj(:,:,2)
!    call H(.false.,.true.,Vx_x2(:,:,1),sigmaxz_x2(:,:,2))
!    call H(.false.,.true.,Vx_z2(:,:,1),sigmaxz_x2(:,:,2))
    !sigmaxz_z
!    sigmaxz_z2(:,:,2) = sigmaxz_z2_adj(:,:,2)
!    call J(.false.,.true.,Vz_x2(:,:,1),sigmaxz_z2(:,:,2))
!    call J(.false.,.true.,Vz_z2(:,:,1),sigmaxz_z2(:,:,2))

    !it=3 until it=nt
!    do it=3,nt
      !Vx_x
!      Vx_x2(:,:,it) = Vx_x2_adj(:,:,it) + Vx_x2(:,:,(it-2))
!      call Ax(.false.,.true.,sigmaxx_x2(:,:,(it-1)),Vx_x2(:,:,it))
!      call Ax(.false.,.true.,sigmaxx_z2(:,:,(it-1)),Vx_x2(:,:,it))
      !Vx_z
!      Vx_z2(:,:,it) = Vx_z2_adj(:,:,it) + Vx_z2(:,:,(it-2))
!      call Bx(.false.,.true.,sigmaxz_x2(:,:,(it-1)),Vx_z2(:,:,it))
!      call Bx(.false.,.true.,sigmaxz_z2(:,:,(it-1)),Vx_z2(:,:,it))
      !Vz_x
!      Vz_x2(:,:,it) = Vz_x2_adj(:,:,it) + Vz_x2(:,:,(it-2))
!      call Az(.false.,.true.,sigmaxz_x2(:,:,(it-1)),Vz_x2(:,:,it))
!      call Az(.false.,.true.,sigmaxz_z2(:,:,(it-1)),Vz_x2(:,:,it))
      !Vz_z
!      Vz_z2(:,:,it) = Vz_z2_adj(:,:,it) + Vz_z2(:,:,(it-2))
!      call Bz(.false.,.true.,sigmazz_x2(:,:,(it-1)),Vz_z2(:,:,it))
!      call Bz(.false.,.true.,sigmazz_z2(:,:,(it-1)),Vz_z2(:,:,it))
      !sigmaxx_x
!      sigmaxx_x2(:,:,it) = sigmaxx_x2_adj(:,:,it) + sigmaxx_x2(:,:,(it-2))
!      call D(.false.,.true.,Vx_x2(:,:,(it-1)),sigmaxx_x2(:,:,it))
!      call D(.false.,.true.,Vx_z2(:,:,(it-1)),sigmaxx_x2(:,:,it))
      !sigmaxx_z
!      sigmaxx_z2(:,:,it) = sigmaxx_z2_adj(:,:,it) + sigmaxx_z2(:,:,(it-2))
!      call E(.false.,.true.,Vz_x2(:,:,(it-1)),sigmaxx_z2(:,:,it))
!      call E(.false.,.true.,Vz_z2(:,:,(it-1)),sigmaxx_z2(:,:,it))
      !sigmazz_x
!      sigmazz_x2(:,:,it) = sigmazz_x2_adj(:,:,it) + sigmazz_x2(:,:,(it-2))
!      call F(.false.,.true.,Vx_x2(:,:,(it-1)),sigmazz_x2(:,:,it))
!      call F(.false.,.true.,Vx_z2(:,:,(it-1)),sigmazz_x2(:,:,it))
      !sigmazz_z
!      sigmazz_z2(:,:,it) = sigmazz_z2_adj(:,:,it) + sigmazz_z2(:,:,(it-2))
!      call G(.false.,.true.,Vz_x2(:,:,(it-1)),sigmazz_z2(:,:,it))
!      call G(.false.,.true.,Vz_z2(:,:,(it-1)),sigmazz_z2(:,:,it))
      !sigmaxz_x
!      sigmaxz_x2(:,:,it) = sigmaxz_x2_adj(:,:,it) + sigmaxz_x2(:,:,(it-2))
!      call H(.false.,.true.,Vx_x2(:,:,(it-1)),sigmaxz_x2(:,:,it))
!      call H(.false.,.true.,Vx_z2(:,:,(it-1)),sigmaxz_x2(:,:,it))
      !sigmaxz_z
!      sigmaxz_z2(:,:,it) = sigmaxz_z2_adj(:,:,it) + sigmaxz_z2(:,:,(it-2))
!      call J(.false.,.true.,Vz_x2(:,:,(it-1)),sigmaxz_z2(:,:,it))
!      call J(.false.,.true.,Vz_z2(:,:,(it-1)),sigmaxz_z2(:,:,it))
!    enddo

!    data1 = sigmaxx_x2 + sigmaxx_z2 + sigmazz_x2 + sigmazz_z2

!    call to_history("n1", nz, "Dat/data2.H")
!    call to_history("o1", 0., "Dat/data2.H")
!    call to_history("d1", dz, "Dat/data2.H")
!    call to_history("n2", nx, "Dat/data2.H")
!    call to_history("o2", 0., "Dat/data2.H")
!    call to_history("d2", dx, "Dat/data2.H")
!    call to_history("n3", nt, "Dat/data2.H")
!    call to_history("o3", 0., "Dat/data2.H")
!    call to_history("d3", dt, "Dat/data2.H")
         
!    call sep_write(data1,"Dat/data2.H")

    sumvx=0.
    sumvz=0.
    sumsigmaxx=0.
    sumsigmazz=0.
    sumsigmaxz=0.

    write(0,*)"============================================================="
    write(0,*)"TEST RANDOM INPUT DATA 2"
    write(0,*)"============================================================="
    do it=1,nt
      do ix=1,nx
        do iz=1,nz
          sumvx = sumvx + Vx_x2(iz,ix,it) + Vx_z2(iz,ix,it)
          sumvz = sumvz + Vz_x2(iz,ix,it) + Vz_z2(iz,ix,it)
          sumsigmaxx = sumsigmaxx + sigmaxx_x2(iz,ix,it) + sigmaxx_z2(iz,ix,it)
          sumsigmazz = sumsigmazz + sigmazz_x2(iz,ix,it) + sigmazz_z2(iz,ix,it)
          sumsigmaxz = sumsigmaxz + sigmaxz_x2(iz,ix,it) + sigmaxz_z2(iz,ix,it)
        enddo
      enddo
    enddo

    write(0,*)sumvx,sumvz,sumsigmaxx,sumsigmazz,sumsigmaxz

    Vx_x2_adj       =0.0
    Vx_z2_adj       =0.0
    Vz_x2_adj       =0.0
    Vz_z2_adj       =0.0
    sigmaxx_x2_adj  =0.0
    sigmaxx_z2_adj  =0.0
    sigmazz_x2_adj  =0.0
    sigmazz_z2_adj  =0.0
    sigmaxz_x2_adj  =0.0
    sigmaxz_z2_adj  =0.0

!###################################################
! Adjoint application
!###################################################
! B*m2 = d2
! Last time steps
    !it=nt
    Vx_x2_adj(:,:,nt) = Vx_x2(:,:,nt)
    Vx_z2_adj(:,:,nt) = Vx_z2(:,:,nt)
    Vz_x2_adj(:,:,nt) = Vz_x2(:,:,nt)
    Vz_z2_adj(:,:,nt) = Vz_z2(:,:,nt)
    sigmaxx_x2_adj(:,:,nt) = sigmaxx_x2(:,:,nt)
    sigmaxx_z2_adj(:,:,nt) = sigmaxx_z2(:,:,nt)
    sigmazz_x2_adj(:,:,nt) = sigmazz_x2(:,:,nt)
    sigmazz_z2_adj(:,:,nt) = sigmazz_z2(:,:,nt)
    sigmaxz_x2_adj(:,:,nt) = sigmaxz_x2(:,:,nt)
    sigmaxz_z2_adj(:,:,nt) = sigmaxz_z2(:,:,nt)

    !it=nt-1
    !Vx_x
    Vx_x2_adj(:,:,nt-1) = Vx_x2(:,:,nt-1)
    call D(.true.,.true.,Vx_x2_adj(:,:,(nt-1)),sigmaxx_x2_adj(:,:,nt))    
    call F(.true.,.true.,Vx_x2_adj(:,:,(nt-1)),sigmazz_x2_adj(:,:,nt))
    call H(.true.,.true.,Vx_x2_adj(:,:,(nt-1)),sigmaxz_x2_adj(:,:,nt))
    !Vx_z
    Vx_z2_adj(:,:,nt-1) = Vx_z2(:,:,nt-1)
    call D(.true.,.true.,Vx_z2_adj(:,:,(nt-1)),sigmaxx_x2_adj(:,:,nt))    
    call F(.true.,.true.,Vx_z2_adj(:,:,(nt-1)),sigmazz_x2_adj(:,:,nt))
    call H(.true.,.true.,Vx_z2_adj(:,:,(nt-1)),sigmaxz_x2_adj(:,:,nt))
    !Vz_x
    Vz_x2_adj(:,:,nt-1) = Vz_x2(:,:,nt-1)
    call E(.true.,.true.,Vz_x2_adj(:,:,(nt-1)),sigmaxx_z2_adj(:,:,nt))    
    call G(.true.,.true.,Vz_x2_adj(:,:,(nt-1)),sigmazz_z2_adj(:,:,nt))
    call J(.true.,.true.,Vz_x2_adj(:,:,(nt-1)),sigmaxz_z2_adj(:,:,nt))
    !Vz_z
    Vz_z2_adj(:,:,nt-1) = Vz_z2(:,:,nt-1)
    call E(.true.,.true.,Vz_z2_adj(:,:,(nt-1)),sigmaxx_z2_adj(:,:,nt))    
    call G(.true.,.true.,Vz_z2_adj(:,:,(nt-1)),sigmazz_z2_adj(:,:,nt))
    call J(.true.,.true.,Vz_z2_adj(:,:,(nt-1)),sigmaxz_z2_adj(:,:,nt))
    !sigmaxx_x
    sigmaxx_x2_adj(:,:,nt-1) = sigmaxx_x2(:,:,nt-1)
    call Ax(.true.,.true.,sigmaxx_x2_adj(:,:,(nt-1)),Vx_x2_adj(:,:,nt))
    !sigmaxx_z
    sigmaxx_z2_adj(:,:,nt-1) = sigmaxx_z2(:,:,nt-1)
    call Ax(.true.,.true.,sigmaxx_z2_adj(:,:,(nt-1)),Vx_x2_adj(:,:,nt))
    !sigmazz_x
    sigmazz_x2_adj(:,:,nt-1) = sigmazz_x2(:,:,nt-1)
    call Bz(.true.,.true.,sigmazz_x2_adj(:,:,(nt-1)),Vz_z2_adj(:,:,nt))
    !sigmazz_z
    sigmazz_z2_adj(:,:,nt-1) = sigmazz_z2(:,:,nt-1)
    call Bz(.true.,.true.,sigmazz_z2_adj(:,:,(nt-1)),Vz_z2_adj(:,:,nt))
    !sigmaxz_x
    sigmaxz_x2_adj(:,:,nt-1) = sigmaxz_x2(:,:,nt-1)
    call Bx(.true.,.true.,sigmaxz_x2_adj(:,:,(nt-1)),Vx_z2_adj(:,:,nt))
    call Az(.true.,.true.,sigmaxz_x2_adj(:,:,(nt-1)),Vz_x2_adj(:,:,nt))
    !sigmaxz_z
    sigmaxz_z2_adj(:,:,nt-1) = sigmaxz_z2(:,:,nt-1)
    call Bx(.true.,.true.,sigmaxz_z2_adj(:,:,(nt-1)),Vx_z2_adj(:,:,nt))
    call Az(.true.,.true.,sigmaxz_z2_adj(:,:,(nt-1)),Vz_x2_adj(:,:,nt))

    !it=nt-2 until it=1
    do it=nt-2,1,-1
      !Vx_x
      Vx_x2_adj(:,:,it) = Vx_x2(:,:,it) + Vx_x2_adj(:,:,(it+2))
      call D(.true.,.true.,Vx_x2_adj(:,:,it),sigmaxx_x2_adj(:,:,(it+1)))    
      call F(.true.,.true.,Vx_x2_adj(:,:,it),sigmazz_x2_adj(:,:,(it+1)))
      call H(.true.,.true.,Vx_x2_adj(:,:,it),sigmaxz_x2_adj(:,:,(it+1)))
      !Vx_z
      Vx_z2_adj(:,:,it) = Vx_z2(:,:,it) + Vx_z2_adj(:,:,(it+2))
      call D(.true.,.true.,Vx_z2_adj(:,:,it),sigmaxx_x2_adj(:,:,(it+1)))    
      call F(.true.,.true.,Vx_z2_adj(:,:,it),sigmazz_x2_adj(:,:,(it+1)))
      call H(.true.,.true.,Vx_z2_adj(:,:,it),sigmaxz_x2_adj(:,:,(it+1)))
      !Vz_x
      Vz_x2_adj(:,:,it) = Vz_x2(:,:,it) + Vz_x2_adj(:,:,(it+2))
      call E(.true.,.true.,Vz_x2_adj(:,:,it),sigmaxx_z2_adj(:,:,(it+1)))    
      call G(.true.,.true.,Vz_x2_adj(:,:,it),sigmazz_z2_adj(:,:,(it+1)))
      call J(.true.,.true.,Vz_x2_adj(:,:,it),sigmaxz_z2_adj(:,:,(it+1)))
      !Vz_z
      Vz_z2_adj(:,:,it) = Vz_z2(:,:,it) + Vz_z2_adj(:,:,(it+2))
      call E(.true.,.true.,Vz_z2_adj(:,:,it),sigmaxx_z2_adj(:,:,(it+1)))    
      call G(.true.,.true.,Vz_z2_adj(:,:,it),sigmazz_z2_adj(:,:,(it+1)))
      call J(.true.,.true.,Vz_z2_adj(:,:,it),sigmaxz_z2_adj(:,:,(it+1)))
      !sigmaxx_x
      sigmaxx_x2_adj(:,:,it) = sigmaxx_x2(:,:,it) + sigmaxx_x2_adj(:,:,(it+2))
      call Ax(.true.,.true.,sigmaxx_x2_adj(:,:,it),Vx_x2_adj(:,:,(it+1)))
      !sigmaxx_z
      sigmaxx_z2_adj(:,:,it) = sigmaxx_z2(:,:,it) + sigmaxx_z2_adj(:,:,(it+2))
      call Ax(.true.,.true.,sigmaxx_z2_adj(:,:,it),Vx_x2_adj(:,:,(it+1)))
      !sigmazz_x
      sigmazz_x2_adj(:,:,it) = sigmazz_x2(:,:,it) + sigmazz_x2_adj(:,:,(it+2))
      call Bz(.true.,.true.,sigmazz_x2_adj(:,:,it),Vz_z2_adj(:,:,(it+1)))
      !sigmazz_z
      sigmazz_z2_adj(:,:,it) = sigmazz_z2(:,:,it) + sigmazz_z2_adj(:,:,(it+2))
      call Bz(.true.,.true.,sigmazz_z2_adj(:,:,it),Vz_z2_adj(:,:,(it+1)))
      !sigmaxz_x
      sigmaxz_x2_adj(:,:,it) = sigmaxz_x2(:,:,it) + sigmaxz_x2_adj(:,:,(it+2))
      call Bx(.true.,.true.,sigmaxz_x2_adj(:,:,it),Vx_z2_adj(:,:,(it+1)))
      call Az(.true.,.true.,sigmaxz_x2_adj(:,:,it),Vz_x2_adj(:,:,(it+1)))
      !sigmaxz_z
      sigmaxz_z2_adj(:,:,it) = sigmaxz_z2(:,:,it) + sigmaxz_z2_adj(:,:,(it+2))
      call Bx(.true.,.true.,sigmaxz_z2_adj(:,:,it),Vx_z2_adj(:,:,(it+1)))
      call Az(.true.,.true.,sigmaxz_z2_adj(:,:,it),Vz_x2_adj(:,:,(it+1)))
    enddo
    
    sumvx=0.
    sumvz=0.
    sumsigmaxx=0.
    sumsigmazz=0.
    sumsigmaxz=0.

    write(0,*)"============================================================="
    write(0,*)"TEST RANDOM OUTPUT MODEL 2"
    write(0,*)"============================================================="
    do it=1,nt
      do ix=1,nx
        do iz=1,nz
          sumvx = sumvx + Vx_x2_adj(iz,ix,it) + Vx_z2_adj(iz,ix,it)
          sumvz = sumvz + Vz_x2_adj(iz,ix,it) + Vz_z2_adj(iz,ix,it)
          sumsigmaxx = sumsigmaxx + sigmaxx_x2_adj(iz,ix,it) + sigmaxx_z2_adj(iz,ix,it)
          sumsigmazz = sumsigmazz + sigmazz_x2_adj(iz,ix,it) + sigmazz_z2_adj(iz,ix,it)
          sumsigmaxz = sumsigmaxz + sigmaxz_x2_adj(iz,ix,it) + sigmaxz_z2_adj(iz,ix,it)
        enddo
      enddo
    enddo

    write(0,*)sumvx,sumvz,sumsigmaxx,sumsigmazz,sumsigmaxz

    write(0,*)"============================================================="
    write(0,*)"DOT PRODUCT TEST"
    write(0,*)"============================================================="
    write(0,*)"Dot product of data 1 and data 2"
    do it=1,nt
      do ix=1,nx
        do iz=1,nz
          dp1 = dp1 + Vx_x2(iz,ix,it)*Vx_x_f(iz,ix,it)+&
                      Vx_z2(iz,ix,it)*Vx_z_f(iz,ix,it)+&
                      Vz_x2(iz,ix,it)*Vz_x_f(iz,ix,it)+&
                      Vz_z2(iz,ix,it)*Vz_z_f(iz,ix,it)+&
                      sigmaxx_x2(iz,ix,it)*sigmaxx_x_f(iz,ix,it)+&
                      sigmaxx_z2(iz,ix,it)*sigmaxx_z_f(iz,ix,it)+&
                      sigmazz_x2(iz,ix,it)*sigmazz_x_f(iz,ix,it)+&
                      sigmazz_z2(iz,ix,it)*sigmazz_z_f(iz,ix,it)+&
                      sigmaxz_x2(iz,ix,it)*sigmaxz_x_f(iz,ix,it)+&
                      sigmaxz_z2(iz,ix,it)*sigmaxz_z_f(iz,ix,it)
        enddo
      enddo
    enddo
    write(0,*)"dp1=",dp1

    write(0,*)"Dot product of model 1 and model 2"
    do it=1,nt
      do ix=1,nx
        do iz=1,nz
          dp2 = dp2 + Vx_x(iz,ix,it)*Vx_x2_adj(iz,ix,it)+&
                      Vx_z(iz,ix,it)*Vx_z2_adj(iz,ix,it)+&
                      Vz_x(iz,ix,it)*Vz_x2_adj(iz,ix,it)+&
                      Vz_z(iz,ix,it)*Vz_z2_adj(iz,ix,it)+&
                      sigmaxx_x(iz,ix,it)*sigmaxx_x2_adj(iz,ix,it)+&
                      sigmaxx_z(iz,ix,it)*sigmaxx_z2_adj(iz,ix,it)+&
                      sigmazz_x(iz,ix,it)*sigmazz_x2_adj(iz,ix,it)+&
                      sigmazz_z(iz,ix,it)*sigmazz_z2_adj(iz,ix,it)+&
                      sigmaxz_x(iz,ix,it)*sigmaxz_x2_adj(iz,ix,it)+&
                      sigmaxz_z(iz,ix,it)*sigmaxz_z2_adj(iz,ix,it)
        enddo
      enddo
    enddo
    write(0,*)"dp2=",dp2

    error = (dp1-dp2)/dp2
    write(0,*)"error=",error

!    data1 = 0.
!    data1 = sigmaxx_x_f + sigmaxx_z_f + sigmazz_x_f + sigmazz_z_f

!    call to_history("n1", nz, "Dat/data1.H")
!    call to_history("o1", 0., "Dat/data1.H")
!    call to_history("d1", dz, "Dat/data1.H")
!    call to_history("n2", nx, "Dat/data1.H")
!    call to_history("o2", 0., "Dat/data1.H")
!    call to_history("d2", dx, "Dat/data1.H")
!    call to_history("n3", nt, "Dat/data1.H")
!    call to_history("o3", 0., "Dat/data1.H")
!    call to_history("d3", dt, "Dat/data1.H")
         
!    call sep_write(data1,"Dat/data1.H")

!    data1 = 0.
!    data1 = sigmaxx_x2_adj + sigmaxx_z2_adj + sigmazz_x2_adj + sigmazz_z2_adj

!    call to_history("n1", nz, "Dat/data3.H")
!    call to_history("o1", 0., "Dat/data3.H")
!    call to_history("d1", dz, "Dat/data3.H")
!    call to_history("n2", nx, "Dat/data3.H")
!    call to_history("o2", 0., "Dat/data3.H")
!    call to_history("d2", dx, "Dat/data3.H")
!    call to_history("n3", nt, "Dat/data3.H")
!    call to_history("o3", 0., "Dat/data3.H")
!    call to_history("d3", dt, "Dat/data3.H")
         
!    call sep_write(data1,"Dat/data3.H")

  end subroutine

end module
