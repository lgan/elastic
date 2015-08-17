! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 12th, 2015

! RTM code that is supposed to take a forward/adjoint
! input, so that it can be easily tested/put into an
! iterative solver.

module rtm_mod_adj

  use model_mod
  use rtmlib_mod
  use sep
  use source_mod
  use abc_mod
  use image_mod_rtm

  implicit none
  integer, private :: nt, nx, nz, xsource, zsource, abc
  real, private    :: dt, dx, dz, fp, aux1
  logical, private :: free_surface, lame
  integer, private :: zrec, xrecbeg, xrecend, xreclen

  real, pointer, dimension(:,:),private :: m1,m2,m3
  real, allocatable, dimension(:,:), private :: m1_x,m1_z
  real, allocatable, dimension(:), private   :: source

  !Density variation for x and z derivatives (Born case only)
  real, allocatable, dimension(:,:),private  :: dm1_x, dm1_z
  !Data
  real, pointer, dimension(:,:),private :: Vx,Vz,sigmaxx,sigmazz,sigmaxz
  !Model
  real, pointer, dimension(:,:),private :: dm1,dm2,dm3
  !Temporary wavefields
  real, allocatable, dimension(:,:,:),private   :: Vx0,Vz0
  real, allocatable, dimension(:,:,:),private   :: sigmaxx0,sigmazz0
  real, allocatable, dimension(:,:,:),private   :: sigmaxz0

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine rtm_init(nz_in,dz_in,nx_in,dx_in,dt_in,xsource_in,zsource_in,zrec_in,fp_in,nt_in,abc_in,free_surface_in,m1_in,m2_in,m3_in,lame_in)

    integer             :: nz_in,nx_in, xsource_in, zsource_in, nt_in, abc_in,zrec_in
    real                :: dz_in,dx_in,dt_in, fp_in
    logical             :: free_surface_in, lame_in
    real, dimension(:,:),target  :: m1_in,m2_in,m3_in

    nz = nz_in
    dz = dz_in
    dt = dt_in
    nx = nx_in
    dx = dx_in
    xsource = xsource_in
    zsource = zsource_in
    fp = fp_in
    nt = nt_in
    abc = abc_in
    free_surface = free_surface_in
    zrec = zrec_in

    lame = lame_in

    !Input background model
    m1 => m1_in
    m2 => m2_in
    m3 => m3_in

    !Allocate lame parameters
    allocate(m1_x(nz,nx))
    allocate(m1_z(nz,nx))

    !Allocate source
    allocate(source(nt))

    !Allocate Velocities and stresses
    !Temporary wavefields
    allocate(Vx0(nz,nx,nt))
    allocate(Vz0(nz,nx,nt))
    allocate(sigmaxx0(nz,nx,nt))
    allocate(sigmazz0(nz,nx,nt))
    allocate(sigmaxz0(nz,nx,nt))

    !Initialize values
    !Background properties
    m1_x     =0.0
    m1_z     =0.0

    !Source
    source   =0.0

    !Initialize subroutines
    call model_op(nz,nx,dx,dt,m1,m2,m3,m1_x,m1_z)
    call source_init(fp,dt,dx,250,nz,nx,nt)
    call source_op2(source)

    !Initialize absorbing boundary
    call abc_init(nz,nx,dx,dt,abc)

    !Initialize derivatives
    call rtmmod_init(nz,nx,dx,dt,m1,m1_x,m1_z,m2,m3,lame)

  end subroutine

  ! Function that is called by Main
  function rtm_op(adj,dm1_in,dm2_in,dm3_in,Vx_in,Vz_in,sigmaxx_in,sigmazz_in,sigmaxz_in) result(stat)
    logical,intent(in)          :: adj!,add !adj=true is RTM, false is Born
    real,dimension(:,:),target    :: dm1_in,dm2_in,dm3_in
    real,dimension(:,:),target    :: Vx_in,Vz_in,sigmaxx_in,sigmazz_in,sigmaxz_in
    integer                     :: stat

    stat=1
    call rtm_op2(adj,dm1_in,dm2_in,dm3_in,Vx_in,Vz_in,sigmaxx_in,sigmazz_in,sigmaxz_in)
    stat=0
  end function

  ! Subroutine used by function
  subroutine rtm_op2(adj,dm1_in,dm2_in,dm3_in,Vx_in,Vz_in,sigmaxx_in,sigmazz_in,sigmaxz_in)
    logical,intent(in)                      :: adj!,add
    real,dimension(:,:),target            :: dm1_in,dm2_in,dm3_in
    real,dimension(:,:),target            :: Vx_in,Vz_in,sigmaxx_in,sigmazz_in,sigmaxz_in

    !General variables
    integer                                 :: it, ix, iz

    !Model
    dm1 => dm1_in
    dm2 => dm2_in
    dm3 => dm3_in
    write(0,*) "Model pointers ==> OK"

    !Data
    Vx => Vx_in
    Vz => Vz_in
    sigmaxx => sigmaxx_in
    sigmazz => sigmazz_in
    sigmaxz => sigmaxz_in
    write(0,*) "Wavefield pointers ==> OK"

    if (.not.adj) then
      allocate(dm1_x(nz,nx))
      allocate(dm1_z(nz,nx))
      dm1_x = 0.
      dm1_z = 0.
    endif
    write(0,*) "Density perturbation allocation ==> OK"

    if (.not.adj) then
      call model_delta(nz,nx,dm1,dm1_x,dm1_z,dm2,dm3)
    endif
    write(0,*) "Model perturbation initialization ==> OK"

!###################################################
! Forward background data d0
!###################################################
    !Time loop for 2nd order time
    write(0,*)"============================================================="
    write(0,*)"BEGIN OF TIME MODELING FOR BACKGROUND WAVEFIELD"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

    !Zero arrays
    Vx0       = 0.
    Vz0       = 0.
    sigmaxx0  = 0.
    sigmazz0  = 0.
    sigmaxz0  = 0.

! First time steps
    !it=1
    Vx0(:,:,1) = 0.
    Vz0(:,:,1) = 0.
    sigmaxx0(:,:,1) = 0.
    sigmazz0(:,:,1) = 0.
    sigmaxz0(:,:,1) = 0.

    !it=2
    Vx0(:,:,2) = 0.
    Vz0(:,:,2) = 0.
    call explosive_source(source(1),xsource,zsource,sigmaxx0(:,:,2))
    call explosive_source(source(1),xsource,zsource,sigmazz0(:,:,2))
    sigmaxz0(:,:,2) = 0.

    !it=3
    !Vx
    call Ax(.false.,.false.,sigmaxx0(:,:,2),Vx0(:,:,3))
    call Bx(.false.,.true.,sigmaxz0(:,:,2),Vx0(:,:,3))
    call K(Vx0(:,:,3))
    !Vz
    call Az(.false.,.false.,sigmaxz0(:,:,2),Vz0(:,:,3))
    call Bz(.false.,.true.,sigmazz0(:,:,2),Vz0(:,:,3))
    call K(Vz0(:,:,3))
    !sigmaxx
    call explosive_source(source(2),xsource,zsource,sigmaxx0(:,:,3))
    call D(.false.,.true.,Vx0(:,:,2),sigmaxx0(:,:,3))
    call E(.false.,.true.,Vz0(:,:,2),sigmaxx0(:,:,3))
    call K(sigmaxx0(:,:,3))
    !sigmazz
    call explosive_source(source(2),xsource,zsource,sigmazz0(:,:,3))
    call F(.false.,.true.,Vx0(:,:,2),sigmazz0(:,:,3))
    call G(.false.,.true.,Vz0(:,:,2),sigmazz0(:,:,3))
    call K(sigmazz0(:,:,3))
    !sigmaxz
    call H(.false.,.false.,Vx0(:,:,2),sigmaxz0(:,:,3))
    call J(.false.,.true.,Vz0(:,:,2),sigmaxz0(:,:,3))
    call K(sigmaxz0(:,:,3))

    !it=3 until it=nt
    do it=4,nt
      !Vx
      Vx0(:,:,it) = Vx0(:,:,(it-2))
      call Ax(.false.,.true.,sigmaxx0(:,:,(it-1)),Vx0(:,:,it))
      call Bx(.false.,.true.,sigmaxz0(:,:,(it-1)),Vx0(:,:,it))
      call K(Vx0(:,:,it))
      !Vz
      Vz0(:,:,it) = Vz0(:,:,(it-2))
      call Az(.false.,.true.,sigmaxz0(:,:,(it-1)),Vz0(:,:,it))
      call Bz(.false.,.true.,sigmazz0(:,:,(it-1)),Vz0(:,:,it))
      call K(Vz0(:,:,it))
      !sigmaxx
      sigmaxx0(:,:,it) = sigmaxx0(:,:,(it-2))
      call explosive_source(source(it-1),xsource,zsource,sigmaxx0(:,:,it))
      call D(.false.,.true.,Vx0(:,:,(it-1)),sigmaxx0(:,:,it))
      call E(.false.,.true.,Vz0(:,:,(it-1)),sigmaxx0(:,:,it))
      call K(sigmaxx0(:,:,it))
      !sigmazz
      sigmazz0(:,:,it) = sigmazz0(:,:,(it-2))
      call explosive_source(source(it-1),xsource,zsource,sigmazz0(:,:,it))
      call F(.false.,.true.,Vx0(:,:,(it-1)),sigmazz0(:,:,it))
      call G(.false.,.true.,Vz0(:,:,(it-1)),sigmazz0(:,:,it))
      call K(sigmazz0(:,:,it))
      !sigmaxz_x
      sigmaxz0(:,:,it) = sigmaxz0(:,:,(it-2))
      call H(.false.,.true.,Vx0(:,:,(it-1)),sigmaxz0(:,:,it))
      call J(.false.,.true.,Vz0(:,:,(it-1)),sigmaxz0(:,:,it))
      call K(sigmaxz0(:,:,it))

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif

    enddo
    write(0,'(a)',advance="yes")"]"

!###################################################
! Initialization of the imaging condition operator
!###################################################

    call image_init(nx,dx,nz,dz,nt,dt,Vx0,Vz0,lame)
    if (.not.lame) call image_init2(m1,m2,m3)

!###################################################
! BORN MODELING
!###################################################
    if(.not.adj) then

    call image_op2(.false.,.false.,dm1,dm2,dm3,Vx0,Vz0,sigmaxx0,sigmazz0,sigmaxz0)

    write(0,*)"============================================================="
    write(0,*)"BEGIN OF TIME MODELING FOR SCATTERED WAVEFIELD"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

    !Zero arrays
    !Vx0       = 0.
    !Vz0       = 0.
    !sigmaxx0  = 0.
    !sigmazz0  = 0.
    !sigmaxz0  = 0.

! First time steps
    !it=1
    call K(Vx0(:,:,1))
    call K(Vz0(:,:,1))
    call K(sigmaxx0(:,:,1))
    call K(sigmazz0(:,:,1))
    call K(sigmaxz0(:,:,1))

    !it=2
    !Since data at it=1 is zero, some terms are supressed
    !Vx
    call Ax(.false.,.true.,sigmaxx0(:,:,1),Vx0(:,:,2))
    call Bx(.false.,.true.,sigmaxz0(:,:,1),Vx0(:,:,2))
    call K(Vx0(:,:,2))
    !Vz
    call Az(.false.,.true.,sigmaxz0(:,:,1),Vz0(:,:,2))
    call Bz(.false.,.true.,sigmazz0(:,:,1),Vz0(:,:,2))
    call K(Vz0(:,:,2))
    !sigmaxx
    call D(.false.,.true.,Vx0(:,:,1),sigmaxx0(:,:,2))
    call E(.false.,.true.,Vz0(:,:,1),sigmaxx0(:,:,2))
    call K(sigmaxx0(:,:,2))
    !sigmazz
    call F(.false.,.true.,Vx0(:,:,1),sigmazz0(:,:,2))
    call G(.false.,.true.,Vz0(:,:,1),sigmazz0(:,:,2))
    call K(sigmazz0(:,:,2))
    !sigmaxz
    call H(.false.,.true.,Vx0(:,:,1),sigmaxz0(:,:,2))
    call J(.false.,.true.,Vz0(:,:,1),sigmaxz0(:,:,2))
    call K(sigmaxz0(:,:,2))

    !it=2 until it=nt
    do it=3,nt
      !Vx
      Vx0(:,:,it) = Vx0(:,:,it) + Vx0(:,:,(it-2))
      call Ax(.false.,.true.,sigmaxx0(:,:,(it-1)),Vx0(:,:,it))
      call Bx(.false.,.true.,sigmaxz0(:,:,(it-1)),Vx0(:,:,it))
      call K(Vx0(:,:,it))
      !Vz
      Vz0(:,:,it) =  Vz0(:,:,it) + Vz0(:,:,(it-2))
      call Az(.false.,.true.,sigmaxz0(:,:,(it-1)),Vz0(:,:,it))
      call Bz(.false.,.true.,sigmazz0(:,:,(it-1)),Vz0(:,:,it))
      call K(Vz0(:,:,it))
      !sigmaxx
      sigmaxx0(:,:,it) = sigmaxx0(:,:,it) + sigmaxx0(:,:,(it-2))
      call D(.false.,.true.,Vx0(:,:,(it-1)),sigmaxx0(:,:,it))
      call E(.false.,.true.,Vz0(:,:,(it-1)),sigmaxx0(:,:,it))
      call K(sigmaxx0(:,:,it))
      !sigmazz
      sigmazz0(:,:,it) = sigmazz0(:,:,it) + sigmazz0(:,:,(it-2))
      call F(.false.,.true.,Vx0(:,:,(it-1)),sigmazz0(:,:,it))
      call G(.false.,.true.,Vz0(:,:,(it-1)),sigmazz0(:,:,it))
      call K(sigmazz0(:,:,it))
      !sigmaxz
      sigmaxz0(:,:,it) = sigmaxz0(:,:,it) + sigmaxz0(:,:,(it-2))
      call H(.false.,.true.,Vx0(:,:,(it-1)),sigmaxz0(:,:,it))
      call J(.false.,.true.,Vz0(:,:,(it-1)),sigmaxz0(:,:,it))
      call K(sigmaxz0(:,:,it))

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif

    enddo
    write(0,'(a)',advance="yes")"]"
!###################################################
! Truncate scattered data to receiver positions
!###################################################
    Vx = 0.
    Vz = 0.
    sigmaxx = 0.
    sigmazz = 0.
    sigmaxz = 0.

    do ix=1,nx
      do it=1,nt
        Vx(it,ix) = Vx0(zrec,ix,it)
        Vz(it,ix) = Vz0(zrec,ix,it) 
        sigmaxx(it,ix) = sigmaxx0(zrec,ix,it)
        sigmazz(it,ix) = sigmazz0(zrec,ix,it)
        sigmaxz(it,ix) = sigmaxz0(zrec,ix,it)
      enddo
    enddo
!###################################################
! RTM modeling
!###################################################
    else

    write(0,*)"============================================================="
    write(0,*)"BEGIN OF ADJOINT PROPAGATION OF RECEIVER DATA"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

    !Zero arrays
    Vx0       = 0.
    Vz0       = 0.
    sigmaxx0  = 0.
    sigmazz0  = 0.
    sigmaxz0  = 0.

    ! Last time steps
    !it=nt
    Vx0(zrec,:,nt) = Vx(nt,:)
    Vz0(zrec,:,nt) = Vz(nt,:)
    sigmaxx0(zrec,:,nt) = sigmaxx(nt,:)
    sigmazz0(zrec,:,nt) = sigmazz(nt,:)
    sigmaxz0(zrec,:,nt) = sigmaxz(nt,:)

    !it=nt-1
!    call K(Vx0(:,:,nt))
!    call K(Vz0(:,:,nt))
!    call K(sigmaxx0(:,:,nt))
!    call K(sigmazz0(:,:,nt))
!    call K(sigmaxz0(:,:,nt))
!
    !Vx
    Vx0(zrec,:,nt-1) = Vx((nt-1),:)
    call D(.true.,.true.,Vx0(:,:,(nt-1)),sigmaxx0(:,:,nt))    
    call F(.true.,.true.,Vx0(:,:,(nt-1)),sigmazz0(:,:,nt))
    call H(.true.,.true.,Vx0(:,:,(nt-1)),sigmaxz0(:,:,nt))
    !Vz
    Vz0(zrec,:,nt-1) = Vz((nt-1),:)
    call E(.true.,.true.,Vz0(:,:,(nt-1)),sigmaxx0(:,:,nt))    
    call G(.true.,.true.,Vz0(:,:,(nt-1)),sigmazz0(:,:,nt))
    call J(.true.,.true.,Vz0(:,:,(nt-1)),sigmaxz0(:,:,nt))
    !sigmaxx
    sigmaxx0(zrec,:,nt-1) = sigmaxx((nt-1),:)
    call Ax(.true.,.true.,sigmaxx0(:,:,(nt-1)),Vx0(:,:,nt))
    !sigmazz
    sigmazz0(zrec,:,nt-1) = sigmazz((nt-1),:)
    call Bz(.true.,.true.,sigmazz0(:,:,(nt-1)),Vz0(:,:,nt))
    !sigmaxz
    sigmaxz0(zrec,:,nt-1) = sigmaxz((nt-1),:)
    call Bx(.true.,.true.,sigmaxz0(:,:,(nt-1)),Vx0(:,:,nt))
    call Az(.true.,.true.,sigmaxz0(:,:,(nt-1)),Vz0(:,:,nt))

    !it=nt-2 until it=1
    do it=nt-2,1,-1

!      call K(Vx0(:,:,it+1))
!      call K(Vz0(:,:,it+1))
!      call K(sigmaxx0(:,:,it+1))
!      call K(sigmazz0(:,:,it+1))
!      call K(sigmaxz0(:,:,it+1))
!
!      call K(Vx0(:,:,it+2))
!      call K(Vz0(:,:,it+2))
!      call K(sigmaxx0(:,:,it+2))
!      call K(sigmazz0(:,:,it+2))
!      call K(sigmaxz0(:,:,it+2))
!
      !Vx
      Vx0(zrec,:,it) = Vx(it,:)
      Vz0(zrec,:,it) = Vz(it,:)
      sigmaxx0(zrec,:,it) = sigmaxx(it,:)
      sigmazz0(zrec,:,it) = sigmazz(it,:)
      sigmaxz0(zrec,:,it) = sigmaxz(it,:)



      Vx0(:,:,it) =  Vx0(:,:,it) + Vx0(:,:,(it+2))
      Vz0(:,:,it) = Vz0(:,:,it) + Vz0(:,:,(it+2))
      sigmaxx0(:,:,it) = sigmaxx0(:,:,it) + sigmaxx0(:,:,(it+2))
      sigmazz0(:,:,it) = sigmazz0(:,:,it) + sigmazz0(:,:,(it+2))
      sigmaxz0(:,:,it) = sigmaxz0(:,:,it) + sigmaxz0(:,:,(it+2))



      call D(.true.,.true.,Vx0(:,:,it),sigmaxx0(:,:,(it+1)))    
      call F(.true.,.true.,Vx0(:,:,it),sigmazz0(:,:,(it+1)))
      call H(.true.,.true.,Vx0(:,:,it),sigmaxz0(:,:,(it+1)))
      !Vz
      call E(.true.,.true.,Vz0(:,:,it),sigmaxx0(:,:,(it+1)))    
      call G(.true.,.true.,Vz0(:,:,it),sigmazz0(:,:,(it+1)))
      call J(.true.,.true.,Vz0(:,:,it),sigmaxz0(:,:,(it+1)))
      !sigmaxx
      call Ax(.true.,.true.,sigmaxx0(:,:,it),Vx0(:,:,(it+1)))
      !sigmazz
      call Bz(.true.,.true.,sigmazz0(:,:,it),Vz0(:,:,(it+1)))
      !sigmaxz
      call Bx(.true.,.true.,sigmaxz0(:,:,it),Vx0(:,:,(it+1)))
      call Az(.true.,.true.,sigmaxz0(:,:,it),Vz0(:,:,(it+1)))

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif

    enddo

    write(0,'(a)',advance="yes")"]"

!###################################################
! Apply imaging condition
!###################################################
    write(0,*)"============================================================="
    write(0,*)"Imaging condition"
    write(0,*)"============================================================="

    call image_op2(.true.,.false.,dm1,dm2,dm3,Vx0,Vz0,sigmaxx0,sigmazz0,sigmaxz0)

    !call snap_movie(nz,dz,nx,dx,nt,dt,abc,Vx0t,Vz0t,Vx0dx,Vz0dz,Vx0dz,Vz0dx,Vx0,Vz0,sigmaxx0,sigmazz0,sigmaxz0)
    endif

  end subroutine

end module
