! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 12th, 2015

! RTM code that is supposed to take a forward/adjoint
! input, so that it can be easily tested/put into an
! iterative solver.

module elastic_mod

  use model_mod
  use elasticlib_mod
  use sep
  use source_mod
  use abc_mod
  use snap_mod_rtm

  implicit none
  integer, private :: nt, nx, nz, xsource, zsource, abc
  real, private    :: dt, dx, dz, fp, aux1
  logical, private :: free_surface
  integer, private :: zrec, xrecbeg, xrecend, xreclen

  real, pointer, dimension(:,:),private :: lambda,mu,rho
  real, allocatable, dimension(:,:), private :: rhox,rhoz,l2mu
  real, allocatable, dimension(:), private   :: source

  !Data
  real, pointer, dimension(:,:),private :: Vx_rec,Vz_rec,sigmaxx_rec,sigmazz_rec,sigmaxz_rec

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine elastic_init(nz_in,dz_in,nx_in,dx_in,dt_in,xsource_in,zsource_in,zrec_in,fp_in,nt_in,abc_in,free_surface_in,lambda_in,mu_in,rho_in)

    integer             :: nz_in,nx_in, xsource_in, zsource_in, nt_in, abc_in,zrec_in
    real                :: dz_in,dx_in,dt_in, fp_in
    logical             :: free_surface_in
    !Background properties
    real, dimension(:,:),target           :: lambda_in,mu_in,rho_in

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

    !Input background model
    lambda => lambda_in
    mu => mu_in
    rho => rho_in

    !Allocate lame parameters
    allocate(l2mu(nz,nx))
    allocate(rhox(nz,nx))
    allocate(rhoz(nz,nx))

    !Allocate source
    allocate(source(nt))

    !Initialize values
    !Background properties
    l2mu     =0.0
    rhox     =0.0
    rhoz     =0.0

    !Source
    source   =0.0

    !Initialize subroutines
    call model_op(nz,nx,abc,dx,dt,lambda,mu,rho,l2mu,rhox,rhoz,free_surface)
    call source_init(fp,dt,dx,250,nz,nx,nt)
    call source_op2(source)

  end subroutine

  ! Function that is called by Main
  function elastic_op(Vx_rec_in,Vz_rec_in,sigmaxx_rec_in,sigmazz_rec_in,sigmaxz_rec_in) result(stat)
    real,dimension(:,:),target  :: Vx_rec_in,Vz_rec_in,sigmaxx_rec_in,sigmazz_rec_in,sigmaxz_rec_in
    integer                     :: stat

    stat=1
    call elastic_op2(Vx_rec_in,Vz_rec_in,sigmaxx_rec_in,sigmazz_rec_in,sigmaxz_rec_in)
    stat=0
  end function

  ! Subroutine used by function
  subroutine elastic_op2(Vx_rec_in,Vz_rec_in,sigmaxx_rec_in,sigmazz_rec_in,sigmaxz_rec_in)
    real,dimension(:,:),target              :: Vx_rec_in,Vz_rec_in,sigmaxx_rec_in,sigmazz_rec_in,sigmaxz_rec_in

    !General variables
    integer                                 :: it, ix, iz

    !Temporary wavefields
    real, allocatable, dimension(:,:,:)   :: Vx0,Vz0
    real, allocatable, dimension(:,:,:)   :: sigmaxx0,sigmazz0
    real, allocatable, dimension(:,:,:)   :: sigmaxz0

    Vx_rec => Vx_rec_in
    Vz_rec => Vz_rec_in
    sigmaxx_rec => sigmaxx_rec_in
    sigmazz_rec => sigmazz_rec_in
    sigmaxz_rec => sigmaxz_rec_in

    !Allocate Velocities and stresses
    !Temporary wavefields
    allocate(Vx0(nz,nx,nt))
    allocate(Vz0(nz,nx,nt))
    allocate(sigmaxx0(nz,nx,nt))
    allocate(sigmazz0(nz,nx,nt))
    allocate(sigmaxz0(nz,nx,nt))

    !Initialize derivatives
    call elasticmod_init(nz,nx,dx,dt,rhox,rhoz,lambda,mu,l2mu)

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
! Truncate scattered data to receiver positions
!###################################################
    Vx_rec = 0.
    Vz_rec = 0.
    sigmaxx_rec = 0.
    sigmazz_rec = 0.
    sigmaxz_rec = 0.

    do ix=1,nx
      do it=1,nt
        Vx_rec(it,ix) = Vx0(zrec,ix,it)
        Vz_rec(it,ix) = Vz0(zrec,ix,it) 
        sigmaxx_rec(it,ix) = sigmaxx0(zrec,ix,it)
        sigmazz_rec(it,ix) = sigmazz0(zrec,ix,it)
        sigmaxz_rec(it,ix) = sigmaxz0(zrec,ix,it)
      enddo
    enddo

    call snap_wfld(nz,dz,nx,dx,nt,dt,abc,Vx0,Vz0,sigmaxx0,sigmazz0,sigmaxz0)

  end subroutine

end module
