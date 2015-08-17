! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Jan, 22th, 2015

! Forward Born modeling program, based on the velocity-stress
! elastic modeling program

! FD elastic operator module
! This is the main subroutine that calls
! the FD derivatives, sources, etc.
module born_mod

  use source_mod
  use seis_mod_b
  use bornlib_mod
  use snap_mod_b
  use pml_mod
  use omp_lib
  !use clib_mod

  implicit none
  integer, private :: nt, nx, nz, abc, zrec, ix
  integer, private :: xsource, zsource, nsource, nsnap, nseis
  real, private    :: dt, dx, dz, fp
  !real*16, private :: source
  real, private :: source
  logical, private :: free_surface

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine elastic_init(nz_in,dz_in,nx_in,dx_in,nt_in,dt_in,xsource_in,  &
                          zsource_in,nsource_in, fp_in, nsnap_in, nseis_in,&
                          abc_in, zrec_in, free_surface_in)

    integer             :: nz_in,nx_in,nt_in,xsource_in,zsource_in,nsource_in
    integer             :: nsnap_in, nseis_in, abc_in, zrec_in
    real                :: dz_in,dx_in,dt_in, fp_in
    logical             :: free_surface_in

    nz = nz_in
    dz = dz_in
    nt = nt_in
    dt = dt_in
    nx = nx_in
    dx = dx_in
    fp = fp_in
    xsource = xsource_in
    zsource = zsource_in
    nsource = nsource_in
    nsnap   = nsnap_in
    nseis = nseis_in
    abc = abc_in
    zrec = zrec_in
    free_surface = free_surface_in

  end subroutine

  ! Function that is called by Main
  function elastic_op(vp0, vs0, rho0, dvp, dvs, drho, seis_Vx, seis_Vz, seis_P, &
                      seis_Vx_b, seis_Vz_b, seis_P_b, snap_Vx, snap_Vz, snap_P, &
                      snap_Vx_b, snap_Vz_b, snap_P_b, lame) result(stat)
    real,dimension(:,:)         :: vp0, vs0, rho0
    real,dimension(:,:)         :: dvp, dvs, drho
    real,dimension(:,:)         :: seis_Vx, seis_Vz, seis_P
    real,dimension(:,:,:)       :: snap_Vx, snap_Vz, snap_P
    real,dimension(:,:)         :: seis_Vx_b, seis_Vz_b, seis_P_b
    real,dimension(:,:,:)       :: snap_Vx_b, snap_Vz_b, snap_P_b
    integer                     :: stat
    logical,intent(in),optional :: lame

    stat=1
    call elastic_op2(vp0, vs0, rho0, dvp, dvs, drho, seis_Vx, seis_Vz, seis_P, &
                     seis_Vx_b, seis_Vz_b, seis_P_b, snap_Vx, snap_Vz, snap_P, &
                     snap_Vx_b, snap_Vz_b, snap_P_b, lame)
    stat=0
  end function

  ! Subroutine used by function
  subroutine elastic_op2(vp0, vs0, rho0, dvp, dvs, drho, seis_Vx, seis_Vz, seis_P, &
                         seis_Vx_b, seis_Vz_b, seis_P_b, snap_Vx, snap_Vz, snap_P, &
                         snap_Vx_b, snap_Vz_b, snap_P_b, lame)
    integer                             :: it, ix, iz, dsnap, dseis
    integer                             :: countsnap,countseis, temp
    real, dimension(nz,nx)              :: vp0,vs0,rho0
    real, dimension(nz,nx)              :: dvp,dvs,drho

    !Outputs
    real, dimension(nseis,nx-2*abc)     :: seis_Vx, seis_Vz, seis_P
    real, dimension(nz,nx,nsnap)        :: snap_Vx, snap_Vz, snap_P
    real, dimension(nseis,nx-2*abc)     :: seis_Vx_b, seis_Vz_b, seis_P_b
    real, dimension(nz,nx,nsnap)        :: snap_Vx_b, snap_Vz_b, snap_P_b

    !Background propagation
    real, allocatable, dimension(:,:)   :: mi0,l0,l2mi0
    real, allocatable, dimension(:,:),target   :: Vx_x1,Vx_z1,Vz_x1,Vz_z1
    real, allocatable, dimension(:,:),target   :: sigmaxx_x1,sigmaxx_z1,sigmazz_x1,sigmazz_z1
    real, allocatable, dimension(:,:),target   :: sigmaxz_x1,sigmaxz_z1
    real, allocatable, dimension(:,:),target   :: Vx_x2,Vx_z2,Vz_x2,Vz_z2
    real, allocatable, dimension(:,:),target   :: sigmaxx_x2,sigmaxx_z2,sigmazz_x2,sigmazz_z2
    real, allocatable, dimension(:,:),target   :: sigmaxz_x2,sigmaxz_z2
    real, allocatable, dimension(:,:),target   :: Vx_x3,Vx_z3,Vz_x3,Vz_z3
    real, allocatable, dimension(:,:),target   :: sigmaxx_x3,sigmaxx_z3,sigmazz_x3,sigmazz_z3
    real, allocatable, dimension(:,:),target   :: sigmaxz_x3,sigmaxz_z3
    real, pointer, dimension(:,:)       :: Vx_xnew,Vx_znew,Vz_xnew,Vz_znew
    real, pointer, dimension(:,:)       :: sigmaxx_xnew,sigmaxx_znew,sigmazz_xnew,sigmazz_znew
    real, pointer, dimension(:,:)       :: sigmaxz_xnew,sigmaxz_znew
    real, pointer, dimension(:,:)       :: Vx_xcur,Vx_zcur,Vz_xcur,Vz_zcur
    real, pointer, dimension(:,:)       :: sigmaxx_xcur,sigmaxx_zcur,sigmazz_xcur,sigmazz_zcur
    real, pointer, dimension(:,:)       :: sigmaxz_xcur,sigmaxz_zcur
    real, pointer, dimension(:,:)       :: Vx_xold,Vx_zold,Vz_xold,Vz_zold
    real, pointer, dimension(:,:)       :: sigmaxx_xold,sigmaxx_zold,sigmazz_xold,sigmazz_zold
    real, pointer, dimension(:,:)       :: sigmaxz_xold,sigmaxz_zold
    real, pointer, dimension(:,:)       :: Vx_xtemp,Vx_ztemp,Vz_xtemp,Vz_ztemp
    real, pointer, dimension(:,:)       :: sigmaxx_xtemp,sigmaxx_ztemp,sigmazz_xtemp,sigmazz_ztemp
    real, pointer, dimension(:,:)       :: sigmaxz_xtemp,sigmaxz_ztemp
    real, allocatable, dimension(:,:)   :: rho0_Vx,rho0_Vz
    real, allocatable, dimension(:,:)   :: drho_Vx,drho_Vz
    real                                :: aux1,t,nsource2

    !Scatterer propagation
    real, allocatable, dimension(:,:)   :: dmi, dl, dl2mi
    real, allocatable, dimension(:,:),target   :: Vx_x1_b,Vx_z1_b,Vz_x1_b,Vz_z1_b
    real, allocatable, dimension(:,:),target   :: sigmaxx_x1_b,sigmaxx_z1_b,sigmazz_x1_b,sigmazz_z1_b
    real, allocatable, dimension(:,:),target   :: sigmaxz_x1_b,sigmaxz_z1_b
    real, allocatable, dimension(:,:),target   :: Vx_x2_b,Vx_z2_b,Vz_x2_b,Vz_z2_b
    real, allocatable, dimension(:,:),target   :: sigmaxx_x2_b,sigmaxx_z2_b,sigmazz_x2_b,sigmazz_z2_b
    real, allocatable, dimension(:,:),target   :: sigmaxz_x2_b,sigmaxz_z2_b
    real, allocatable, dimension(:,:),target   :: Vx_x3_b,Vx_z3_b,Vz_x3_b,Vz_z3_b
    real, allocatable, dimension(:,:),target   :: sigmaxx_x3_b,sigmaxx_z3_b,sigmazz_x3_b,sigmazz_z3_b
    real, allocatable, dimension(:,:),target   :: sigmaxz_x3_b,sigmaxz_z3_b
    real, pointer, dimension(:,:)       :: Vx_xnew_b,Vx_znew_b,Vz_xnew_b,Vz_znew_b
    real, pointer, dimension(:,:)       :: sigmaxx_xnew_b,sigmaxx_znew_b,sigmazz_xnew_b,sigmazz_znew_b
    real, pointer, dimension(:,:)       :: sigmaxz_xnew_b,sigmaxz_znew_b
    real, pointer, dimension(:,:)       :: Vx_xcur_b,Vx_zcur_b,Vz_xcur_b,Vz_zcur_b
    real, pointer, dimension(:,:)       :: sigmaxx_xcur_b,sigmaxx_zcur_b,sigmazz_xcur_b,sigmazz_zcur_b
    real, pointer, dimension(:,:)       :: sigmaxz_xcur_b,sigmaxz_zcur_b
    real, pointer, dimension(:,:)       :: Vx_xold_b,Vx_zold_b,Vz_xold_b,Vz_zold_b
    real, pointer, dimension(:,:)       :: sigmaxx_xold_b,sigmaxx_zold_b,sigmazz_xold_b,sigmazz_zold_b
    real, pointer, dimension(:,:)       :: sigmaxz_xold_b,sigmaxz_zold_b
    real, pointer, dimension(:,:)       :: Vx_xtemp_b,Vx_ztemp_b,Vz_xtemp_b,Vz_ztemp_b
    real, pointer, dimension(:,:)       :: sigmaxx_xtemp_b,sigmaxx_ztemp_b,sigmazz_xtemp_b,sigmazz_ztemp_b
    real, pointer, dimension(:,:)       :: sigmaxz_xtemp_b,sigmaxz_ztemp_b
    integer                             :: stat
    logical,intent(in),optional         :: lame

    !Memory allocation

    !Allocate lame parameters
    allocate(mi0(nz,nx))
    allocate(l0(nz,nx))
    allocate(l2mi0(nz,nx))
    allocate(rho0_Vx(nz,nx))
    allocate(rho0_Vz(nz,nx))
    allocate(dmi(nz,nx))
    allocate(dl(nz,nx))
    allocate(dl2mi(nz,nx))
    allocate(drho_Vx(nz,nx))
    allocate(drho_Vz(nz,nx))

    !Allocate Velocities and stresses
    allocate(Vx_x1(nz,nx))
    allocate(Vx_z1(nz,nx))
    allocate(Vz_x1(nz,nx))
    allocate(Vz_z1(nz,nx))
    allocate(sigmaxx_x1(nz,nx))
    allocate(sigmaxx_z1(nz,nx))
    allocate(sigmazz_x1(nz,nx))
    allocate(sigmazz_z1(nz,nx))
    allocate(sigmaxz_x1(nz,nx))
    allocate(sigmaxz_z1(nz,nx))
    allocate(Vx_x2(nz,nx))
    allocate(Vx_z2(nz,nx))
    allocate(Vz_x2(nz,nx))
    allocate(Vz_z2(nz,nx))
    allocate(sigmaxx_x2(nz,nx))
    allocate(sigmaxx_z2(nz,nx))
    allocate(sigmazz_x2(nz,nx))
    allocate(sigmazz_z2(nz,nx))
    allocate(sigmaxz_x2(nz,nx))
    allocate(sigmaxz_z2(nz,nx))
    allocate(Vx_x3(nz,nx))
    allocate(Vx_z3(nz,nx))
    allocate(Vz_x3(nz,nx))
    allocate(Vz_z3(nz,nx))
    allocate(sigmaxx_x3(nz,nx))
    allocate(sigmaxx_z3(nz,nx))
    allocate(sigmazz_x3(nz,nx))
    allocate(sigmazz_z3(nz,nx))
    allocate(sigmaxz_x3(nz,nx))
    allocate(sigmaxz_z3(nz,nx))

    !Scatterer wavefield
    allocate(Vx_x1_b(nz,nx))
    allocate(Vx_z1_b(nz,nx))
    allocate(Vz_x1_b(nz,nx))
    allocate(Vz_z1_b(nz,nx))
    allocate(sigmaxx_x1_b(nz,nx))
    allocate(sigmaxx_z1_b(nz,nx))
    allocate(sigmazz_x1_b(nz,nx))
    allocate(sigmazz_z1_b(nz,nx))
    allocate(sigmaxz_x1_b(nz,nx))
    allocate(sigmaxz_z1_b(nz,nx))
    allocate(Vx_x2_b(nz,nx))
    allocate(Vx_z2_b(nz,nx))
    allocate(Vz_x2_b(nz,nx))
    allocate(Vz_z2_b(nz,nx))
    allocate(sigmaxx_x2_b(nz,nx))
    allocate(sigmaxx_z2_b(nz,nx))
    allocate(sigmazz_x2_b(nz,nx))
    allocate(sigmazz_z2_b(nz,nx))
    allocate(sigmaxz_x2_b(nz,nx))
    allocate(sigmaxz_z2_b(nz,nx))
    allocate(Vx_x3_b(nz,nx))
    allocate(Vx_z3_b(nz,nx))
    allocate(Vz_x3_b(nz,nx))
    allocate(Vz_z3_b(nz,nx))
    allocate(sigmaxx_x3_b(nz,nx))
    allocate(sigmaxx_z3_b(nz,nx))
    allocate(sigmazz_x3_b(nz,nx))
    allocate(sigmazz_z3_b(nz,nx))
    allocate(sigmaxz_x3_b(nz,nx))
    allocate(sigmaxz_z3_b(nz,nx))

    !Initialize values
    !Background properties
    mi0         =0.0
    l0          =0.0
    l2mi0       =0.0
    rho0_Vx     =0.0
    rho0_Vz     =0.0

    !Perturbation properties
    dmi         =0.0
    dl          =0.0
    dl2mi       =0.0
    drho_Vx     =0.0
    drho_Vz     =0.0

    !Background wavefield
    Vx_x1       =0.0
    Vx_z1       =0.0
    Vz_x1       =0.0
    Vz_z1       =0.0
    sigmaxx_x1  =0.0
    sigmaxx_z1  =0.0
    sigmazz_x1  =0.0
    sigmazz_z1  =0.0
    sigmaxz_x1  =0.0
    sigmaxz_z1  =0.0
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
    Vx_x3       =0.0
    Vx_z3       =0.0
    Vz_x3       =0.0
    Vz_z3       =0.0
    sigmaxx_x3  =0.0
    sigmaxx_z3  =0.0
    sigmazz_x3  =0.0
    sigmazz_z3  =0.0
    sigmaxz_x3  =0.0
    sigmaxz_z3  =0.0

    !Scatterer wavefield
    Vx_x1_b       =0.0
    Vx_z1_b       =0.0
    Vz_x1_b       =0.0
    Vz_z1_b       =0.0
    sigmaxx_x1_b  =0.0
    sigmaxx_z1_b  =0.0
    sigmazz_x1_b  =0.0
    sigmazz_z1_b  =0.0
    sigmaxz_x1_b  =0.0
    sigmaxz_z1_b  =0.0
    Vx_x2_b       =0.0
    Vx_z2_b       =0.0
    Vz_x2_b       =0.0
    Vz_z2_b       =0.0
    sigmaxx_x2_b  =0.0
    sigmaxx_z2_b  =0.0
    sigmazz_x2_b  =0.0
    sigmazz_z2_b  =0.0
    sigmaxz_x2_b  =0.0
    sigmaxz_z2_b  =0.0
    Vx_x3_b       =0.0
    Vx_z3_b       =0.0
    Vz_x3_b       =0.0
    Vz_z3_b       =0.0
    sigmaxx_x3_b  =0.0
    sigmaxx_z3_b  =0.0
    sigmazz_x3_b  =0.0
    sigmazz_z3_b  =0.0
    sigmaxz_x3_b  =0.0
    sigmaxz_z3_b  =0.0

    dseis      =(nt/nseis)+1
    dsnap      =nt/nsnap
    countsnap  =1
    countseis  =1
    nsource2   =2*nint(sqrt(6.0)/3.14159265/fp/dt)

    !Define Lame parameters
    !These model properties are still in the original grid, centered at i,j
    do ix=1,nx
      do iz=1,nz
        if(lame)then
          l0(iz,ix)=vp0(iz,ix)
          mi0(iz,ix)=vs0(iz,ix)

          dl(iz,ix)=dvp(iz,ix)
          dmi(iz,ix)=dvs(iz,ix)
        else
          mi0(iz,ix)=vs0(iz,ix)*vs0(iz,ix)*rho0(iz,ix)
          l0(iz,ix)=vp0(iz,ix)*vp0(iz,ix)*rho0(iz,ix)-2.0*mi0(iz,ix)
          l2mi0(iz,ix)=l0(iz,ix)+2.0*mi0(iz,ix)

          dmi(iz,ix)=dvs(iz,ix)*dvs(iz,ix)*drho(iz,ix)
          dl(iz,ix)=dvp(iz,ix)*dvp(iz,ix)*drho(iz,ix)-2.0*dmi(iz,ix)
          dl2mi(iz,ix)=dl(iz,ix)+2.0*dmi(iz,ix)
        endif
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
        drho_Vx(iz,ix)=-(drho(iz+1,ix)/rho0(iz+1,ix) +&
                         drho(iz  ,ix)/rho0(iz  ,ix))
      enddo
    enddo
    !Rho_vz
    do ix=1,nx-1
      do iz=1,nz
        rho0_Vz(iz,ix)=aux1*1.0/(0.5*(rho0(iz,ix+1) + rho0(iz,ix)))
        drho_Vz(iz,ix)=-(drho(iz,ix+1)/rho0(iz,ix+1) +&
                         drho(iz,ix  )/rho0(iz,ix  ))
      enddo
    enddo
    !Right side
    do iz=1,nz
      rho0_Vz(iz,nx)=aux1*(1.0/rho0(iz,nx))
      drho_Vz(iz,nx)=-drho(iz,nx)/rho0(iz,nx)
    enddo
    !Bottom
    do ix=1,nx
      rho0_Vx(nz,ix)=aux1*(1.0/rho0(nz,ix))
      drho_Vx(nz,ix)=-drho(nz,ix)/rho0(nz,ix)
    enddo
    
    !Calculate weigthed average for the lambda variable
    !These put lambda and (lambda+2mi) in their respective staggered grids,
    !using the harmonic average proposed in Moczo, 2002
    do ix=2,nx-1
      do iz=2,nz-1
        l0(iz,ix)=aux1*4.0/((1.0/l0(iz,ix))+(1.0/l0(iz+1,ix))+(1.0/l0(iz,ix+1))+(1.0/l0(iz+1,ix+1)))
        if (mi0(iz,ix)==0.or.mi0(iz+1,ix)==0.or.mi0(iz,ix+1)==0.or.mi0(iz+1,ix+1)==0) then
          mi0(iz,ix)=0.0
        else
          mi0(iz,ix)=aux1*4.0/((1.0/mi0(iz,ix))+(1.0/mi0(iz+1,ix))+(1.0/mi0(iz,ix+1))+(1.0/mi0(iz+1,ix+1)))
        endif
        if (dl(iz,ix)<1e-2.or.dl(iz+1,ix)<1e-2.or.dl(iz,ix+1)<1e-2.or.dl(iz+1,ix+1)<1e-2) then
          dl(iz,ix)=0.0
        else
          dl(iz,ix)=aux1*4.0/((1.0/dl(iz,ix))+(1.0/dl(iz+1,ix))+(1.0/dl(iz,ix+1))+(1.0/dl(iz+1,ix+1)))
        endif
        if (dmi(iz,ix)<1e-2.or.dmi(iz+1,ix)<1e-2.or.dmi(iz,ix+1)<1e-2.or.dmi(iz+1,ix+1)<1e-2) then
          dmi(iz,ix)=0.0
        else
          dmi(iz,ix)=aux1*4.0/((1.0/dmi(iz,ix))+(1.0/dmi(iz+1,ix))+(1.0/dmi(iz,ix+1))+(1.0/dmi(iz+1,ix+1)))
        endif
        l2mi0(iz,ix)=l0(iz,ix)+2.0*mi0(iz,ix)
        dl2mi(iz,ix)=dl(iz,ix)+2.0*dmi(iz,ix)
      enddo
    enddo

    !Initialize subroutines
    call source_init(fp,dt,dx,nsource,nz,nx)
    call seis_init(nz,nx, abc, zrec, seis_Vx, seis_Vz, seis_P, seis_Vx_b, seis_Vz_b, seis_P_b)
    call bornmod_init(nz,nx,rho0_Vx,rho0_Vz,l0,mi0,l2mi0,drho_vx,drho_vz,dl,dmi,dl2mi)
    call pml_init(nz,nx,dx,dt,abc,vp0,free_surface)
    call snap_init(nz,nx,abc,snap_Vx,snap_Vz, snap_P, snap_Vx_b, snap_Vz_b, snap_P_b)
    
    !Time loop for 2nd order time
    write(0,*)"============================================================="
    write(0,*)"BEGIN OF TIME MODELING"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

    !Initialize pointers
    Vx_xnew => Vx_x1
    Vx_znew => Vx_z1
    Vz_xnew => Vz_x1
    Vz_znew => Vz_z1
    sigmaxx_xnew => sigmaxx_x1
    sigmaxx_znew => sigmaxx_z1
    sigmazz_xnew => sigmazz_x1
    sigmazz_znew => sigmazz_z1
    sigmaxz_xnew => sigmaxz_x1
    sigmaxz_znew => sigmaxz_z1
    Vx_xcur => Vx_x2
    Vx_zcur => Vx_z2
    Vz_xcur => Vz_x2
    Vz_zcur => Vz_z2
    sigmaxx_xcur => sigmaxx_x2
    sigmaxx_zcur => sigmaxx_z2
    sigmazz_xcur => sigmazz_x2
    sigmazz_zcur => sigmazz_z2
    sigmaxz_xcur => sigmaxz_x2
    sigmaxz_zcur => sigmaxz_z2
    Vx_xold => Vx_x3
    Vx_zold => Vx_z3
    Vz_xold => Vz_x3
    Vz_zold => Vz_z3
    sigmaxx_xold => sigmaxx_x3
    sigmaxx_zold => sigmaxx_z3
    sigmazz_xold => sigmazz_x3
    sigmazz_zold => sigmazz_z3
    sigmaxz_xold => sigmaxz_x3
    sigmaxz_zold => sigmaxz_z3
    Vx_xnew_b => Vx_x1_b
    Vx_znew_b => Vx_z1_b
    Vz_xnew_b => Vz_x1_b
    Vz_znew_b => Vz_z1_b
    sigmaxx_xnew_b => sigmaxx_x1_b
    sigmaxx_znew_b => sigmaxx_z1_b
    sigmazz_xnew_b => sigmazz_x1_b
    sigmazz_znew_b => sigmazz_z1_b
    sigmaxz_xnew_b => sigmaxz_x1_b
    sigmaxz_znew_b => sigmaxz_z1_b
    Vx_xcur_b => Vx_x2_b
    Vx_zcur_b => Vx_z2_b
    Vz_xcur_b => Vz_x2_b
    Vz_zcur_b => Vz_z2_b
    sigmaxx_xcur_b => sigmaxx_x2_b
    sigmaxx_zcur_b => sigmaxx_z2_b
    sigmazz_xcur_b => sigmazz_x2_b
    sigmazz_zcur_b => sigmazz_z2_b
    sigmaxz_xcur_b => sigmaxz_x2_b
    sigmaxz_zcur_b => sigmaxz_z2_b
    Vx_xold_b => Vx_x3_b
    Vx_zold_b => Vx_z3_b
    Vz_xold_b => Vz_x3_b
    Vz_zold_b => Vz_z3_b
    sigmaxx_xold_b => sigmaxx_x3_b
    sigmaxx_zold_b => sigmaxx_z3_b
    sigmazz_xold_b => sigmazz_x3_b
    sigmazz_zold_b => sigmazz_z3_b
    sigmaxz_xold_b => sigmaxz_x3_b
    sigmaxz_zold_b => sigmaxz_z3_b

    do it=1,nt
      if (it<nsource) then
        call source_op2(it, source)
      else
        source=0.0
      endif
      source = 0.7E7*source*(dt/(dx*dx))*0.25

!###################################################
! Propagation of the wavefields in the background
!###################################################
      !Vx_x
      call Ax(.false.,.false.,sigmaxx_xcur,Vx_xnew)
      call Ax(.false.,.true.,sigmaxx_zcur,Vx_xnew)
      call L(Vx_xold)
      Vx_xnew = Vx_xnew + Vx_xold
      !call explosive_source(.false.,.false.,source,xsource,zsource,temp1)
      !call C(.false.,.true.,temp1,Vx_xnew)
      call K(Vx_xnew)

      !Vx_z
      call Bx(.false.,.false.,sigmaxz_xcur,Vx_znew)
      call Bx(.false.,.true.,sigmaxz_zcur,Vx_znew)
      call N(Vx_zold)
      Vx_znew = Vx_znew + Vx_zold
      !call explosive_source(.false.,.false.,source,xsource,zsource,temp1)
      !call C(.false.,.true.,temp1,Vx_znew)
      call M(Vx_znew)

      !Vz_x
      call Az(.false.,.false.,sigmaxz_xcur,Vz_xnew)
      call Az(.false.,.true.,sigmaxz_zcur,Vz_xnew)
      call P(Vz_xold)
      Vz_xnew = Vz_xnew + Vz_xold
      !call explosive_source(.false.,.false.,source,xsource,zsource,temp1)
      !call C(.false.,.true.,temp1,Vz_xnew)
      call O(Vz_xnew)

      !Vz_z
      call Bz(.false.,.false.,sigmazz_xcur,Vz_znew)
      call Bz(.false.,.true.,sigmazz_zcur,Vz_znew)
      call R(Vz_zold)
      Vz_znew = Vz_znew + Vz_zold
      !call explosive_source(.false.,.false.,source,xsource,zsource,temp1)
      !call C(.false.,.true.,temp1,Vz_znew)
      call Q(Vz_znew)

      !sigmaxx_x
      call D(.false.,.false.,Vx_xcur,sigmaxx_xnew)
      call D(.false.,.true.,Vx_zcur,sigmaxx_xnew)
      call P(sigmaxx_xold)
      sigmaxx_xnew = sigmaxx_xnew + sigmaxx_xold
      call explosive_sourcex(.false.,.true.,source,xsource,zsource,sigmaxx_xnew)
      call O(sigmaxx_xnew)

      !sigmaxx_z
      call E(.false.,.false.,Vz_xcur,sigmaxx_znew)
      call E(.false.,.true.,Vz_zcur,sigmaxx_znew)
      call R(sigmaxx_zold)
      sigmaxx_znew = sigmaxx_znew + sigmaxx_zold
      call explosive_sourcez(.false.,.true.,source,xsource,zsource,sigmaxx_znew)
      call Q(sigmaxx_znew)

      !sigmazz_x
      call F(.false.,.false.,Vx_xcur,sigmazz_xnew)
      call F(.false.,.true.,Vx_zcur,sigmazz_xnew)
      call P(sigmazz_xold)
      sigmazz_xnew = sigmazz_xnew + sigmazz_xold
      call explosive_sourcex(.false.,.true.,source,xsource,zsource,sigmazz_xnew)
      call O(sigmazz_xnew)

      !sigmazz_z
      call G(.false.,.false.,Vz_xcur,sigmazz_znew)
      call G(.false.,.true.,Vz_zcur,sigmazz_znew)
      call R(sigmazz_zold)
      sigmazz_znew = sigmazz_znew + sigmazz_zold
      call explosive_sourcez(.false.,.true.,source,xsource,zsource,sigmazz_znew)
      call Q(sigmazz_znew)

      !sigmaxz_x
      call H(.false.,.false.,Vx_xcur,sigmaxz_xnew)
      call H(.false.,.true.,Vx_zcur,sigmaxz_xnew)
      call L(sigmaxz_xold)
      sigmaxz_xnew = sigmaxz_xnew + sigmaxz_xold
      !call explosive_source(.false.,.true.,source,xsource,zsource,sigmaxz_xnew)
      call K(sigmaxz_xnew)

      !sigmaxz_z
      call J(.false.,.false.,Vz_xcur,sigmaxz_znew)
      call J(.false.,.true.,Vz_zcur,sigmaxz_znew)
      call N(sigmaxz_zold)
      sigmaxz_znew = sigmaxz_znew + sigmaxz_zold
      !call explosive_source(.false.,.true.,source,xsource,zsource,sigmaxz_znew)
      call M(sigmaxz_znew)

!###################################################
! Propagation of the scattered wavefield
!###################################################
      !Vx_x
      call Ax(.false.,.false.,sigmaxx_xcur_b,Vx_xnew_b)
      call Ax(.false.,.true.,sigmaxx_zcur_b,Vx_xnew_b)
      call L(Vx_xold_b)
      Vx_xnew_b = Vx_xnew_b + Vx_xold_b
      call S1(.false.,.true.,Vx_xnew,Vx_xcur,Vx_xnew_b)
      !call C(.false.,.true.,temp1,Vx_xnew_b)
      call K(Vx_xnew_b)

      !Vx_z
      call Bx(.false.,.false.,sigmaxz_xcur_b,Vx_znew_b)
      call Bx(.false.,.true.,sigmaxz_zcur_b,Vx_znew_b)
      call N(Vx_zold_b)
      Vx_znew_b = Vx_znew_b + Vx_zold_b
      call S1(.false.,.true.,Vx_znew,Vx_zcur,Vx_znew_b)
      !call C(.false.,.true.,temp1,Vx_znew_b)
      call M(Vx_znew_b)

      !Vz_x
      call Az(.false.,.false.,sigmaxz_xcur_b,Vz_xnew_b)
      call Az(.false.,.true.,sigmaxz_zcur_b,Vz_xnew_b)
      call P(Vz_xold_b)
      Vz_xnew_b = Vz_xnew_b + Vz_xold_b
      call S2(.false.,.true.,Vz_xnew,Vz_xcur,Vz_xnew_b)
      !call C(.false.,.true.,temp1,Vz_xnew_b)
      call O(Vz_xnew_b)

      !Vz_z
      call Bz(.false.,.false.,sigmazz_xcur_b,Vz_znew_b)
      call Bz(.false.,.true.,sigmazz_zcur_b,Vz_znew_b)
      call R(Vz_zold_b)
      Vz_znew_b = Vz_znew_b + Vz_zold_b
      call S2(.false.,.true.,Vz_znew,Vz_zcur,Vz_znew_b)
      !call C(.false.,.true.,temp1,Vz_znew_b)
      call Q(Vz_znew_b)

      !sigmaxx_x
      call D(.false.,.false.,Vx_xcur_b,sigmaxx_xnew_b)
      call D(.false.,.true.,Vx_zcur_b,sigmaxx_xnew_b)
      call P(sigmaxx_xold_b)
      sigmaxx_xnew_b = sigmaxx_xnew_b + sigmaxx_xold_b
      call S3(.false.,.true.,Vx_xnew,sigmaxx_xnew_b)
      call S3(.false.,.true.,Vx_znew,sigmaxx_xnew_b)
      call O(sigmaxx_xnew_b)

      !sigmaxx_z
      call E(.false.,.false.,Vz_xcur_b,sigmaxx_znew_b)
      call E(.false.,.true.,Vz_zcur_b,sigmaxx_znew_b)
      call R(sigmaxx_zold_b)
      sigmaxx_znew_b = sigmaxx_znew_b + sigmaxx_zold_b
      call S4(.false.,.true.,Vz_xnew,sigmaxx_znew_b)
      call S4(.false.,.true.,Vz_znew,sigmaxx_znew_b)
      call Q(sigmaxx_znew_b)

      !sigmazz_x
      call F(.false.,.false.,Vx_xcur_b,sigmazz_xnew_b)
      call F(.false.,.true.,Vx_zcur_b,sigmazz_xnew_b)
      call P(sigmazz_xold_b)
      sigmazz_xnew_b = sigmazz_xnew_b + sigmazz_xold_b
      call S5(.false.,.true.,Vx_xnew,sigmazz_xnew_b)
      call S5(.false.,.true.,Vx_znew,sigmazz_xnew_b)
      call O(sigmazz_xnew_b)

      !sigmazz_z
      call G(.false.,.false.,Vz_xcur_b,sigmazz_znew_b)
      call G(.false.,.true.,Vz_zcur_b,sigmazz_znew_b)
      call R(sigmazz_zold_b)
      sigmazz_znew_b = sigmazz_znew_b + sigmazz_zold_b
      call S6(.false.,.true.,Vz_xnew,sigmazz_znew_b)
      call S6(.false.,.true.,Vz_znew,sigmazz_znew_b)
      call Q(sigmazz_znew_b)

      !sigmaxz_x
      call H(.false.,.false.,Vx_xcur_b,sigmaxz_xnew_b)
      call H(.false.,.true.,Vx_zcur_b,sigmaxz_xnew_b)
      call L(sigmaxz_xold_b)
      sigmaxz_xnew_b = sigmaxz_xnew_b + sigmaxz_xold_b
      call S7(.false.,.true.,Vx_xnew,sigmaxz_xnew_b)
      call S7(.false.,.true.,Vx_znew,sigmaxz_xnew_b)
      call K(sigmaxz_xnew_b)

      !sigmaxz_z
      call J(.false.,.false.,Vz_xcur_b,sigmaxz_znew_b)
      call J(.false.,.true.,Vz_zcur_b,sigmaxz_znew_b)
      call N(sigmaxz_zold_b)
      sigmaxz_znew_b = sigmaxz_znew_b + sigmaxz_zold_b
      call S8(.false.,.true.,Vz_xnew,sigmaxz_znew_b)
      call S8(.false.,.true.,Vz_znew,sigmaxz_znew_b)
      call M(sigmaxz_znew_b)

!###################################################
! Outputs and time updates
!###################################################

      !Write seismograms
      if (it.gt.nsource2.and.mod(it,dseis)==0) then
        call seis_op2(countseis,Vx_xnew,Vx_znew,Vz_xnew,Vz_znew)
        call seis_op3(countseis,sigmaxx_xnew,sigmaxx_znew,sigmazz_xnew,sigmazz_znew)
        call seis_op4(countseis,Vx_xnew_b,Vx_znew_b,Vz_xnew_b,Vz_znew_b)
        call seis_op5(countseis,sigmaxx_xnew_b,sigmaxx_znew_b,sigmazz_xnew_b,sigmazz_znew_b)
        countseis = countseis + 1
        !write(0,*)"Seismogram time sample"
      endif

      !Write snapshots
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        call snap_op2(countsnap,Vx_xnew,Vx_znew,Vz_xnew,Vz_znew,Vx_xnew_b,Vx_znew_b,Vz_xnew_b,Vz_znew_b)
        call snap_op3(countsnap,sigmaxx_xnew,sigmaxx_znew,sigmazz_xnew,sigmazz_znew,sigmaxx_xnew_b,sigmaxx_znew_b,sigmazz_xnew_b,sigmazz_znew_b)
        countsnap = countsnap + 1
      endif

      !Update pointers
      Vx_xtemp => Vx_xold
      Vx_ztemp => Vx_zold
      Vz_xtemp => Vz_xold
      Vz_ztemp => Vz_zold
      sigmaxx_xtemp => sigmaxx_xold
      sigmaxx_ztemp => sigmaxx_zold
      sigmazz_xtemp => sigmazz_xold
      sigmazz_ztemp => sigmazz_zold
      sigmaxz_xtemp => sigmaxz_xold
      sigmaxz_ztemp => sigmaxz_zold

      Vx_xold => Vx_xcur
      Vx_zold => Vx_zcur
      Vz_xold => Vz_xcur
      Vz_zold => Vz_zcur
      sigmaxx_xold => sigmaxx_xcur
      sigmaxx_zold => sigmaxx_zcur
      sigmazz_xold => sigmazz_xcur
      sigmazz_zold => sigmazz_zcur
      sigmaxz_xold => sigmaxz_xcur
      sigmaxz_zold => sigmaxz_zcur

      Vx_xcur => Vx_xnew
      Vx_zcur => Vx_znew
      Vz_xcur => Vz_xnew
      Vz_zcur => Vz_znew
      sigmaxx_xcur => sigmaxx_xnew
      sigmaxx_zcur => sigmaxx_znew
      sigmazz_xcur => sigmazz_xnew
      sigmazz_zcur => sigmazz_znew
      sigmaxz_xcur => sigmaxz_xnew
      sigmaxz_zcur => sigmaxz_znew

      Vx_xnew => Vx_xtemp
      Vx_znew => Vx_ztemp
      Vz_xnew => Vz_xtemp
      Vz_znew => Vz_ztemp
      sigmaxx_xnew => sigmaxx_xtemp
      sigmaxx_znew => sigmaxx_ztemp
      sigmazz_xnew => sigmazz_xtemp
      sigmazz_znew => sigmazz_ztemp
      sigmaxz_xnew => sigmaxz_xtemp
      sigmaxz_znew => sigmaxz_ztemp

      Vx_xtemp_b => Vx_xold_b
      Vx_ztemp_b => Vx_zold_b
      Vz_xtemp_b => Vz_xold_b
      Vz_ztemp_b => Vz_zold_b
      sigmaxx_xtemp_b => sigmaxx_xold_b
      sigmaxx_ztemp_b => sigmaxx_zold_b
      sigmazz_xtemp_b => sigmazz_xold_b
      sigmazz_ztemp_b => sigmazz_zold_b
      sigmaxz_xtemp_b => sigmaxz_xold_b
      sigmaxz_ztemp_b => sigmaxz_zold_b

      Vx_xold_b => Vx_xcur_b
      Vx_zold_b => Vx_zcur_b
      Vz_xold_b => Vz_xcur_b
      Vz_zold_b => Vz_zcur_b
      sigmaxx_xold_b => sigmaxx_xcur_b
      sigmaxx_zold_b => sigmaxx_zcur_b
      sigmazz_xold_b => sigmazz_xcur_b
      sigmazz_zold_b => sigmazz_zcur_b
      sigmaxz_xold_b => sigmaxz_xcur_b
      sigmaxz_zold_b => sigmaxz_zcur_b

      Vx_xcur_b => Vx_xnew_b
      Vx_zcur_b => Vx_znew_b
      Vz_xcur_b => Vz_xnew_b
      Vz_zcur_b => Vz_znew_b
      sigmaxx_xcur_b => sigmaxx_xnew_b
      sigmaxx_zcur_b => sigmaxx_znew_b
      sigmazz_xcur_b => sigmazz_xnew_b
      sigmazz_zcur_b => sigmazz_znew_b
      sigmaxz_xcur_b => sigmaxz_xnew_b
      sigmaxz_zcur_b => sigmaxz_znew_b

      Vx_xnew_b => Vx_xtemp_b
      Vx_znew_b => Vx_ztemp_b
      Vz_xnew_b => Vz_xtemp_b
      Vz_znew_b => Vz_ztemp_b
      sigmaxx_xnew_b => sigmaxx_xtemp_b
      sigmaxx_znew_b => sigmaxx_ztemp_b
      sigmazz_xnew_b => sigmazz_xtemp_b
      sigmazz_znew_b => sigmazz_ztemp_b
      sigmaxz_xnew_b => sigmaxz_xtemp_b
      sigmaxz_znew_b => sigmaxz_ztemp_b

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif
    enddo
    write(0,'(a)',advance="yes")"]"
    write(0,*)"Last seismogram sample is:",countseis

  end subroutine

end module
