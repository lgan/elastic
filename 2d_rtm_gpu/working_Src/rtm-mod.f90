! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Jan, 22th, 2015

! Forward Born modeling program, based on the velocity-stress
! elastic modeling program

! FD elastic operator module
! This is the main subroutine that calls
! the FD derivatives, sources, etc.
module rtm_mod

  use rtmlib_mod
  use sep
  use source_mod

  implicit none
  integer, private :: nt, nx, nz, xsource, zsource
  real, private    :: dt, dx, dz, fp

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine elastic_init(nz_in,dz_in,nx_in,dx_in,dt_in,xsource_in,zsource_in,fp_in)

    integer             :: nz_in,nx_in, xsource_in, zsource_in
    real                :: dz_in,dx_in,dt_in, fp_in

    nz = nz_in
    dz = dz_in
    dt = dt_in
    nx = nx_in
    dx = dx_in
    xsource = xsource_in
    zsource = zsource_in
    fp = fp_in

  end subroutine

  ! Function that is called by Main
  function elastic_op(lambda0, mi0, rho0,dlambda,dmi,drho) result(stat)
    real,dimension(:,:)         :: lambda0, mi0, rho0
    real,dimension(:,:)         :: dlambda, dmi, drho
    integer                     :: stat

    stat=1
    call elastic_op2(lambda0,mi0,rho0,dlambda,dmi,drho)
    stat=0
  end function

  ! Subroutine used by function
  subroutine elastic_op2(lambda0, mi0, rho0, dlambda, dmi, drho)
    integer                                 :: it, ix, iz
    real, dimension(nz,nx)                  :: lambda0,mi0,rho0
    real, dimension(nz,nx)                  :: dlambda,dmi,drho
    real                                    :: aux1
    integer                                 :: stat

    !Background properties
    real, allocatable, dimension(:,:)     :: l2mi0, dl2mi
    real, allocatable, dimension(:,:)     :: rho0_Vx,rho0_Vz, drho_Vx, drho_Vz

    !model 1
    real, allocatable, dimension(:,:,:)   :: Vx_x,Vx_z,Vz_x,Vz_z
    real, allocatable, dimension(:,:,:)   :: sigmaxx_x,sigmaxx_z,sigmazz_x,sigmazz_z
    real, allocatable, dimension(:,:,:)   :: sigmaxz_x,sigmaxz_z

    !background data
    real, allocatable, dimension(:,:,:)   :: Vx_x0,Vx_z0,Vz_x0,Vz_z0
    real, allocatable, dimension(:,:,:)   :: sigmaxx_x0,sigmaxx_z0,sigmazz_x0,sigmazz_z0
    real, allocatable, dimension(:,:,:)   :: sigmaxz_x0,sigmaxz_z0

    !Scatterer data
    real, allocatable, dimension(:,:,:)   :: Vx_xB,Vx_zB,Vz_xB,Vz_zB
    real, allocatable, dimension(:,:,:)   :: sigmaxx_xB,sigmaxx_zB,sigmazz_xB,sigmazz_zB
    real, allocatable, dimension(:,:,:)   :: sigmaxz_xB,sigmaxz_zB

    !data 2
!    real, allocatable, dimension(:,:,:)   :: Vx_x2,Vx_z2,Vz_x2,Vz_z2
!    real, allocatable, dimension(:,:,:)   :: sigmaxx_x2,sigmaxx_z2,sigmazz_x2,sigmazz_z2
!    real, allocatable, dimension(:,:,:)   :: sigmaxz_x2,sigmaxz_z2

    !model 2
    real, allocatable, dimension(:,:,:)   :: Vx_xB_adj,Vx_zB_adj,Vz_xB_adj,Vz_zB_adj
    real, allocatable, dimension(:,:,:)   :: sigmaxx_xB_adj,sigmaxx_zB_adj,sigmazz_xB_adj,sigmazz_zB_adj
    real, allocatable, dimension(:,:,:)   :: sigmaxz_xB_adj,sigmaxz_zB_adj

    !model 2 scaled by Y'
    real, allocatable, dimension(:,:,:)   :: Vx_xB_adj2,Vx_zB_adj2,Vz_xB_adj2,Vz_zB_adj2
    real, allocatable, dimension(:,:,:)   :: sigmaxx_xB_adj2,sigmaxx_zB_adj2,sigmazz_xB_adj2,sigmazz_zB_adj2
    real, allocatable, dimension(:,:,:)   :: sigmaxz_xB_adj2,sigmaxz_zB_adj2
 
   !model 2 scaled by theta'
    real, allocatable, dimension(:,:,:)   :: Vx_xB_adj3,Vx_zB_adj3,Vz_xB_adj3,Vz_zB_adj3
    real, allocatable, dimension(:,:,:)   :: sigmaxx_xB_adj3,sigmaxx_zB_adj3,sigmazz_xB_adj3,sigmazz_zB_adj3
    real, allocatable, dimension(:,:,:)   :: sigmaxz_xB_adj3,sigmaxz_zB_adj3

    !Image
    real, allocatable, dimension(:,:)     ::r_rhox, r_rhoz, r_l2mix, r_l2miz, r_lx, r_lz, r_mi

    !Source terms
    real                                   :: source1

    !Output
    real, allocatable, dimension(:,:)        :: data0
    real, allocatable, dimension(:,:,:)      :: data1
    real, allocatable, dimension(:)          :: out1

    nt=1400

    !Allocate lame parameters
    allocate(l2mi0(nz,nx))
    allocate(rho0_Vx(nz,nx))
    allocate(rho0_Vz(nz,nx))
    allocate(dl2mi(nz,nx))
    allocate(drho_Vx(nz,nx))
    allocate(drho_Vz(nz,nx))

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
    allocate(Vx_x0(nz,nx,nt))
    allocate(Vx_z0(nz,nx,nt))
    allocate(Vz_x0(nz,nx,nt))
    allocate(Vz_z0(nz,nx,nt))
    allocate(sigmaxx_x0(nz,nx,nt))
    allocate(sigmaxx_z0(nz,nx,nt))
    allocate(sigmazz_x0(nz,nx,nt))
    allocate(sigmazz_z0(nz,nx,nt))
    allocate(sigmaxz_x0(nz,nx,nt))
    allocate(sigmaxz_z0(nz,nx,nt))

    !data 2
    allocate(Vx_xB(nz,nx,nt))
    allocate(Vx_zB(nz,nx,nt))
    allocate(Vz_xB(nz,nx,nt))
    allocate(Vz_zB(nz,nx,nt))
    allocate(sigmaxx_xB(nz,nx,nt))
    allocate(sigmaxx_zB(nz,nx,nt))
    allocate(sigmazz_xB(nz,nx,nt))
    allocate(sigmazz_zB(nz,nx,nt))
    allocate(sigmaxz_xB(nz,nx,nt))
    allocate(sigmaxz_zB(nz,nx,nt))

    !data 2 adj
    allocate(Vx_xB_adj(nz,nx,nt))
    allocate(Vx_zB_adj(nz,nx,nt))
    allocate(Vz_xB_adj(nz,nx,nt))
    allocate(Vz_zB_adj(nz,nx,nt))
    allocate(sigmaxx_xB_adj(nz,nx,nt))
    allocate(sigmaxx_zB_adj(nz,nx,nt))
    allocate(sigmazz_xB_adj(nz,nx,nt))
    allocate(sigmazz_zB_adj(nz,nx,nt))
    allocate(sigmaxz_xB_adj(nz,nx,nt))
    allocate(sigmaxz_zB_adj(nz,nx,nt))

    !data 2 adj
    allocate(Vx_xB_adj2(nz,nx,nt))
    allocate(Vx_zB_adj2(nz,nx,nt))
    allocate(Vz_xB_adj2(nz,nx,nt))
    allocate(Vz_zB_adj2(nz,nx,nt))
    allocate(sigmaxx_xB_adj2(nz,nx,nt))
    allocate(sigmaxx_zB_adj2(nz,nx,nt))
    allocate(sigmazz_xB_adj2(nz,nx,nt))
    allocate(sigmazz_zB_adj2(nz,nx,nt))
    allocate(sigmaxz_xB_adj2(nz,nx,nt))
    allocate(sigmaxz_zB_adj2(nz,nx,nt))

    !data 3 adj
    allocate(Vx_xB_adj3(nz,nx,nt))
    allocate(Vx_zB_adj3(nz,nx,nt))
    allocate(Vz_xB_adj3(nz,nx,nt))
    allocate(Vz_zB_adj3(nz,nx,nt))
    allocate(sigmaxx_xB_adj3(nz,nx,nt))
    allocate(sigmaxx_zB_adj3(nz,nx,nt))
    allocate(sigmazz_xB_adj3(nz,nx,nt))
    allocate(sigmazz_zB_adj3(nz,nx,nt))
    allocate(sigmaxz_xB_adj3(nz,nx,nt))
    allocate(sigmaxz_zB_adj3(nz,nx,nt))

    !Output
    allocate(r_rhox(nz,nx))
    allocate(r_rhoz(nz,nx))
    allocate(r_l2mix(nz,nx))
    allocate(r_l2miz(nz,nx))
    allocate(r_lx(nz,nx))
    allocate(r_lz(nz,nx))
    allocate(r_mi(nz,nx))
    allocate(data0(nz,nx))
    allocate(data1(nz,nx,nt))
    allocate(out1(nt))

    !Initialize values
    !Background properties
    l2mi0       =0.0
    rho0_Vx     =0.0
    rho0_Vz     =0.0
    dl2mi       =0.0
    drho_Vx     =0.0
    drho_Vz     =0.0

    !model 1
    Vx_x        =0.0
    Vx_z        =0.0
    Vz_x        =0.0
    Vz_z        =0.0
    sigmaxx_x   =0.0
    sigmaxx_z   =0.0
    sigmazz_x   =0.0
    sigmazz_z   =0.0
    sigmaxz_x   =0.0
    sigmaxz_z   =0.0

    !data 1
    Vx_x0       =0.0
    Vx_z0       =0.0
    Vz_x0       =0.0
    Vz_z0       =0.0
    sigmaxx_x0  =0.0
    sigmaxx_z0  =0.0
    sigmazz_x0  =0.0
    sigmazz_z0  =0.0
    sigmaxz_x0  =0.0
    sigmaxz_z0  =0.0

    !data 2
    Vx_xB       =0.0
    Vx_zB       =0.0
    Vz_xB       =0.0
    Vz_zB       =0.0
    sigmaxx_xB  =0.0
    sigmaxx_zB  =0.0
    sigmazz_xB  =0.0
    sigmazz_zB  =0.0
    sigmaxz_xB  =0.0
    sigmaxz_zB  =0.0

    !Output
    data0 = 0.
    data1 = 0.
    out1 = 0.

    !Calculate inverse of density
    !These also calculate the values of density at the velocity staggered grid,
    !using the arithmetic average proposed in Moczo, 2002

    !Rho_vx
    do ix=1,nx
      do iz=1,nz-1
        rho0_Vx(iz,ix)=1.0/(0.5*(rho0(iz+1,ix) + rho0(iz,ix)))
        if (drho(iz+1,ix).le.1.0e-1.or.drho(iz,ix).le.1.0e-1) then
          drho_Vx(iz,ix)=0.0
        else
          drho_Vx(iz,ix)=-(drho(iz+1,ix)/rho0(iz+1,ix) +&
                           drho(iz  ,ix)/rho0(iz  ,ix))
        endif
      enddo
    enddo
    !Rho_vz
    do ix=1,nx-1
      do iz=1,nz
        rho0_Vz(iz,ix)=1.0/(0.5*(rho0(iz,ix+1) + rho0(iz,ix)))
        if (drho(iz,ix+1).le.1.0e-1.or.drho(iz,ix).le.1.0e-1) then
          drho_Vz(iz,ix)=0.0
        else
          drho_Vz(iz,ix)=-(drho(iz,ix+1)/rho0(iz,ix+1) +&
                           drho(iz,ix  )/rho0(iz,ix  ))
        endif
      enddo
    enddo
    !Right side
    do iz=1,nz
      rho0_Vz(iz,nx)=(1.0/rho0(iz,nx))
      if (drho(iz,nx).le.1.0e-1) then
        drho_Vz(iz,nx)=0.0
      else
        drho_Vz(iz,nx)=-drho(iz,nx)/rho0(iz,nx)
      endif
    enddo
    !Bottom
    do ix=1,nx
      rho0_Vx(nz,ix)=(1.0/rho0(nz,ix))
      if (drho(nz,ix).le.1.0e-1) then
        drho_Vx(nz,ix)=0.0
      else
        drho_Vx(nz,ix)=-drho(nz,ix)/rho0(nz,ix)
      endif
    enddo
    
    !Calculate weigthed average for the lambda variable
    !These put lambda and (lambda+2mi) in their respective staggered grids,
    !using the harmonic average proposed in Moczo, 2002
    do ix=2,nx-1
      do iz=2,nz-1
        lambda0(iz,ix)=4.0/((1.0/lambda0(iz,ix))+(1.0/lambda0(iz+1,ix))+(1.0/lambda0(iz,ix+1))+(1.0/lambda0(iz+1,ix+1)))
        if (mi0(iz,ix)==0.or.mi0(iz+1,ix)==0.or.mi0(iz,ix+1)==0.or.mi0(iz+1,ix+1)==0) then
          mi0(iz,ix)=0.0
        else
          mi0(iz,ix)=4.0/((1.0/mi0(iz,ix))+(1.0/mi0(iz+1,ix))+(1.0/mi0(iz,ix+1))+(1.0/mi0(iz+1,ix+1)))
        endif
        if (dlambda(iz,ix)<1e-2.or.dlambda(iz+1,ix)<1e-2.or.dlambda(iz,ix+1)<1e-2.or.dlambda(iz+1,ix+1)<1e-2) then
          dlambda(iz,ix)=0.0
        else
          dlambda(iz,ix)=4.0/((1.0/dlambda(iz,ix))+(1.0/dlambda(iz+1,ix))+(1.0/dlambda(iz,ix+1))+(1.0/dlambda(iz+1,ix+1)))
        endif
        if (dmi(iz,ix)<1e-2.or.dmi(iz+1,ix)<1e-2.or.dmi(iz,ix+1)<1e-2.or.dmi(iz+1,ix+1)<1e-2) then
          dmi(iz,ix)=0.0
        else
          dmi(iz,ix)=4.0/((1.0/dmi(iz,ix))+(1.0/dmi(iz+1,ix))+(1.0/dmi(iz,ix+1))+(1.0/dmi(iz+1,ix+1)))
        endif
      enddo
    enddo

    l2mi0 = lambda0 + 2*mi0

    dl2mi=(dlambda+2.0*dmi)/l2mi0
    dlambda=dlambda/l2mi0
    dmi=dmi/mi0

    aux1=(dt/dx)

    !Applying dt/dx to background properties
    rho0_Vx = aux1*rho0_Vx   
    rho0_Vz = aux1*rho0_Vz   
    lambda0 = aux1*lambda0
    mi0=aux1*mi0
    l2mi0=aux1*l2mi0

    !Initialize subroutines
    call rtmmod_init(nz,nx,dx,rho0_Vx,rho0_Vz,lambda0,mi0,l2mi0,drho_Vx,drho_Vz,dlambda,dmi,dl2mi)
    call source_init(fp,dt,dx,250,nz,nx)

!###################################################
! Create source for model 1
!###################################################
    do it=1,nt
      call source_op2(it,source1)
      call explosive_sourcex(source1,xsource,zsource,sigmaxx_x(:,:,it))
      call explosive_sourcez(source1,xsource,zsource,sigmazz_z(:,:,it))
    enddo

    out1(:) = sigmazz_z(zsource,xsource,:)

    call to_history("n1", nt, "Dat/source1.H")
    call to_history("o1", 0., "Dat/source1.H")
    call to_history("d1", dt, "Dat/source1.H")
         
    call sep_write(out1,"Dat/source1.H")

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

! First time steps
    !it=1
    Vx_x0(:,:,1) = Vx_x(:,:,1)
    Vx_z0(:,:,1) = Vx_z(:,:,1)
    Vz_x0(:,:,1) = Vz_x(:,:,1)
    Vz_z0(:,:,1) = Vz_z(:,:,1)
    sigmaxx_x0(:,:,1) = sigmaxx_x(:,:,1)
    sigmaxx_z0(:,:,1) = sigmaxx_z(:,:,1)
    sigmazz_x0(:,:,1) = sigmazz_x(:,:,1)
    sigmazz_z0(:,:,1) = sigmazz_z(:,:,1)
    sigmaxz_x0(:,:,1) = sigmaxz_x(:,:,1)
    sigmaxz_z0(:,:,1) = sigmaxz_z(:,:,1)

    !it=2
    !Vx_x
    Vx_x0(:,:,2) = Vx_x(:,:,2)
    call Ax(.false.,.true.,sigmaxx_x0(:,:,1),Vx_x0(:,:,2))
    call Ax(.false.,.true.,sigmaxx_z0(:,:,1),Vx_x0(:,:,2))
    !Vx_z
    Vx_z0(:,:,2) = Vx_z(:,:,2)
    call Bx(.false.,.true.,sigmaxz_x0(:,:,1),Vx_z0(:,:,2))
    call Bx(.false.,.true.,sigmaxz_z0(:,:,1),Vx_z0(:,:,2))
    !Vz_x
    Vz_x0(:,:,2) = Vz_x(:,:,2)
    call Az(.false.,.true.,sigmaxz_x0(:,:,1),Vz_x0(:,:,2))
    call Az(.false.,.true.,sigmaxz_z0(:,:,1),Vz_x0(:,:,2))
    !Vz_z
    Vz_z0(:,:,2) = Vz_z(:,:,2)
    call Bz(.false.,.true.,sigmazz_x0(:,:,1),Vz_z0(:,:,2))
    call Bz(.false.,.true.,sigmazz_z0(:,:,1),Vz_z0(:,:,2))
    !sigmaxx_x
    sigmaxx_x0(:,:,2) = sigmaxx_x(:,:,2)
    call D(.false.,.true.,Vx_x0(:,:,1),sigmaxx_x0(:,:,2))
    call D(.false.,.true.,Vx_z0(:,:,1),sigmaxx_x0(:,:,2))
    !sigmaxx_z
    sigmaxx_z0(:,:,2) = sigmaxx_z(:,:,2)
    call E(.false.,.true.,Vz_x0(:,:,1),sigmaxx_z0(:,:,2))
    call E(.false.,.true.,Vz_z0(:,:,1),sigmaxx_z0(:,:,2))
    !sigmazz_x
    sigmazz_x0(:,:,2) = sigmazz_x(:,:,2)
    call F(.false.,.true.,Vx_x0(:,:,1),sigmazz_x0(:,:,2))
    call F(.false.,.true.,Vx_z0(:,:,1),sigmazz_x0(:,:,2))
    !sigmazz_z
    sigmazz_z0(:,:,2) = sigmazz_z(:,:,2)
    call G(.false.,.true.,Vz_x0(:,:,1),sigmazz_z0(:,:,2))
    call G(.false.,.true.,Vz_z0(:,:,1),sigmazz_z0(:,:,2))
    !sigmaxz_x
    sigmaxz_x0(:,:,2) = sigmaxz_x(:,:,2)
    call H(.false.,.true.,Vx_x0(:,:,1),sigmaxz_x0(:,:,2))
    call H(.false.,.true.,Vx_z0(:,:,1),sigmaxz_x0(:,:,2))
    !sigmaxz_z
    sigmaxz_z0(:,:,2) = sigmaxz_z(:,:,2)
    call J(.false.,.true.,Vz_x0(:,:,1),sigmaxz_z0(:,:,2))
    call J(.false.,.true.,Vz_z0(:,:,1),sigmaxz_z0(:,:,2))

    !it=3 until it=nt
    do it=3,nt
      !Vx_x
      Vx_x0(:,:,it) = Vx_x(:,:,it) + Vx_x0(:,:,(it-2))
      call Ax(.false.,.true.,sigmaxx_x0(:,:,(it-1)),Vx_x0(:,:,it))
      call Ax(.false.,.true.,sigmaxx_z0(:,:,(it-1)),Vx_x0(:,:,it))
      !Vx_z
      Vx_z0(:,:,it) = Vx_z(:,:,it) + Vx_z0(:,:,(it-2))
      call Bx(.false.,.true.,sigmaxz_x0(:,:,(it-1)),Vx_z0(:,:,it))
      call Bx(.false.,.true.,sigmaxz_z0(:,:,(it-1)),Vx_z0(:,:,it))
      !Vz_x
      Vz_x0(:,:,it) = Vz_x(:,:,it) + Vz_x0(:,:,(it-2))
      call Az(.false.,.true.,sigmaxz_x0(:,:,(it-1)),Vz_x0(:,:,it))
      call Az(.false.,.true.,sigmaxz_z0(:,:,(it-1)),Vz_x0(:,:,it))
      !Vz_z
      Vz_z0(:,:,it) = Vz_z(:,:,it) + Vz_z0(:,:,(it-2))
      call Bz(.false.,.true.,sigmazz_x0(:,:,(it-1)),Vz_z0(:,:,it))
      call Bz(.false.,.true.,sigmazz_z0(:,:,(it-1)),Vz_z0(:,:,it))
      !sigmaxx_x
      sigmaxx_x0(:,:,it) = sigmaxx_x(:,:,it) + sigmaxx_x0(:,:,(it-2))
      call D(.false.,.true.,Vx_x0(:,:,(it-1)),sigmaxx_x0(:,:,it))
      call D(.false.,.true.,Vx_z0(:,:,(it-1)),sigmaxx_x0(:,:,it))
      !sigmaxx_z
      sigmaxx_z0(:,:,it) = sigmaxx_z(:,:,it) + sigmaxx_z0(:,:,(it-2))
      call E(.false.,.true.,Vz_x0(:,:,(it-1)),sigmaxx_z0(:,:,it))
      call E(.false.,.true.,Vz_z0(:,:,(it-1)),sigmaxx_z0(:,:,it))
      !sigmazz_x
      sigmazz_x0(:,:,it) = sigmazz_x(:,:,it) + sigmazz_x0(:,:,(it-2))
      call F(.false.,.true.,Vx_x0(:,:,(it-1)),sigmazz_x0(:,:,it))
      call F(.false.,.true.,Vx_z0(:,:,(it-1)),sigmazz_x0(:,:,it))
      !sigmazz_z
      sigmazz_z0(:,:,it) = sigmazz_z(:,:,it) + sigmazz_z0(:,:,(it-2))
      call G(.false.,.true.,Vz_x0(:,:,(it-1)),sigmazz_z0(:,:,it))
      call G(.false.,.true.,Vz_z0(:,:,(it-1)),sigmazz_z0(:,:,it))
      !sigmaxz_x
      sigmaxz_x0(:,:,it) = sigmaxz_x(:,:,it) + sigmaxz_x0(:,:,(it-2))
      call H(.false.,.true.,Vx_x0(:,:,(it-1)),sigmaxz_x0(:,:,it))
      call H(.false.,.true.,Vx_z0(:,:,(it-1)),sigmaxz_x0(:,:,it))
      !sigmaxz_z
      sigmaxz_z0(:,:,it) = sigmaxz_z(:,:,it) + sigmaxz_z0(:,:,(it-2))
      call J(.false.,.true.,Vz_x0(:,:,(it-1)),sigmaxz_z0(:,:,it))
      call J(.false.,.true.,Vz_z0(:,:,(it-1)),sigmaxz_z0(:,:,it))

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif

    enddo
    write(0,'(a)',advance="yes")"]"

    data1 = 0.
    data1 = sigmaxx_x0 + sigmaxx_z0 + sigmazz_x0 + sigmazz_z0

    call to_history("n1", nz, "Dat/data0.H")
    call to_history("o1", 0., "Dat/data0.H")
    call to_history("d1", dz, "Dat/data0.H")
    call to_history("n2", nx, "Dat/data0.H")
    call to_history("o2", 0., "Dat/data0.H")
    call to_history("d2", dx, "Dat/data0.H")
    call to_history("n3", nt, "Dat/data0.H")
    call to_history("o3", 0., "Dat/data0.H")
    call to_history("d3", dt, "Dat/data0.H")
         
    call sep_write(data1,"Dat/data0.H")

!###################################################
! Forward scattered data delta_d
!###################################################
    write(0,*)"============================================================="
    write(0,*)"BEGIN OF TIME MODELING FOR SCATTERED WAVEFIELD"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

! First time steps
    !it=1
    Vx_xB(:,:,1) = 0.
    Vx_zB(:,:,1) = 0.
    Vz_xB(:,:,1) = 0.
    Vz_zB(:,:,1) = 0.
    sigmaxx_xB(:,:,1) = 0.
    sigmaxx_zB(:,:,1) = 0.
    sigmazz_xB(:,:,1) = 0.
    sigmazz_zB(:,:,1) = 0.
    sigmaxz_xB(:,:,1) = 0.
    sigmaxz_zB(:,:,1) = 0.

    !it=2
    !Vx_x
    call Ax(.false.,.false.,sigmaxx_xB(:,:,1),Vx_xB(:,:,2))
    call Ax(.false.,.true.,sigmaxx_zB(:,:,1),Vx_xB(:,:,2))
    !Vx_z
    call Bx(.false.,.false.,sigmaxz_xB(:,:,1),Vx_zB(:,:,2))
    call Bx(.false.,.true.,sigmaxz_zB(:,:,1),Vx_zB(:,:,2))
    !Vz_x
    call Az(.false.,.false.,sigmaxz_xB(:,:,1),Vz_xB(:,:,2))
    call Az(.false.,.true.,sigmaxz_zB(:,:,1),Vz_xB(:,:,2))
    !Vz_z
    call Bz(.false.,.false.,sigmazz_xB(:,:,1),Vz_zB(:,:,2))
    call Bz(.false.,.true.,sigmazz_zB(:,:,1),Vz_zB(:,:,2))
    !sigmaxx_x
    call S3(.false.,.false.,sigmaxx_x0(:,:,1:2),sigmaxx_xB(:,:,2))
    call S4(.false.,.true.,sigmazz_x0(:,:,1:2),sigmaxx_xB(:,:,2))
    call D(.false.,.true.,Vx_xB(:,:,1),sigmaxx_xB(:,:,2))
    call D(.false.,.true.,Vx_zB(:,:,1),sigmaxx_xB(:,:,2))
    !sigmaxx_z
    call S3(.false.,.false.,sigmaxx_z0(:,:,1:2),sigmaxx_zB(:,:,2))
    call S4(.false.,.true.,sigmazz_z0(:,:,1:2),sigmaxx_zB(:,:,2))
    call E(.false.,.true.,Vz_xB(:,:,1),sigmaxx_zB(:,:,2))
    call E(.false.,.true.,Vz_zB(:,:,1),sigmaxx_zB(:,:,2))
    !sigmazz_x
    call S3(.false.,.false.,sigmazz_x0(:,:,1:2),sigmazz_xB(:,:,2))
    call S4(.false.,.true.,sigmaxx_x0(:,:,1:2),sigmazz_xB(:,:,2))
    call F(.false.,.true.,Vx_xB(:,:,1),sigmazz_xB(:,:,2))
    call F(.false.,.true.,Vx_zB(:,:,1),sigmazz_xB(:,:,2))
    !sigmazz_z
    call S3(.false.,.false.,sigmazz_z0(:,:,1:2),sigmazz_zB(:,:,2))
    call S4(.false.,.true.,sigmaxx_z0(:,:,1:2),sigmazz_zB(:,:,2))
    call G(.false.,.true.,Vz_xB(:,:,1),sigmazz_zB(:,:,2))
    call G(.false.,.true.,Vz_zB(:,:,1),sigmazz_zB(:,:,2))
    !sigmaxz_x
    call S5(.false.,.false.,sigmaxz_x0(:,:,1:2),sigmaxz_xB(:,:,2))
    call H(.false.,.true.,Vx_xB(:,:,1),sigmaxz_xB(:,:,2))
    call H(.false.,.true.,Vx_zB(:,:,1),sigmaxz_xB(:,:,2))
    !sigmaxz_z
    call S5(.false.,.false.,sigmaxz_z0(:,:,1:2),sigmaxz_zB(:,:,2))
    call J(.false.,.true.,Vz_xB(:,:,1),sigmaxz_zB(:,:,2))
    call J(.false.,.true.,Vz_zB(:,:,1),sigmaxz_zB(:,:,2))

    !it=3 until it=nt
    do it=3,nt
      !Vx_x
      Vx_xB(:,:,it) = Vx_xB(:,:,(it-2))
      call S1(.false.,.true.,Vx_x0(:,:,(it-2):(it-1)),Vx_xB(:,:,it))
      call Ax(.false.,.true.,sigmaxx_xB(:,:,(it-1)),Vx_xB(:,:,it))
      call Ax(.false.,.true.,sigmaxx_zB(:,:,(it-1)),Vx_xB(:,:,it))
      !Vx_z
      Vx_zB(:,:,it) = Vx_zB(:,:,(it-2))
      call S1(.false.,.true.,Vx_z0(:,:,(it-2):(it-1)),Vx_zB(:,:,it))
      call Bx(.false.,.true.,sigmaxz_xB(:,:,(it-1)),Vx_zB(:,:,it))
      call Bx(.false.,.true.,sigmaxz_zB(:,:,(it-1)),Vx_zB(:,:,it))
      !Vz_x
      Vz_xB(:,:,it) = Vz_xB(:,:,(it-2))
      call S2(.false.,.true.,Vz_x0(:,:,(it-2):(it-1)),Vz_xB(:,:,it))
      call Az(.false.,.true.,sigmaxz_xB(:,:,(it-1)),Vz_xB(:,:,it))
      call Az(.false.,.true.,sigmaxz_zB(:,:,(it-1)),Vz_xB(:,:,it))
      !Vz_z
      Vz_zB(:,:,it) = Vz_zB(:,:,(it-2))
      call S2(.false.,.true.,Vz_z0(:,:,(it-2):(it-1)),Vz_zB(:,:,it))
      call Bz(.false.,.true.,sigmazz_xB(:,:,(it-1)),Vz_zB(:,:,it))
      call Bz(.false.,.true.,sigmazz_zB(:,:,(it-1)),Vz_zB(:,:,it))
      !sigmaxx_x
      sigmaxx_xB(:,:,it) = sigmaxx_xB(:,:,(it-2))
      call S3(.false.,.true.,sigmaxx_x0(:,:,(it-1):it),sigmaxx_xB(:,:,it))
      call S4(.false.,.true.,sigmazz_x0(:,:,(it-1):it),sigmaxx_xB(:,:,it))
      call D(.false.,.true.,Vx_xB(:,:,(it-1)),sigmaxx_xB(:,:,it))
      call D(.false.,.true.,Vx_zB(:,:,(it-1)),sigmaxx_xB(:,:,it))
      !sigmaxx_z
      sigmaxx_zB(:,:,it) = sigmaxx_zB(:,:,(it-2))
      call S3(.false.,.true.,sigmaxx_z0(:,:,(it-1):it),sigmaxx_zB(:,:,it))
      call S4(.false.,.true.,sigmazz_z0(:,:,(it-1):it),sigmaxx_zB(:,:,it))
      call E(.false.,.true.,Vz_xB(:,:,(it-1)),sigmaxx_zB(:,:,it))
      call E(.false.,.true.,Vz_zB(:,:,(it-1)),sigmaxx_zB(:,:,it))
      !sigmazz_x
      sigmazz_xB(:,:,it) = sigmazz_xB(:,:,(it-2))
      call S3(.false.,.true.,sigmazz_x0(:,:,(it-1):it),sigmazz_xB(:,:,it))
      call S4(.false.,.true.,sigmaxx_x0(:,:,(it-1):it),sigmazz_xB(:,:,it))
      call F(.false.,.true.,Vx_xB(:,:,(it-1)),sigmazz_xB(:,:,it))
      call F(.false.,.true.,Vx_zB(:,:,(it-1)),sigmazz_xB(:,:,it))
      !sigmazz_z
      sigmazz_zB(:,:,it) = sigmazz_zB(:,:,(it-2))
      call S3(.false.,.true.,sigmazz_z0(:,:,(it-1):it),sigmazz_zB(:,:,it))
      call S4(.false.,.true.,sigmaxx_z0(:,:,(it-1):it),sigmazz_zB(:,:,it))
      call G(.false.,.true.,Vz_xB(:,:,(it-1)),sigmazz_zB(:,:,it))
      call G(.false.,.true.,Vz_zB(:,:,(it-1)),sigmazz_zB(:,:,it))
      !sigmaxz_x
      sigmaxz_xB(:,:,it) = sigmaxz_xB(:,:,(it-2))
      call S5(.false.,.true.,sigmaxz_x0(:,:,(it-1):it),sigmaxz_xB(:,:,it))
      call H(.false.,.true.,Vx_xB(:,:,(it-1)),sigmaxz_xB(:,:,it))
      call H(.false.,.true.,Vx_zB(:,:,(it-1)),sigmaxz_xB(:,:,it))
      !sigmaxz_z
      sigmaxz_zB(:,:,it) = sigmaxz_zB(:,:,(it-2))
      call S5(.false.,.true.,sigmaxz_z0(:,:,(it-1):it),sigmaxz_zB(:,:,it))
      call J(.false.,.true.,Vz_xB(:,:,(it-1)),sigmaxz_zB(:,:,it))
      call J(.false.,.true.,Vz_zB(:,:,(it-1)),sigmaxz_zB(:,:,it))

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif

    enddo
    write(0,'(a)',advance="yes")"]"

    data1 = 0.
    data1 = sigmaxx_xB + sigmaxx_zB + sigmazz_xB + sigmazz_zB

    call to_history("n1", nz, "Dat/dataB.H")
    call to_history("o1", 0., "Dat/dataB.H")
    call to_history("d1", dz, "Dat/dataB.H")
    call to_history("n2", nx, "Dat/dataB.H")
    call to_history("o2", 0., "Dat/dataB.H")
    call to_history("d2", dx, "Dat/dataB.H")
    call to_history("n3", nt, "Dat/dataB.H")
    call to_history("o3", 0., "Dat/dataB.H")
    call to_history("d3", dt, "Dat/dataB.H")
         
    call sep_write(data1,"Dat/dataB.H")

!###################################################
! Adjoint application to delta d
!###################################################
    write(0,*)"============================================================="
    write(0,*)"BEGIN OF ADJOINT TIME MODELING FOR SCATTERED WAVEFIELD"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

    !data 2 adj
    Vx_xB_adj       =0.0
    Vx_zB_adj       =0.0
    Vz_xB_adj       =0.0
    Vz_zB_adj       =0.0
    sigmaxx_xB_adj  =0.0
    sigmaxx_zB_adj  =0.0
    sigmazz_xB_adj  =0.0
    sigmazz_zB_adj  =0.0
    sigmaxz_xB_adj  =0.0
    sigmaxz_zB_adj  =0.0

! Last time steps
    !it=nt
    Vx_xB_adj(:,:,nt) = Vx_xB(:,:,nt)
    Vx_zB_adj(:,:,nt) = Vx_zB(:,:,nt)
    Vz_xB_adj(:,:,nt) = Vz_xB(:,:,nt)
    Vz_zB_adj(:,:,nt) = Vz_zB(:,:,nt)
    sigmaxx_xB_adj(:,:,nt) = sigmaxx_xB(:,:,nt)
    sigmaxx_zB_adj(:,:,nt) = sigmaxx_zB(:,:,nt)
    sigmazz_xB_adj(:,:,nt) = sigmazz_xB(:,:,nt)
    sigmazz_zB_adj(:,:,nt) = sigmazz_zB(:,:,nt)
    sigmaxz_xB_adj(:,:,nt) = sigmaxz_xB(:,:,nt)
    sigmaxz_zB_adj(:,:,nt) = sigmaxz_zB(:,:,nt)

    !it=nt-1
    !Vx_x
    Vx_xB_adj(:,:,nt-1) = Vx_xB(:,:,nt-1)
    call D(.true.,.true.,Vx_xB_adj(:,:,(nt-1)),sigmaxx_xB_adj(:,:,nt))    
    call F(.true.,.true.,Vx_xB_adj(:,:,(nt-1)),sigmazz_xB_adj(:,:,nt))
    call H(.true.,.true.,Vx_xB_adj(:,:,(nt-1)),sigmaxz_xB_adj(:,:,nt))
    !Vx_z
    Vx_zB_adj(:,:,nt-1) = Vx_zB(:,:,nt-1)
    call D(.true.,.true.,Vx_zB_adj(:,:,(nt-1)),sigmaxx_xB_adj(:,:,nt))    
    call F(.true.,.true.,Vx_zB_adj(:,:,(nt-1)),sigmazz_xB_adj(:,:,nt))
    call H(.true.,.true.,Vx_zB_adj(:,:,(nt-1)),sigmaxz_xB_adj(:,:,nt))
    !Vz_x
    Vz_xB_adj(:,:,nt-1) = Vz_xB(:,:,nt-1)
    call E(.true.,.true.,Vz_xB_adj(:,:,(nt-1)),sigmaxx_zB_adj(:,:,nt))    
    call G(.true.,.true.,Vz_xB_adj(:,:,(nt-1)),sigmazz_zB_adj(:,:,nt))
    call J(.true.,.true.,Vz_xB_adj(:,:,(nt-1)),sigmaxz_zB_adj(:,:,nt))
    !Vz_z
    Vz_zB_adj(:,:,nt-1) = Vz_zB(:,:,nt-1)
    call E(.true.,.true.,Vz_zB_adj(:,:,(nt-1)),sigmaxx_zB_adj(:,:,nt))    
    call G(.true.,.true.,Vz_zB_adj(:,:,(nt-1)),sigmazz_zB_adj(:,:,nt))
    call J(.true.,.true.,Vz_zB_adj(:,:,(nt-1)),sigmaxz_zB_adj(:,:,nt))
    !sigmaxx_x
    sigmaxx_xB_adj(:,:,nt-1) = sigmaxx_xB(:,:,nt-1)
    call Ax(.true.,.true.,sigmaxx_xB_adj(:,:,(nt-1)),Vx_xB_adj(:,:,nt))
    !sigmaxx_z
    sigmaxx_zB_adj(:,:,nt-1) = sigmaxx_zB(:,:,nt-1)
    call Ax(.true.,.true.,sigmaxx_zB_adj(:,:,(nt-1)),Vx_xB_adj(:,:,nt))
    !sigmazz_x
    sigmazz_xB_adj(:,:,nt-1) = sigmazz_xB(:,:,nt-1)
    call Bz(.true.,.true.,sigmazz_xB_adj(:,:,(nt-1)),Vz_zB_adj(:,:,nt))
    !sigmazz_z
    sigmazz_zB_adj(:,:,nt-1) = sigmazz_zB(:,:,nt-1)
    call Bz(.true.,.true.,sigmazz_zB_adj(:,:,(nt-1)),Vz_zB_adj(:,:,nt))
    !sigmaxz_x
    sigmaxz_xB_adj(:,:,nt-1) = sigmaxz_xB(:,:,nt-1)
    call Bx(.true.,.true.,sigmaxz_xB_adj(:,:,(nt-1)),Vx_zB_adj(:,:,nt))
    call Az(.true.,.true.,sigmaxz_xB_adj(:,:,(nt-1)),Vz_xB_adj(:,:,nt))
    !sigmaxz_z
    sigmaxz_zB_adj(:,:,nt-1) = sigmaxz_zB(:,:,nt-1)
    call Bx(.true.,.true.,sigmaxz_zB_adj(:,:,(nt-1)),Vx_zB_adj(:,:,nt))
    call Az(.true.,.true.,sigmaxz_zB_adj(:,:,(nt-1)),Vz_xB_adj(:,:,nt))

    !it=nt-2 until it=1
    do it=nt-2,1,-1
      !Vx_x
      Vx_xB_adj(:,:,it) = Vx_xB(:,:,it) + Vx_xB_adj(:,:,(it+2))
      call D(.true.,.true.,Vx_xB_adj(:,:,it),sigmaxx_xB_adj(:,:,(it+1)))    
      call F(.true.,.true.,Vx_xB_adj(:,:,it),sigmazz_xB_adj(:,:,(it+1)))
      call H(.true.,.true.,Vx_xB_adj(:,:,it),sigmaxz_xB_adj(:,:,(it+1)))
      !Vx_z
      Vx_zB_adj(:,:,it) = Vx_zB(:,:,it) + Vx_zB_adj(:,:,(it+2))
      call D(.true.,.true.,Vx_zB_adj(:,:,it),sigmaxx_xB_adj(:,:,(it+1)))    
      call F(.true.,.true.,Vx_zB_adj(:,:,it),sigmazz_xB_adj(:,:,(it+1)))
      call H(.true.,.true.,Vx_zB_adj(:,:,it),sigmaxz_xB_adj(:,:,(it+1)))
      !Vz_x
      Vz_xB_adj(:,:,it) = Vz_xB(:,:,it) + Vz_xB_adj(:,:,(it+2))
      call E(.true.,.true.,Vz_xB_adj(:,:,it),sigmaxx_zB_adj(:,:,(it+1)))    
      call G(.true.,.true.,Vz_xB_adj(:,:,it),sigmazz_zB_adj(:,:,(it+1)))
      call J(.true.,.true.,Vz_xB_adj(:,:,it),sigmaxz_zB_adj(:,:,(it+1)))
      !Vz_z
      Vz_zB_adj(:,:,it) = Vz_zB(:,:,it) + Vz_zB_adj(:,:,(it+2))
      call E(.true.,.true.,Vz_zB_adj(:,:,it),sigmaxx_zB_adj(:,:,(it+1)))    
      call G(.true.,.true.,Vz_zB_adj(:,:,it),sigmazz_zB_adj(:,:,(it+1)))
      call J(.true.,.true.,Vz_zB_adj(:,:,it),sigmaxz_zB_adj(:,:,(it+1)))
      !sigmaxx_x
      sigmaxx_xB_adj(:,:,it) = sigmaxx_xB(:,:,it) + sigmaxx_xB_adj(:,:,(it+2))
      call Ax(.true.,.true.,sigmaxx_xB_adj(:,:,it),Vx_xB_adj(:,:,(it+1)))
      !sigmaxx_z
      sigmaxx_zB_adj(:,:,it) = sigmaxx_zB(:,:,it) + sigmaxx_zB_adj(:,:,(it+2))
      call Ax(.true.,.true.,sigmaxx_zB_adj(:,:,it),Vx_xB_adj(:,:,(it+1)))
      !sigmazz_x
      sigmazz_xB_adj(:,:,it) = sigmazz_xB(:,:,it) + sigmazz_xB_adj(:,:,(it+2))
      call Bz(.true.,.true.,sigmazz_xB_adj(:,:,it),Vz_zB_adj(:,:,(it+1)))
      !sigmazz_z
      sigmazz_zB_adj(:,:,it) = sigmazz_zB(:,:,it) + sigmazz_zB_adj(:,:,(it+2))
      call Bz(.true.,.true.,sigmazz_zB_adj(:,:,it),Vz_zB_adj(:,:,(it+1)))
      !sigmaxz_x
      sigmaxz_xB_adj(:,:,it) = sigmaxz_xB(:,:,it) + sigmaxz_xB_adj(:,:,(it+2))
      call Bx(.true.,.true.,sigmaxz_xB_adj(:,:,it),Vx_zB_adj(:,:,(it+1)))
      call Az(.true.,.true.,sigmaxz_xB_adj(:,:,it),Vz_xB_adj(:,:,(it+1)))
      !sigmaxz_z
      sigmaxz_zB_adj(:,:,it) = sigmaxz_zB(:,:,it) + sigmaxz_zB_adj(:,:,(it+2))
      call Bx(.true.,.true.,sigmaxz_zB_adj(:,:,it),Vx_zB_adj(:,:,(it+1)))
      call Az(.true.,.true.,sigmaxz_zB_adj(:,:,it),Vz_xB_adj(:,:,(it+1)))

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif

    enddo

    write(0,'(a)',advance="yes")"]"

!###################################################
! Scaling of adjoint
!###################################################
    write(0,*)"============================================================="
    write(0,*)"BEGIN OF SCALING OF ADJOINT BY Y'"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

    !data 2 adj
    Vx_xB_adj2       =0.0
    Vx_zB_adj2       =0.0
    Vz_xB_adj2       =0.0
    Vz_zB_adj2       =0.0
    sigmaxx_xB_adj2  =0.0
    sigmaxx_zB_adj2  =0.0
    sigmazz_xB_adj2  =0.0
    sigmazz_zB_adj2  =0.0
    sigmaxz_xB_adj2  =0.0
    sigmaxz_zB_adj2  =0.0

! Last time steps
    !it=nt
    Vx_xB_adj2(:,:,nt) = 0.
    Vx_zB_adj2(:,:,nt) = 0.
    Vz_xB_adj2(:,:,nt) = 0.
    Vz_zB_adj2(:,:,nt) = 0.
    sigmaxx_xB_adj2(:,:,nt) = sigmaxx_xB_adj(:,:,nt)
    sigmaxx_zB_adj2(:,:,nt) = sigmaxx_zB_adj(:,:,nt)
    sigmazz_xB_adj2(:,:,nt) = sigmazz_xB_adj(:,:,nt)
    sigmazz_zB_adj2(:,:,nt) = sigmazz_zB_adj(:,:,nt)
    sigmaxz_xB_adj2(:,:,nt) = sigmaxz_xB_adj(:,:,nt)
    sigmaxz_zB_adj2(:,:,nt) = sigmaxz_zB_adj(:,:,nt)
    
    !it=nt-1 until it=1
    do it=nt-1,1,-1
      !Vx_x
      call Cx(.true.,.false., Vx_xB_adj2(:,:,it),Vx_xB_adj(:,:,it+1))
      !Vx_z
      call Cx(.true.,.false., Vx_zB_adj2(:,:,it),Vx_zB_adj(:,:,it+1))
      !Vz_x
      call Cz(.true.,.false., Vz_xB_adj2(:,:,it),Vz_xB_adj(:,:,it+1))
      !Vz_z
      call Cz(.true.,.false., Vz_zB_adj2(:,:,it),Vz_zB_adj(:,:,it+1))
      !sigmaxx_x
      sigmaxx_xB_adj2(:,:,it) = sigmaxx_xB_adj(:,:,it)
      !sigmaxx_z
      sigmaxx_zB_adj2(:,:,it) = sigmaxx_zB_adj(:,:,it)
      !sigmazz_x
      sigmazz_xB_adj2(:,:,it) = sigmazz_xB_adj(:,:,it)
      !sigmazz_z
      sigmazz_zB_adj2(:,:,it) = sigmazz_zB_adj(:,:,it)
      !sigmaxz_x
      sigmaxz_xB_adj2(:,:,it) = sigmaxz_xB_adj(:,:,it)
      !sigmaxz_z
      sigmaxz_zB_adj2(:,:,it) = sigmaxz_zB_adj(:,:,it)

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif

    enddo

    write(0,'(a)',advance="yes")"]"

    data1 = 0.
    data1 = Vx_xB_adj2 + Vx_zB_adj2

    call to_history("n1", nz, "Dat/data-vx.H")
    call to_history("o1", 0., "Dat/data-vx.H")
    call to_history("d1", dz, "Dat/data-vx.H")
    call to_history("n2", nx, "Dat/data-vx.H")
    call to_history("o2", 0., "Dat/data-vx.H")
    call to_history("d2", dx, "Dat/data-vx.H")
    call to_history("n3", nt, "Dat/data-vx.H")
    call to_history("o3", 0., "Dat/data-vx.H")
    call to_history("d3", dt, "Dat/data-vx.H")
         
    call sep_write(data1,"Dat/data-vx.H")

    data1 = 0.
    data1 = Vz_xB_adj2 + Vz_zB_adj2

    call to_history("n1", nz, "Dat/data-vz.H")
    call to_history("o1", 0., "Dat/data-vz.H")
    call to_history("d1", dz, "Dat/data-vz.H")
    call to_history("n2", nx, "Dat/data-vz.H")
    call to_history("o2", 0., "Dat/data-vz.H")
    call to_history("d2", dx, "Dat/data-vz.H")
    call to_history("n3", nt, "Dat/data-vz.H")
    call to_history("o3", 0., "Dat/data-vz.H")
    call to_history("d3", dt, "Dat/data-vz.H")
         
    call sep_write(data1,"Dat/data-vz.H")

    write(0,*)"============================================================="
    write(0,*)"BEGIN OF SCALING OF ADJOINT BY Theta'"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

    !data 3 adj
    Vx_xB_adj3       =0.0
    Vx_zB_adj3       =0.0
    Vz_xB_adj3       =0.0
    Vz_zB_adj3       =0.0
    sigmaxx_xB_adj3  =0.0
    sigmaxx_zB_adj3  =0.0
    sigmazz_xB_adj3  =0.0
    sigmazz_zB_adj3  =0.0
    sigmaxz_xB_adj3  =0.0
    sigmaxz_zB_adj3  =0.0

    !it=1 until it=nt
    do it=1,nt-1
      call T(.true.,.false.,Vx_xB_adj3(:,:,it:(it+1)),Vx_xB_adj2(:,:,it))
      call T(.true.,.false.,Vx_zB_adj3(:,:,it:(it+1)),Vx_zB_adj2(:,:,it))
      call T(.true.,.false.,Vz_xB_adj3(:,:,it:(it+1)),Vz_xB_adj2(:,:,it))
      call T(.true.,.false.,Vz_zB_adj3(:,:,it:(it+1)),Vz_zB_adj2(:,:,it))
      call T(.true.,.false.,sigmaxx_xB_adj3(:,:,it:(it+1)),sigmaxx_xB_adj2(:,:,it))
      call T(.true.,.false.,sigmaxx_zB_adj3(:,:,it:(it+1)),sigmaxx_zB_adj2(:,:,it))
      call T(.true.,.false.,sigmazz_xB_adj3(:,:,it:(it+1)),sigmazz_xB_adj2(:,:,it))
      call T(.true.,.false.,sigmazz_zB_adj3(:,:,it:(it+1)),sigmazz_zB_adj2(:,:,it))
      call T(.true.,.false.,sigmaxz_xB_adj3(:,:,it:(it+1)),sigmaxz_xB_adj2(:,:,it))
      call T(.true.,.false.,sigmaxz_zB_adj3(:,:,it:(it+1)),sigmaxz_zB_adj2(:,:,it))

      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif
    enddo

    write(0,'(a)',advance="yes")"]"

    data1 = 0.
    data1 = Vx_xB_adj3 + Vx_zB_adj3

    call to_history("n1", nz, "Dat/data2-vx.H")
    call to_history("o1", 0., "Dat/data2-vx.H")
    call to_history("d1", dz, "Dat/data2-vx.H")
    call to_history("n2", nx, "Dat/data2-vx.H")
    call to_history("o2", 0., "Dat/data2-vx.H")
    call to_history("d2", dx, "Dat/data2-vx.H")
    call to_history("n3", nt, "Dat/data2-vx.H")
    call to_history("o3", 0., "Dat/data2-vx.H")
    call to_history("d3", dt, "Dat/data2-vx.H")
         
    call sep_write(data1,"Dat/data2-vx.H")

    data1 = 0.
    data1 = Vz_xB_adj3 + Vz_zB_adj3

    call to_history("n1", nz, "Dat/data2-vz.H")
    call to_history("o1", 0., "Dat/data2-vz.H")
    call to_history("d1", dz, "Dat/data2-vz.H")
    call to_history("n2", nx, "Dat/data2-vz.H")
    call to_history("o2", 0., "Dat/data2-vz.H")
    call to_history("d2", dx, "Dat/data2-vz.H")
    call to_history("n3", nt, "Dat/data2-vz.H")
    call to_history("o3", 0., "Dat/data2-vz.H")
    call to_history("d3", dt, "Dat/data2-vz.H")
         
    call sep_write(data1,"Dat/data2-vz.H")

!###################################################
! Apply imaging condition
!###################################################
    write(0,*)"============================================================="
    write(0,*)"Imaging condition"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["

    r_rhox = 0.
    r_rhoz = 0.
    r_l2mix = 0.
    r_l2miz = 0.
    r_lx = 0.
    r_lz = 0.
    r_mi = 0.

    do it=1,nt
      r_rhox(:,:) = r_rhox(:,:) + (Vx_x0(:,:,it)+Vx_z0(:,:,it))*&
                                  (Vx_xB_adj3(:,:,it)+Vx_zB_adj3(:,:,it))
      r_rhoz(:,:) = r_rhoz(:,:) + (Vz_x0(:,:,it)+Vz_z0(:,:,it))*&
                                  (Vz_xB_adj3(:,:,it)+Vz_zB_adj3(:,:,it))
      r_l2mix(:,:) = r_l2mix(:,:) + (sigmaxx_x0(:,:,it)+sigmaxx_z0(:,:,it))*&
                                    (sigmaxx_xB_adj3(:,:,it)+sigmaxx_zB_adj3(:,:,it))
      r_l2miz(:,:) = r_l2miz(:,:) + (sigmazz_x0(:,:,it)+sigmazz_z0(:,:,it))*&
                                    (sigmazz_xB_adj3(:,:,it)+sigmazz_zB_adj3(:,:,it))
      r_lx(:,:) = r_lx(:,:) + (sigmaxx_x0(:,:,it)+sigmaxx_z0(:,:,it))*&
                              (sigmazz_xB_adj3(:,:,it)+sigmazz_zB_adj3(:,:,it))
      r_lz(:,:) = r_lz(:,:) + (sigmazz_x0(:,:,it)+sigmazz_z0(:,:,it))*&
                              (sigmaxx_xB_adj3(:,:,it)+sigmaxx_zB_adj3(:,:,it))
      r_mi(:,:) = r_mi(:,:) + (sigmaxz_x0(:,:,it)+sigmaxz_z0(:,:,it))*&
                              (sigmaxz_xB_adj3(:,:,it)+sigmaxz_zB_adj3(:,:,it))
      !Write to terminal
      if (mod(it,(nt/59))==0) then
        write(0,'(a)',advance="no")"#"
        close(0)
      endif
    enddo
    write(0,'(a)',advance="yes")"]"

    data0 = 0.
    data0 = r_rhox

    call to_history("n1", nz, "Dat/r-rhox.H")
    call to_history("o1", 0., "Dat/r-rhox.H")
    call to_history("d1", dz, "Dat/r-rhox.H")
    call to_history("n2", nx, "Dat/r-rhox.H")
    call to_history("o2", 0., "Dat/r-rhox.H")
    call to_history("d2", dx, "Dat/r-rhox.H")
         
    call sep_write(data0,"Dat/r-rhox.H")

    data0 = 0.
    data0 = r_rhoz

    call to_history("n1", nz, "Dat/r-rhoz.H")
    call to_history("o1", 0., "Dat/r-rhoz.H")
    call to_history("d1", dz, "Dat/r-rhoz.H")
    call to_history("n2", nx, "Dat/r-rhoz.H")
    call to_history("o2", 0., "Dat/r-rhoz.H")
    call to_history("d2", dx, "Dat/r-rhoz.H")
         
    call sep_write(data0,"Dat/r-rhoz.H")

    data0 = 0.
    data0 = r_l2mix

    call to_history("n1", nz, "Dat/r-l2mix.H")
    call to_history("o1", 0., "Dat/r-l2mix.H")
    call to_history("d1", dz, "Dat/r-l2mix.H")
    call to_history("n2", nx, "Dat/r-l2mix.H")
    call to_history("o2", 0., "Dat/r-l2mix.H")
    call to_history("d2", dx, "Dat/r-l2mix.H")
         
    call sep_write(data0,"Dat/r-l2mix.H")

    data0 = 0.
    data0 = r_l2miz

    call to_history("n1", nz, "Dat/r-l2miz.H")
    call to_history("o1", 0., "Dat/r-l2miz.H")
    call to_history("d1", dz, "Dat/r-l2miz.H")
    call to_history("n2", nx, "Dat/r-l2miz.H")
    call to_history("o2", 0., "Dat/r-l2miz.H")
    call to_history("d2", dx, "Dat/r-l2miz.H")
         
    call sep_write(data0,"Dat/r-l2miz.H")

    data0 = 0.
    data0 = r_lx

    call to_history("n1", nz, "Dat/r-lx.H")
    call to_history("o1", 0., "Dat/r-lx.H")
    call to_history("d1", dz, "Dat/r-lx.H")
    call to_history("n2", nx, "Dat/r-lx.H")
    call to_history("o2", 0., "Dat/r-lx.H")
    call to_history("d2", dx, "Dat/r-lx.H")
         
    call sep_write(data0,"Dat/r-lx.H")

    data0 = 0.
    data0 = r_lz

    call to_history("n1", nz, "Dat/r-lz.H")
    call to_history("o1", 0., "Dat/r-lz.H")
    call to_history("d1", dz, "Dat/r-lz.H")
    call to_history("n2", nx, "Dat/r-lz.H")
    call to_history("o2", 0., "Dat/r-lz.H")
    call to_history("d2", dx, "Dat/r-lz.H")
         
    call sep_write(data0,"Dat/r-lz.H")

    data0 = 0.
    data0 = r_mi

    call to_history("n1", nz, "Dat/r-mi.H")
    call to_history("o1", 0., "Dat/r-mi.H")
    call to_history("d1", dz, "Dat/r-mi.H")
    call to_history("n2", nx, "Dat/r-mi.H")
    call to_history("o2", 0., "Dat/r-mi.H")
    call to_history("d2", dx, "Dat/r-mi.H")
         
    call sep_write(data0,"Dat/r-mi.H")

  end subroutine

end module
