! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Mar, 24th, 2015

! finite difference module
module rtmlib_mod

  use sep
  use omp_lib
  use rtmlib_modF

  implicit none
  integer, private   :: nx, nz, iz, ix
  real, private      :: dx, dt, aux
  real, private      :: c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0
  real, pointer, dimension(:,:),private :: m1,m1_x,m1_z,m2,m3
  real, allocatable, dimension(:,:), private :: aux_m2,aux_m3,aux_m2m3
  logical,private    :: lame

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine rtmmod_init(nz_in,nx_in,dx_in, dt_in, m1_x_in,m1_z_in,m2_in,m3_in)
    integer        :: nx_in, nz_in
    real           :: dx_in, dt_in
    real, dimension(:,:),target :: m1_x_in, m1_z_in, m2_in, m3_in
    logical        :: lame_in
    !OMP settings
!    integer                               :: node, nodes

    !OMP
    !nodes=16
    !node = omp_get_num_procs()
    !call omp_set_num_threads(nodes)

    nx = nx_in
    nz = nz_in
    dx = dx_in
    dt = dt_in
 !   m1 => m1_in
    m1_x => m1_x_in
    m1_z => m1_z_in
    m2 => m2_in
    m3 => m3_in
 !   lame = lame_in

    allocate(aux_m2(nz,nx))
    allocate(aux_m3(nz,nx))
    allocate(aux_m2m3(nz,nx))

    aux = dt/dx

    m1_x = aux*m1_x
    m1_z = aux*m1_z
!   if (lame) then
      aux_m2 = aux*m2
      aux_m3 = aux*m3
      aux_m2m3 = aux*(m2+2.0*m3)
!    else
!      aux_m2 = aux*m1*(m2*m2-2.0*m3*m3)
!      aux_m3 = aux*m1*m3*m3
!      aux_m2m3 = aux*m1*m2*m2
!    endif

  endsubroutine

  subroutine Test(flag)
  integer flag
  write(0,*) "This is a test module"
  end subroutine
    
  subroutine rtm_op3(nt, nz, nx, zrec, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, Vx, Vz, sigmaxx, sigmazz, sigmaxz, m1_x, m1_z,aux_m2_c, aux_m3_c, aux_m2m3_c)
  
   integer :: zrec, nt, nz, nx, it
!  real, pointer, dimension(:,:),private :: Vx_in, Vz_in, sigmaxx_in, sigmazz_in, sigmaxz_in
   real,dimension(nt,nx),intent(inout)   :: Vx, Vz, sigmaxx, sigmazz, sigmaxz
   real,dimension(nz,nx)   :: m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c
   real,dimension(nz, nx,nt),intent(inout)   :: Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0! model,data,m1_x
  real*8                                :: start_time_rtm,end_time_rtm
 
    write(0,*)"============================================================="
    write(0,*)"BEGIN OF ADJOINT PROPAGATION OF RECEIVER DATA"
    write(0,*)"============================================================="
    write(0,*)"Processing:"
    write(0,*)"0%           25%            50%           75%            100%"
    write(0,'(a)',advance="no")"["
!    call rtm_op3_c(nt, nz, nx, zrec, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, Vx, Vz, sigmaxx, sigmazz, sigmaxz, m1_x, m1_z,aux_m2_c, aux_m3_c, aux_m2m3_c)
 !   call Test(1)
!
    !Zero arrays
    Vx0       = 0.
    Vz0       = 0.
    sigmaxx0  = 0.
    sigmazz0  = 0.
    sigmaxz0  = 0.

!    ! Last time steps
!    !it=nt
    Vx0(zrec,:,nt) = Vx(nt,:)
   ! call Test(1)
    Vz0(zrec,:,nt) = Vz(nt,:)
  !  call Test(1)
    sigmaxx0(zrec,:,nt) = sigmaxx(nt,:)
    sigmazz0(zrec,:,nt) = sigmazz(nt,:)
    sigmaxz0(zrec,:,nt) = sigmaxz(nt,:)
!
!    !it=nt-1
!!    call K(Vx0(:,:,nt))
!!    call K(Vz0(:,:,nt))
!!    call K(sigmaxx0(:,:,nt))
!!    call K(sigmazz0(:,:,nt))
!!    call K(sigmaxz0(:,:,nt))
!!    !Vx
   ! call Test(1)
    Vx0(zrec,:,nt-1) = Vx((nt-1),:)
    !write (0,*), "It goes here in Fortran"
    call D_c(nz, nx, Vx0(:,:,(nt-1)),sigmaxx0(:,:,nt), aux_m2m3_c(:,:))    
    call F_c(nz, nx, Vx0(:,:,(nt-1)),sigmazz0(:,:,nt), aux_m2_c(:,:))
    call H_c(nz, nx, Vx0(:,:,(nt-1)),sigmaxz0(:,:,nt), aux_m3_c(:,:))
    !Vz
    Vz0(zrec,:,nt-1) = Vz((nt-1),:)
    call E_c(nz, nx, Vz0(:,:,(nt-1)),sigmaxx0(:,:,nt), aux_m2_c(:,:))    
    call G_c(nz, nx, Vz0(:,:,(nt-1)),sigmazz0(:,:,nt), aux_m2m3_c(:,:))
    call J_c(nz, nx, Vz0(:,:,(nt-1)),sigmaxz0(:,:,nt), aux_m3_c(:,:))
    !sigmaxx
    sigmaxx0(zrec,:,nt-1) = sigmaxx((nt-1),:)
    call Ax_c(nz, nx, sigmaxx0(:,:,(nt-1)),Vx0(:,:,nt), m1_x(:,:))
    !sigmazz
    sigmazz0(zrec,:,nt-1) = sigmazz((nt-1),:)
    call Bz_c(nz, nx, sigmazz0(:,:,(nt-1)),Vz0(:,:,nt), m1_z(:,:))
    !sigmaxz
    sigmaxz0(zrec,:,nt-1) = sigmaxz((nt-1),:)
    call Bx_c(nz, nx, sigmaxz0(:,:,(nt-1)),Vx0(:,:,nt), m1_x(:,:))
    call Az_c(nz, nx, sigmaxz0(:,:,(nt-1)),Vz0(:,:,nt), m1_z(:,:))
   ! call Test(1)

    do it=nt-2,1,-1
      Vx0(zrec,:,it) = Vx(it,:)
      Vz0(zrec,:,it) = Vz(it,:)
      sigmaxx0(zrec,:,it) = sigmaxx(it,:)
      sigmazz0(zrec,:,it) = sigmazz(it,:)
      sigmaxz0(zrec,:,it) = sigmaxz(it,:)
!      call Az_c(nz, nx, sigmaxz0(:,:,it),Vz0(:,:,(it+1)), m1_z(:,:))


    enddo   

    
    !it=nt-2 until it=1
!    do it=nt-2,1,-1

!      Vx0(:,:,it) =  Vx0(:,:,it) + Vx0(:,:,(it+2))
!      Vz0(:,:,it) = Vz0(:,:,it) + Vz0(:,:,(it+2))
!      sigmaxx0(:,:,it) = sigmaxx0(:,:,it) + sigmaxx0(:,:,(it+2))
!      sigmazz0(:,:,it) = sigmazz0(:,:,it) + sigmazz0(:,:,(it+2))
!      sigmaxz0(:,:,it) = sigmaxz0(:,:,it) + sigmaxz0(:,:,(it+2))
  start_time_rtm = omp_get_wtime()
  call rtm_op3_c(nt, nz, nx, zrec,Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0, m1_x, m1_z, aux_m2_c, aux_m3_c, aux_m2m3_c)
  end_time_rtm = omp_get_wtime()
  write(0,*)"============================================================="
  write(0,*)"RTM ELAPSED TIME, MEASURED IN (RTMLIB-MOD.F90)"
  write(0,*)"Total=",end_time_rtm - start_time_rtm, "seconds."
  write(0,*)"============================================================="

   
!      call D_c(nz, nx, Vx0(:,:,it),sigmaxx0(:,:,(it+1)), aux_m2m3_c(:,:))    
!      call F_c(nz, nx, Vx0(:,:,it),sigmazz0(:,:,(it+1)), aux_m2_c(:,:))
!      call H_c(nz, nx, Vx0(:,:,it),sigmaxz0(:,:,(it+1)), aux_m3_c(:,:))
      !Vz
      !Vz0(zrec,:,it) = Vz(it,:)
!      call E_c(nz, nx, Vz0(:,:,it),sigmaxx0(:,:,(it+1)), aux_m2_c(:,:))    
!      call G_c(nz, nx, Vz0(:,:,it),sigmazz0(:,:,(it+1)), aux_m2m3_c(:,:))
!      call J_c(nz, nx, Vz0(:,:,it),sigmaxz0(:,:,(it+1)), aux_m3_c(:,:))
      !sigmaxx
      !sigmaxx0(zrec,:,it) = sigmaxx(it,:)
      !write(0,*) "It goes here in 00000 Fortran"

!      call Ax_c(nz, nx, sigmaxx0(:,:,it),Vx0(:,:,(it+1)), m1_x(:,:))
      !sigmazz
      !sigmazz0(zrec,:,it) = sigmazz(it,:)
!      call Bz_c(nz, nx, sigmazz0(:,:,it),Vz0(:,:,(it+1)), m1_z)
      !sigmaxz
      !sigmaxz0(zrec,:,it) = sigmaxz(it,:)
!      call Bx_c(nz, nx, sigmaxz0(:,:,it),Vx0(:,:,(it+1)), m1_x(:,:))
!      call Az_c(nz, nx, sigmaxz0(:,:,it),Vz0(:,:,(it+1)), m1_z(:,:))

      !Write to terminal
!      if (mod(it,(nt/59))==0) then
!        write(0,'(a)',advance="no")"#"
!        close(0)
!      endif
!
!    enddo
!    write(0,'(a)',advance="yes")"]"

endsubroutine

!=============================================================================
!!!! 10TH ORDER OPERATOR
  subroutine Ax(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data
    integer :: ith

    if(adj)then
      if(.not.add) model = 0.
 !     write(0,*) m1_x(15,15)

      do ix=6,nx-5
!!$OMP PARALLEL DO SHARED(model,data,m1_x,ix) PRIVATE(ix,iz)
!        write(0,*) omp_get_thread_num()
        do iz=5,nz-5
          model(iz,ix+4) = model(iz,ix+4) + m1_x(iz,ix)*c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + m1_x(iz,ix)*c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + m1_x(iz,ix)*c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + m1_x(iz,ix)*c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + m1_x(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - m1_x(iz,ix)*c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - m1_x(iz,ix)*c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - m1_x(iz,ix)*c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - m1_x(iz,ix)*c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - m1_x(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data,m1_x,c1,c2,c3,c4,c5,nx,nz) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
        data(iz,ix) = data(iz,ix) + m1_x(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
                                                 c2*(model(iz,ix+3) - model(iz,ix-4))+&
                                                 c3*(model(iz,ix+2) - model(iz,ix-3))+&
                                                 c4*(model(iz,ix+1) - model(iz,ix-2))+&
                                                 c5*(model(iz,ix  ) - model(iz,ix-1)))        
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

 
  subroutine Bx(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data,m1_x) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + m1_x(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + m1_x(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + m1_x(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + m1_x(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + m1_x(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - m1_x(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - m1_x(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - m1_x(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - m1_x(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - m1_x(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data,m1_x) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + m1_x(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                   c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                   c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                   c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                   c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine Az(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          model(iz,ix+5) = model(iz,ix+5) + m1_z(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + m1_z(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + m1_z(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + m1_z(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + m1_z(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - m1_z(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - m1_z(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - m1_z(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - m1_z(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - m1_z(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + m1_z(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                   c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                   c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                   c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                   c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Bz(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          model(iz+4,ix) = model(iz+4,ix) + m1_z(iz,ix)*c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + m1_z(iz,ix)*c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + m1_z(iz,ix)*c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + m1_z(iz,ix)*c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + m1_z(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - m1_z(iz,ix)*c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - m1_z(iz,ix)*c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - m1_z(iz,ix)*c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - m1_z(iz,ix)*c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - m1_z(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + m1_z(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                                   c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                                   c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                                   c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                                   c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine Cx(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          model(iz,ix) = model(iz,ix) + dt*m1_x(iz,ix)*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          data(iz,ix) = data(iz,ix) + dt*m1_x(iz,ix)*model(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Cz(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          model(iz,ix) = model(iz,ix) + dt*m1_z(iz,ix)*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          data(iz,ix) = data(iz,ix) + dt*m1_z(iz,ix)*model(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine D(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz,ix+5) = model(iz,ix+5) + aux_m2m3(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + aux_m2m3(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m2m3(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m2m3(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m2m3(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - aux_m2m3(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m2m3(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m2m3(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m2m3(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m2m3(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + aux_m2m3(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                       c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                       c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                       c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                       c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine E(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + aux_m2(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + aux_m2(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m2(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m2(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m2(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - aux_m2(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m2(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m2(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m2(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m2(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + aux_m2(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                     c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                     c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                     c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                     c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine F(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz,ix+5) = model(iz,ix+5) + aux_m2(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + aux_m2(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m2(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m2(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m2(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - aux_m2(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m2(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m2(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m2(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m2(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + aux_m2(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                     c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                     c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                     c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                     c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine G(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + aux_m2m3(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + aux_m2m3(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m2m3(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m2m3(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m2m3(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - aux_m2m3(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m2m3(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m2m3(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m2m3(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m2m3(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + aux_m2m3(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                       c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                       c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                       c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                       c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine H(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          model(iz+4,ix) = model(iz+4,ix) + aux_m3(iz,ix)*c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m3(iz,ix)*c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m3(iz,ix)*c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m3(iz,ix)*c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + aux_m3(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m3(iz,ix)*c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m3(iz,ix)*c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m3(iz,ix)*c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m3(iz,ix)*c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - aux_m3(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + aux_m3(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                                     c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                                     c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                                     c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                                     c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine J(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj)then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          model(iz,ix+4) = model(iz,ix+4) + aux_m3(iz,ix)*c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m3(iz,ix)*c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m3(iz,ix)*c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m3(iz,ix)*c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + aux_m3(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m3(iz,ix)*c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m3(iz,ix)*c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m3(iz,ix)*c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m3(iz,ix)*c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - aux_m3(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
        data(iz,ix) = data(iz,ix) + aux_m3(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
                                                   c2*(model(iz,ix+3) - model(iz,ix-4))+&
                                                   c3*(model(iz,ix+2) - model(iz,ix-3))+&
                                                   c4*(model(iz,ix+1) - model(iz,ix-2))+&
                                                   c5*(model(iz,ix  ) - model(iz,ix-1)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine
 
  subroutine Ax_c(nz, nx, model,data,m1_x)
    !logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data,m1_x
    integer :: nz, nx

      do ix=6,nx-5
!!$OMP PARALLEL DO SHARED(model,data,m1_x,ix) PRIVATE(ix,iz)
!        write(0,*) omp_get_thread_num()
        do iz=5,nz-5
          model(iz,ix+4) = model(iz,ix+4) + m1_x(iz,ix)*c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + m1_x(iz,ix)*c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + m1_x(iz,ix)*c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + m1_x(iz,ix)*c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + m1_x(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - m1_x(iz,ix)*c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - m1_x(iz,ix)*c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - m1_x(iz,ix)*c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - m1_x(iz,ix)*c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - m1_x(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
  endsubroutine

  subroutine Bx_c(nz, nx, model,data, m1_x)
    !logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, m1_x
    !write(0,*) "Bx"
    integer :: nz, nx

      do ix=6,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + m1_x(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + m1_x(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + m1_x(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + m1_x(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + m1_x(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - m1_x(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - m1_x(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - m1_x(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - m1_x(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - m1_x(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
  endsubroutine

  subroutine Az_c(nz, nx, model,data, m1_z)
    !logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, m1_z
    !write(0,*) "Az"
    integer :: nz, nx

!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          model(iz,ix+5) = model(iz,ix+5) + m1_z(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + m1_z(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + m1_z(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + m1_z(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + m1_z(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - m1_z(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - m1_z(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - m1_z(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - m1_z(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - m1_z(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    !else
    !  if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
    !  do ix=5,nx-5
    !    do iz=6,nz-5
    !      data(iz,ix) = data(iz,ix) + m1_z(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
    !                                               c2*(model(iz,ix+4) - model(iz,ix-3))+&
    !                                               c3*(model(iz,ix+3) - model(iz,ix-2))+&
    !                                               c4*(model(iz,ix+2) - model(iz,ix-1))+&
    !                                               c5*(model(iz,ix+1) - model(iz,ix  )))
    !    enddo
    !  enddo
!!!$OMP END PARALLEL DO
    !endif
  endsubroutine

  subroutine Bz_c(nz, nx, model,data, m1_z)
   ! logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, m1_z
    !write(0,*) "Bz"
    integer :: nz, nx

    !if(adj) then
      !if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          model(iz+4,ix) = model(iz+4,ix) + m1_z(iz,ix)*c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + m1_z(iz,ix)*c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + m1_z(iz,ix)*c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + m1_z(iz,ix)*c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + m1_z(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - m1_z(iz,ix)*c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - m1_z(iz,ix)*c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - m1_z(iz,ix)*c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - m1_z(iz,ix)*c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - m1_z(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
   !else
   !  if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
   !  do ix=5,nx-5
   !    do iz=6,nz-5
   !      data(iz,ix) = data(iz,ix) + m1_z(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
   !                                               c2*(model(iz+3,ix) - model(iz-4,ix))+&
   !                                               c3*(model(iz+2,ix) - model(iz-3,ix))+&
   !                                               c4*(model(iz+1,ix) - model(iz-2,ix))+&
   !                                               c5*(model(iz  ,ix) - model(iz-1,ix)))
   !    enddo
   !  enddo
!!$OMP END PARALLEL DO
   !endif 
  endsubroutine


  subroutine D_c(nz, nx, model,data, aux_m2m3_c)
    !logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, aux_m2m3_c 
    !write(0,*) "Dc"
    integer :: nz, nx

    !if(adj) then
     ! if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz,ix+5) = model(iz,ix+5) + aux_m2m3_c(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + aux_m2m3_c(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m2m3_c(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m2m3_c(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m2m3_c(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - aux_m2m3_c(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m2m3_c(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m2m3_c(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m2m3_c(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m2m3_c(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
   ! else
   !   if(.not.add) data = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
   !   do ix=5,nx-5
   !     do iz=5,nz-5
   !       data(iz,ix) = data(iz,ix) + aux_m2m3(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
   !                                                    c2*(model(iz,ix+4) - model(iz,ix-3))+&
   !                                                    c3*(model(iz,ix+3) - model(iz,ix-2))+&
   !                                                    c4*(model(iz,ix+2) - model(iz,ix-1))+&
   !                                                    c5*(model(iz,ix+1) - model(iz,ix  )))
   !     enddo
   !   enddo
!!!!$OMP END PARALLEL DO
   ! endif
  endsubroutine

  subroutine E_c(nz, nx, model,data, aux_m2_c)
   ! logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, aux_m2_c
    !write(0,*) "Ec"
    integer :: nz, nx

    !if(adj) then
      !if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + aux_m2_c(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + aux_m2_c(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m2_c(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m2_c(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m2_c(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - aux_m2_c(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m2_c(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m2_c(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m2_c(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m2_c(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
  !  else
  !    if(.not.add) data = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
  !    do ix=5,nx-5
  !      do iz=5,nz-5
  !        data(iz,ix) = data(iz,ix) + aux_m2(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
  !                                                   c2*(model(iz+4,ix) - model(iz-3,ix))+&
  !                                                   c3*(model(iz+3,ix) - model(iz-2,ix))+&
  !                                                   c4*(model(iz+2,ix) - model(iz-1,ix))+&
  !                                                   c5*(model(iz+1,ix) - model(iz  ,ix)))
  !      enddo
  !    enddo
!!!!$OMP END PARALLEL DO
  !  endif
  endsubroutine

  subroutine F_c(nz, nx, model,data, aux_m2_c)
    !logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, aux_m2_c
    !write(0,*) "Fc"
    integer :: nz, nx

   ! if(adj) then
      !if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz,ix+5) = model(iz,ix+5) + aux_m2_c(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + aux_m2_c(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m2_c(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m2_c(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m2_c(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - aux_m2_c(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m2_c(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m2_c(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m2_c(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m2_c(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
   ! else
   !   if(.not.add) data = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
   !   do ix=5,nx-5
   !     do iz=5,nz-5
   !       data(iz,ix) = data(iz,ix) + aux_m2(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
   !                                                  c2*(model(iz,ix+4) - model(iz,ix-3))+&
   !                                                  c3*(model(iz,ix+3) - model(iz,ix-2))+&
   !                                                  c4*(model(iz,ix+2) - model(iz,ix-1))+&
   !                                                  c5*(model(iz,ix+1) - model(iz,ix  )))
   !     enddo
   !   enddo
!!!!$OMP END PARALLEL DO
   ! endif
  endsubroutine

  subroutine G_c(nz, nx, model,data, aux_m2m3_c)
    !logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, aux_m2m3_c
    !write(0,*) "Gc"
    integer :: nz, nx

    !if(adj) then
     ! if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + aux_m2m3_c(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + aux_m2m3_c(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m2m3_c(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m2m3_c(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m2m3_c(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - aux_m2m3_c(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m2m3_c(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m2m3_c(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m2m3_c(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m2m3_c(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
   ! else
   !   if(.not.add) data = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
   !   do ix=5,nx-5
   !     do iz=5,nz-5
   !       data(iz,ix) = data(iz,ix) + aux_m2m3(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
   !                                                    c2*(model(iz+4,ix) - model(iz-3,ix))+&
   !                                                    c3*(model(iz+3,ix) - model(iz-2,ix))+&
   !                                                    c4*(model(iz+2,ix) - model(iz-1,ix))+&
   !                                                    c5*(model(iz+1,ix) - model(iz  ,ix)))
   !     enddo
   !   enddo
!!!!$OMP END PARALLEL DO
   ! endif 
  endsubroutine

  subroutine H_c(nz, nx, model,data, aux_m3_c)
    !logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, aux_m3_c
    !write(0,*) "Hc"
    integer :: nz, nx

    !if(adj) then
     ! if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          model(iz+4,ix) = model(iz+4,ix) + aux_m3_c(iz,ix)*c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m3_c(iz,ix)*c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m3_c(iz,ix)*c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m3_c(iz,ix)*c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + aux_m3_c(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m3_c(iz,ix)*c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m3_c(iz,ix)*c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m3_c(iz,ix)*c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m3_c(iz,ix)*c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - aux_m3_c(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
   ! else
   !   if(.not.add) data = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
   !   do ix=6,nx-5
   !     do iz=6,nz-5
   !       data(iz,ix) = data(iz,ix) + aux_m3(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
   !                                                  c2*(model(iz+3,ix) - model(iz-4,ix))+&
   !                                                  c3*(model(iz+2,ix) - model(iz-3,ix))+&
   !                                                  c4*(model(iz+1,ix) - model(iz-2,ix))+&
   !                                                  c5*(model(iz  ,ix) - model(iz-1,ix)))
   !     enddo
   !   enddo
!!!!$OMP END PARALLEL DO
   ! endif
  endsubroutine

  subroutine J_c(nz, nx, model,data, aux_m3_c)
    !logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data, aux_m3_c
    !write(0,*) "Hc"
    integer :: nz, nx

    !if(adj)then
     ! if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          model(iz,ix+4) = model(iz,ix+4) + aux_m3_c(iz,ix)*c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m3_c(iz,ix)*c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m3_c(iz,ix)*c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m3_c(iz,ix)*c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + aux_m3_c(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m3_c(iz,ix)*c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m3_c(iz,ix)*c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m3_c(iz,ix)*c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m3_c(iz,ix)*c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - aux_m3_c(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
   ! else
   !   if(.not.add) data = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
   !   do ix=6,nx-5
   !     do iz=6,nz-5
   !     data(iz,ix) = data(iz,ix) + aux_m3(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
   !                                                c2*(model(iz,ix+3) - model(iz,ix-4))+&
   !                                                c3*(model(iz,ix+2) - model(iz,ix-3))+&
   !                                                c4*(model(iz,ix+1) - model(iz,ix-2))+&
   !                                                c5*(model(iz,ix  ) - model(iz,ix-1)))
   !     enddo
   !   enddo
!!!!$OMP END PARALLEL DO
   ! endif
  endsubroutine

end module
