! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Mar, 24th, 2015

!Modified by Lin Gan at 
! the Stanford Exploration Project
! Jul, 15th, 2015


! Fortran Interfaces to C
module rtmlib_modF

  use iso_c_binding


  
  implicit none
!=============================================================================
!!! 10TH ORDER OPERATOR

!  subroutine Test(flag)
!  integer flag
!  write(0,*) "This is a test module"
!  end subroutine


interface
 
  subroutine rtm_op3_c(nt, nz, nx, zrec, Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0,  m1_x, m1_z,aux_m2_c, aux_m3_c, aux_m2m3_c) bind(c, name='rtm_op3_c_')
  import
  integer(c_int), value        :: nt, nz, nx, zrec
  real(c_float),dimension(*)   :: Vx0, Vz0, sigmaxx0, sigmazz0, sigmaxz0 
!  real(c_float),dimension(*)   :: Vx, Vz, sigmaxx, sigmazz, sigmaxz
  real(c_float),dimension(*)   :: m1_x, m1_z,aux_m2_c, aux_m3_c, aux_m2m3_c 
  write(0,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" 
  write(0,*) "This is the output from the Fortran interface parts" 
!      do ix=6,nx-5
!        do iz=5,nz-5
!          model(iz,ix+4) = model(iz,ix+4) + m1_x(iz,ix)*c1*data(iz,ix)
!          model(iz,ix+3) = model(iz,ix+3) + m1_x(iz,ix)*c2*data(iz,ix)
!          model(iz,ix+2) = model(iz,ix+2) + m1_x(iz,ix)*c3*data(iz,ix)
!          model(iz,ix+1) = model(iz,ix+1) + m1_x(iz,ix)*c4*data(iz,ix)
!          model(iz,ix  ) = model(iz,ix  ) + m1_x(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-1) = model(iz,ix-1) - m1_x(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-2) = model(iz,ix-2) - m1_x(iz,ix)*c4*data(iz,ix)
!          model(iz,ix-3) = model(iz,ix-3) - m1_x(iz,ix)*c3*data(iz,ix)
!          model(iz,ix-4) = model(iz,ix-4) - m1_x(iz,ix)*c2*data(iz,ix)
!          model(iz,ix-5) = model(iz,ix-5) - m1_x(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo

     endsubroutine



!  subroutine Ax_c(nz, nx, model,data,m1_x) bind(c, name='Ax_c_')
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,m1_x
!  write(0,*) "FFFFFF" 
!!      do ix=6,nx-5
!!        do iz=5,nz-5
!!          model(iz,ix+4) = model(iz,ix+4) + m1_x(iz,ix)*c1*data(iz,ix)
!!          model(iz,ix+3) = model(iz,ix+3) + m1_x(iz,ix)*c2*data(iz,ix)
!!          model(iz,ix+2) = model(iz,ix+2) + m1_x(iz,ix)*c3*data(iz,ix)
!!          model(iz,ix+1) = model(iz,ix+1) + m1_x(iz,ix)*c4*data(iz,ix)
!!          model(iz,ix  ) = model(iz,ix  ) + m1_x(iz,ix)*c5*data(iz,ix)
!!          model(iz,ix-1) = model(iz,ix-1) - m1_x(iz,ix)*c5*data(iz,ix)
!!          model(iz,ix-2) = model(iz,ix-2) - m1_x(iz,ix)*c4*data(iz,ix)
!!          model(iz,ix-3) = model(iz,ix-3) - m1_x(iz,ix)*c3*data(iz,ix)
!!          model(iz,ix-4) = model(iz,ix-4) - m1_x(iz,ix)*c2*data(iz,ix)
!!          model(iz,ix-5) = model(iz,ix-5) - m1_x(iz,ix)*c1*data(iz,ix)
!!        enddo
!!      enddo
!
!     endsubroutine
!
!  subroutine Bx_c(nz, nx, model,data,m1_x) bind(c, name='Bx_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,m1_x
!  
!  endsubroutine
! 
!  subroutine Az_c(nz, nx, model,data,m1_z) bind(c, name='Az_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,m1_z
!  
!  endsubroutine
!  
!  
!  subroutine Bz_c(nz, nx, model,data,m1_z) bind(c, name='Bz_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,m1_z
!  endsubroutine
! 
!  subroutine D_c(nz, nx, model,data,aux_m2m3_c) bind(c, name='D_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,aux_m2m3_c
!  endsubroutine
! 
! 
!  subroutine E_c(nz, nx, model,data,aux_m2_c) bind(c, name='E_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,aux_m2_c
!  endsubroutine
!  
!  subroutine F_c(nz, nx, model,data,aux_m2_c) bind(c, name='F_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,aux_m2_c
!  endsubroutine
!  
!  subroutine G_c(nz, nx, model,data,aux_m2m3_c) bind(c, name='G_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,aux_m2m3_c
!  endsubroutine
!  
!  subroutine H_c(nz, nx, model,data,aux_m3_c) bind(c, name='H_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,aux_m3_c
!  endsubroutine
!
!   subroutine J_c(nz, nx, model,data,aux_m3_c) bind(c, name='J_c_')
!  
!  import
!  integer(c_int), value        :: nz, nx
!  real(c_float),dimension(*)    :: model,data,aux_m3_c
!  endsubroutine
!  
! 
!
!
!

 endinterface 
 
 !subroutine PPP(nz)
 ! integer ::nz
 ! write(0,*) "adsfadsf"
 ! endsubroutine



 
 
 


!  use sep
!  use omp_lib
!
!  implicit none
!  integer, private   :: nx, nz, iz, ix
!  real, private      :: dx, dt, aux
!  real, private      :: c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0
!  real, pointer, dimension(:,:),private ::m1_x,m1_z,m2,m3, aux_m2F,aux_m3F,aux_m2m3F
! 
! ! real, allocatable, dimension(:,:), private :: aux_m2F,aux_m3F,aux_m2m3F
! ! logical,private    :: lame
!
!  contains
!
!  subroutine rtmmod_initF(nz_in,nx_in,dx_in, dt_in, m2_in,m3_in)
!    integer        :: nx_in, nz_in
!    real           :: dx_in, dt_in
!    real, dimension(:,:),target :: m2_in, m3_in
! !   logical        :: lame_in
!    !OMP settings
!!    integer                               :: node, nodes
!
!    !OMP
!    !nodes=16
!    !node = omp_get_num_procs()
!    !call omp_set_num_threads(nodes)
!
!    nx = nx_in
!    nz = nz_in
!    dx = dx_in
!    dt = dt_in
! !   m1 => m1_in
!  !  m1_x => m1_x_in
! !   m1_z => m1_z_in
!    m2 => m2_in
!    m3 => m3_in
! !   lame = lame_in
!
!    allocate(aux_m2F(nz,nx))
!    allocate(aux_m3F(nz,nx))
!    allocate(aux_m2m3F(nz,nx))
!
!    aux = dt/dx
!
!!   m1_x = aux*m1_x
!!   m1_z = aux*m1_z
!!   if (lame) then
!      aux_m2F = aux*m2
!      aux_m3F = aux*m3
!      aux_m2m3F = aux*(m2+2.0*m3)
!!    else
!!      aux_m2 = aux*m1*(m2*m2-2.0*m3*m3)
!!      aux_m3F = aux*m1*m3*m3
!!      aux_m2m3 = aux*m1*m2*m2
!!    endif
!
!  endsubroutine
!
!
!  subroutine AxF_c(nz, nx,model,data,m1_x)
!    !logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data,m1_x
!    integer :: nz, nx
!
!
!      do ix=6,nx-5
!!!$OMP PARALLEL DO SHARED(model,data,m1_x,ix) PRIVATE(ix,iz)
!!        write(0,*) omp_get_thread_num()
!        do iz=5,nz-5
!          model(iz,ix+4) = model(iz,ix+4) + m1_x(iz,ix)*c1*data(iz,ix)
!          model(iz,ix+3) = model(iz,ix+3) + m1_x(iz,ix)*c2*data(iz,ix)
!          model(iz,ix+2) = model(iz,ix+2) + m1_x(iz,ix)*c3*data(iz,ix)
!          model(iz,ix+1) = model(iz,ix+1) + m1_x(iz,ix)*c4*data(iz,ix)
!          model(iz,ix  ) = model(iz,ix  ) + m1_x(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-1) = model(iz,ix-1) - m1_x(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-2) = model(iz,ix-2) - m1_x(iz,ix)*c4*data(iz,ix)
!          model(iz,ix-3) = model(iz,ix-3) - m1_x(iz,ix)*c3*data(iz,ix)
!          model(iz,ix-4) = model(iz,ix-4) - m1_x(iz,ix)*c2*data(iz,ix)
!          model(iz,ix-5) = model(iz,ix-5) - m1_x(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
! endsubroutine
!
!  subroutine BxF_c(nz, nx, model,data, m1_x)
!    !logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data, m1_x
!    !write(0,*) "Bx"
!    integer :: nz, nx
!
!      do ix=6,nx-5
!        do iz=5,nz-5
!          model(iz+5,ix) = model(iz+5,ix) + m1_x(iz,ix)*c1*data(iz,ix)
!          model(iz+4,ix) = model(iz+4,ix) + m1_x(iz,ix)*c2*data(iz,ix)
!          model(iz+3,ix) = model(iz+3,ix) + m1_x(iz,ix)*c3*data(iz,ix)
!          model(iz+2,ix) = model(iz+2,ix) + m1_x(iz,ix)*c4*data(iz,ix)
!          model(iz+1,ix) = model(iz+1,ix) + m1_x(iz,ix)*c5*data(iz,ix)
!          model(iz  ,ix) = model(iz  ,ix) - m1_x(iz,ix)*c5*data(iz,ix)
!          model(iz-1,ix) = model(iz-1,ix) - m1_x(iz,ix)*c4*data(iz,ix)
!          model(iz-2,ix) = model(iz-2,ix) - m1_x(iz,ix)*c3*data(iz,ix)
!          model(iz-3,ix) = model(iz-3,ix) - m1_x(iz,ix)*c2*data(iz,ix)
!          model(iz-4,ix) = model(iz-4,ix) - m1_x(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
! endsubroutine
!
!  subroutine AzF_c(nz, nx, model,data, m1_z)
!    !logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data, m1_z
!    !write(0,*) "Az"
!    integer :: nz, nx
!
!      do ix=5,nx-5
!        do iz=6,nz-5
!          model(iz,ix+5) = model(iz,ix+5) + m1_z(iz,ix)*c1*data(iz,ix)
!          model(iz,ix+4) = model(iz,ix+4) + m1_z(iz,ix)*c2*data(iz,ix)
!          model(iz,ix+3) = model(iz,ix+3) + m1_z(iz,ix)*c3*data(iz,ix)
!          model(iz,ix+2) = model(iz,ix+2) + m1_z(iz,ix)*c4*data(iz,ix)
!          model(iz,ix+1) = model(iz,ix+1) + m1_z(iz,ix)*c5*data(iz,ix)
!          model(iz,ix  ) = model(iz,ix  ) - m1_z(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-1) = model(iz,ix-1) - m1_z(iz,ix)*c4*data(iz,ix)
!          model(iz,ix-2) = model(iz,ix-2) - m1_z(iz,ix)*c3*data(iz,ix)
!          model(iz,ix-3) = model(iz,ix-3) - m1_z(iz,ix)*c2*data(iz,ix)
!          model(iz,ix-4) = model(iz,ix-4) - m1_z(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
! endsubroutine
!
!  subroutine BzF_c(nz, nx, model,data, m1_z)
!    real,dimension(nz,nx),intent(inout)    :: model,data, m1_z
!    integer :: nz, nx
!
!     do ix=5,nx-5
!        do iz=6,nz-5
!          model(iz+4,ix) = model(iz+4,ix) + m1_z(iz,ix)*c1*data(iz,ix)
!          model(iz+3,ix) = model(iz+3,ix) + m1_z(iz,ix)*c2*data(iz,ix)
!          model(iz+2,ix) = model(iz+2,ix) + m1_z(iz,ix)*c3*data(iz,ix)
!          model(iz+1,ix) = model(iz+1,ix) + m1_z(iz,ix)*c4*data(iz,ix)
!          model(iz  ,ix) = model(iz  ,ix) + m1_z(iz,ix)*c5*data(iz,ix)
!          model(iz-1,ix) = model(iz-1,ix) - m1_z(iz,ix)*c5*data(iz,ix)
!          model(iz-2,ix) = model(iz-2,ix) - m1_z(iz,ix)*c4*data(iz,ix)
!          model(iz-3,ix) = model(iz-3,ix) - m1_z(iz,ix)*c3*data(iz,ix)
!          model(iz-4,ix) = model(iz-4,ix) - m1_z(iz,ix)*c2*data(iz,ix)
!          model(iz-5,ix) = model(iz-5,ix) - m1_z(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
! endsubroutine
! 
! subroutine DF_c(nz, nx, model,data)
!    !logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data
!    integer :: nz, nx
!    !if(adj) then
!     ! if(.not.add) model = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!      do ix=5,nx-5
!        do iz=5,nz-5
!          model(iz,ix+5) = model(iz,ix+5) + aux_m2m3F(iz,ix)*c1*data(iz,ix)
!          model(iz,ix+4) = model(iz,ix+4) + aux_m2m3F(iz,ix)*c2*data(iz,ix)
!          model(iz,ix+3) = model(iz,ix+3) + aux_m2m3F(iz,ix)*c3*data(iz,ix)
!          model(iz,ix+2) = model(iz,ix+2) + aux_m2m3F(iz,ix)*c4*data(iz,ix)
!          model(iz,ix+1) = model(iz,ix+1) + aux_m2m3F(iz,ix)*c5*data(iz,ix)
!          model(iz,ix  ) = model(iz,ix  ) - aux_m2m3F(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-1) = model(iz,ix-1) - aux_m2m3F(iz,ix)*c4*data(iz,ix)
!          model(iz,ix-2) = model(iz,ix-2) - aux_m2m3F(iz,ix)*c3*data(iz,ix)
!          model(iz,ix-3) = model(iz,ix-3) - aux_m2m3F(iz,ix)*c2*data(iz,ix)
!          model(iz,ix-4) = model(iz,ix-4) - aux_m2m3F(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
!!!!$OMP END PARALLEL DO
!   ! else
!   !   if(.not.add) data = 0.
!!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!   !   do ix=5,nx-5
!   !     do iz=5,nz-5
!   !       data(iz,ix) = data(iz,ix) + aux_m2m3(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
!   !                                                    c2*(model(iz,ix+4) - model(iz,ix-3))+&
!   !                                                    c3*(model(iz,ix+3) - model(iz,ix-2))+&
!   !                                                    c4*(model(iz,ix+2) - model(iz,ix-1))+&
!   !                                                    c5*(model(iz,ix+1) - model(iz,ix  )))
!   !     enddo
!   !   enddo
!!!!!$OMP END PARALLEL DO
!   ! endif
!  endsubroutine
!
!  subroutine EF_c(nz, nx, model,data)!,aux_m2F)
!   ! logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data !, aux_m2F
!    !write(0,*) "Ec"
!    integer :: nz, nx
!
!    !if(adj) then
!      !if(.not.add) model = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!      do ix=5,nx-5
!        do iz=5,nz-5
!          model(iz+5,ix) = model(iz+5,ix) + aux_m2F(iz,ix)*c1*data(iz,ix)
!          model(iz+4,ix) = model(iz+4,ix) + aux_m2F(iz,ix)*c2*data(iz,ix)
!          model(iz+3,ix) = model(iz+3,ix) + aux_m2F(iz,ix)*c3*data(iz,ix)
!          model(iz+2,ix) = model(iz+2,ix) + aux_m2F(iz,ix)*c4*data(iz,ix)
!          model(iz+1,ix) = model(iz+1,ix) + aux_m2F(iz,ix)*c5*data(iz,ix)
!          model(iz  ,ix) = model(iz  ,ix) - aux_m2F(iz,ix)*c5*data(iz,ix)
!          model(iz-1,ix) = model(iz-1,ix) - aux_m2F(iz,ix)*c4*data(iz,ix)
!          model(iz-2,ix) = model(iz-2,ix) - aux_m2F(iz,ix)*c3*data(iz,ix)
!          model(iz-3,ix) = model(iz-3,ix) - aux_m2F(iz,ix)*c2*data(iz,ix)
!          model(iz-4,ix) = model(iz-4,ix) - aux_m2F(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
!!!!$OMP END PARALLEL DO
!  !  else
!  !    if(.not.add) data = 0.
!!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!  !    do ix=5,nx-5
!  !      do iz=5,nz-5
!  !        data(iz,ix) = data(iz,ix) + aux_m2(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
!  !                                                   c2*(model(iz+4,ix) - model(iz-3,ix))+&
!  !                                                   c3*(model(iz+3,ix) - model(iz-2,ix))+&
!  !                                                   c4*(model(iz+2,ix) - model(iz-1,ix))+&
!  !                                                   c5*(model(iz+1,ix) - model(iz  ,ix)))
!  !      enddo
!  !    enddo
!!!!!$OMP END PARALLEL DO
!  !  endif
!  endsubroutine
!
!  subroutine FF_c(nz, nx, model,data)!, aux_m2F)
!    !logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data!, aux_m2F
!    !write(0,*) "Fc"
!    integer :: nz, nx
!
!   ! if(adj) then
!      !if(.not.add) model = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!      do ix=5,nx-5
!        do iz=5,nz-5
!          model(iz,ix+5) = model(iz,ix+5) + aux_m2F(iz,ix)*c1*data(iz,ix)
!          model(iz,ix+4) = model(iz,ix+4) + aux_m2F(iz,ix)*c2*data(iz,ix)
!          model(iz,ix+3) = model(iz,ix+3) + aux_m2F(iz,ix)*c3*data(iz,ix)
!          model(iz,ix+2) = model(iz,ix+2) + aux_m2F(iz,ix)*c4*data(iz,ix)
!          model(iz,ix+1) = model(iz,ix+1) + aux_m2F(iz,ix)*c5*data(iz,ix)
!          model(iz,ix  ) = model(iz,ix  ) - aux_m2F(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-1) = model(iz,ix-1) - aux_m2F(iz,ix)*c4*data(iz,ix)
!          model(iz,ix-2) = model(iz,ix-2) - aux_m2F(iz,ix)*c3*data(iz,ix)
!          model(iz,ix-3) = model(iz,ix-3) - aux_m2F(iz,ix)*c2*data(iz,ix)
!          model(iz,ix-4) = model(iz,ix-4) - aux_m2F(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
!!!!$OMP END PARALLEL DO
!   ! else
!   !   if(.not.add) data = 0.
!!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!   !   do ix=5,nx-5
!   !     do iz=5,nz-5
!   !       data(iz,ix) = data(iz,ix) + aux_m2(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
!   !                                                  c2*(model(iz,ix+4) - model(iz,ix-3))+&
!   !                                                  c3*(model(iz,ix+3) - model(iz,ix-2))+&
!   !                                                  c4*(model(iz,ix+2) - model(iz,ix-1))+&
!   !                                                  c5*(model(iz,ix+1) - model(iz,ix  )))
!   !     enddo
!   !   enddo
!!!!!$OMP END PARALLEL DO
!   ! endif
!  endsubroutine
!
!  subroutine GF_c(nz, nx, model,data)!, aux_m2m3F)
!    !logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data!, aux_m2m3F
!    !write(0,*) "Gc"
!    integer :: nz, nx
!
!    !if(adj) then
!     ! if(.not.add) model = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!      do ix=5,nx-5
!        do iz=5,nz-5
!          model(iz+5,ix) = model(iz+5,ix) + aux_m2m3F(iz,ix)*c1*data(iz,ix)
!          model(iz+4,ix) = model(iz+4,ix) + aux_m2m3F(iz,ix)*c2*data(iz,ix)
!          model(iz+3,ix) = model(iz+3,ix) + aux_m2m3F(iz,ix)*c3*data(iz,ix)
!          model(iz+2,ix) = model(iz+2,ix) + aux_m2m3F(iz,ix)*c4*data(iz,ix)
!          model(iz+1,ix) = model(iz+1,ix) + aux_m2m3F(iz,ix)*c5*data(iz,ix)
!          model(iz  ,ix) = model(iz  ,ix) - aux_m2m3F(iz,ix)*c5*data(iz,ix)
!          model(iz-1,ix) = model(iz-1,ix) - aux_m2m3F(iz,ix)*c4*data(iz,ix)
!          model(iz-2,ix) = model(iz-2,ix) - aux_m2m3F(iz,ix)*c3*data(iz,ix)
!          model(iz-3,ix) = model(iz-3,ix) - aux_m2m3F(iz,ix)*c2*data(iz,ix)
!          model(iz-4,ix) = model(iz-4,ix) - aux_m2m3F(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
!!!!$OMP END PARALLEL DO
!   ! else
!   !   if(.not.add) data = 0.
!!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!   !   do ix=5,nx-5
!   !     do iz=5,nz-5
!   !       data(iz,ix) = data(iz,ix) + aux_m2m3(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
!   !                                                    c2*(model(iz+4,ix) - model(iz-3,ix))+&
!   !                                                    c3*(model(iz+3,ix) - model(iz-2,ix))+&
!   !                                                    c4*(model(iz+2,ix) - model(iz-1,ix))+&
!   !                                                    c5*(model(iz+1,ix) - model(iz  ,ix)))
!   !     enddo
!   !   enddo
!!!!!$OMP END PARALLEL DO
!   ! endif 
!  endsubroutine
!
!  subroutine HF_c(nz, nx, model,data)!, aux_m3F)
!    !logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data!, aux_m3F
!    !write(0,*) "Hc"
!    integer :: nz, nx
!
!    !if(adj) then
!     ! if(.not.add) model = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!      do ix=6,nx-5
!        do iz=6,nz-5
!          model(iz+4,ix) = model(iz+4,ix) + aux_m3F(iz,ix)*c1*data(iz,ix)
!          model(iz+3,ix) = model(iz+3,ix) + aux_m3F(iz,ix)*c2*data(iz,ix)
!          model(iz+2,ix) = model(iz+2,ix) + aux_m3F(iz,ix)*c3*data(iz,ix)
!          model(iz+1,ix) = model(iz+1,ix) + aux_m3F(iz,ix)*c4*data(iz,ix)
!          model(iz  ,ix) = model(iz  ,ix) + aux_m3F(iz,ix)*c5*data(iz,ix)
!          model(iz-1,ix) = model(iz-1,ix) - aux_m3F(iz,ix)*c5*data(iz,ix)
!          model(iz-2,ix) = model(iz-2,ix) - aux_m3F(iz,ix)*c4*data(iz,ix)
!          model(iz-3,ix) = model(iz-3,ix) - aux_m3F(iz,ix)*c3*data(iz,ix)
!          model(iz-4,ix) = model(iz-4,ix) - aux_m3F(iz,ix)*c2*data(iz,ix)
!          model(iz-5,ix) = model(iz-5,ix) - aux_m3F(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
!!!!$OMP END PARALLEL DO
!   ! else
!   !   if(.not.add) data = 0.
!!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!   !   do ix=6,nx-5
!   !     do iz=6,nz-5
!   !       data(iz,ix) = data(iz,ix) + aux_m3F(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
!   !                                                  c2*(model(iz+3,ix) - model(iz-4,ix))+&
!   !                                                  c3*(model(iz+2,ix) - model(iz-3,ix))+&
!   !                                                  c4*(model(iz+1,ix) - model(iz-2,ix))+&
!   !                                                  c5*(model(iz  ,ix) - model(iz-1,ix)))
!   !     enddo
!   !   enddo
!!!!!$OMP END PARALLEL DO
!   ! endif
!  endsubroutine
!
!  subroutine JF_c(nz, nx, model,data)!, aux_m3F)
!    !logical,intent(in)          :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: model,data!, aux_m3F
!    !write(0,*) "Hc"
!    integer :: nz, nx
!
!    !if(adj)then
!     ! if(.not.add) model = 0.
!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!      do ix=6,nx-5
!        do iz=6,nz-5
!          model(iz,ix+4) = model(iz,ix+4) + aux_m3F(iz,ix)*c1*data(iz,ix)
!          model(iz,ix+3) = model(iz,ix+3) + aux_m3F(iz,ix)*c2*data(iz,ix)
!          model(iz,ix+2) = model(iz,ix+2) + aux_m3F(iz,ix)*c3*data(iz,ix)
!          model(iz,ix+1) = model(iz,ix+1) + aux_m3F(iz,ix)*c4*data(iz,ix)
!          model(iz,ix  ) = model(iz,ix  ) + aux_m3F(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-1) = model(iz,ix-1) - aux_m3F(iz,ix)*c5*data(iz,ix)
!          model(iz,ix-2) = model(iz,ix-2) - aux_m3F(iz,ix)*c4*data(iz,ix)
!          model(iz,ix-3) = model(iz,ix-3) - aux_m3F(iz,ix)*c3*data(iz,ix)
!          model(iz,ix-4) = model(iz,ix-4) - aux_m3F(iz,ix)*c2*data(iz,ix)
!          model(iz,ix-5) = model(iz,ix-5) - aux_m3F(iz,ix)*c1*data(iz,ix)
!        enddo
!      enddo
!!!!$OMP END PARALLEL DO
!   ! else
!   !   if(.not.add) data = 0.
!!!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
!   !   do ix=6,nx-5
!   !     do iz=6,nz-5
!   !     data(iz,ix) = data(iz,ix) + aux_m3F(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
!   !                                                c2*(model(iz,ix+3) - model(iz,ix-4))+&
!   !                                                c3*(model(iz,ix+2) - model(iz,ix-3))+&
!   !                                                c4*(model(iz,ix+1) - model(iz,ix-2))+&
!   !                                                c5*(model(iz,ix  ) - model(iz,ix-1)))
!   !     enddo
!   !   enddo
!!!!!$OMP END PARALLEL DO
!   ! endif
!  endsubroutine
!
!
  end module
