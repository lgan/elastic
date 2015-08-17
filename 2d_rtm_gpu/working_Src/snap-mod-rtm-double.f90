! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! snapshot module
module snap_mod_double

  implicit none
  integer, private :: nx, nz, abc
  real, dimension(:,:,:), pointer, private :: snap_Vx, snap_Vz, snap_P
  real, dimension(:,:,:), pointer, private :: snap_Vx_b, snap_Vz_b, snap_P_b

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine snap_init(nz_in, nx_in, abc_in, snap_Vx_in, snap_Vz_in, snap_P_in,&
                       snap_Vx_b_in, snap_Vz_b_in, snap_P_b_in)
    integer             :: nx_in, nz_in, abc_in
    real, dimension(:,:,:), target :: snap_Vx_in, snap_Vz_in, snap_P_in
    real, dimension(:,:,:), target :: snap_Vx_b_in, snap_Vz_b_in, snap_P_b_in
    nx = nx_in
    nz = nz_in
    abc= abc_in
    snap_Vx => snap_Vx_in
    snap_Vz => snap_Vz_in
    snap_P => snap_P_in
    snap_Vx_b => snap_Vx_b_in
    snap_Vz_b => snap_Vz_b_in
    snap_P_b => snap_P_b_in
  end subroutine

  subroutine snap_op2(countsnap,Vx_x,Vx_z,Vz_x,Vz_z)
    integer,intent(in)          :: countsnap
    real*8,dimension(nz,nx)       :: Vx_x,Vx_z,Vz_x,Vz_z
    integer                     :: ix,iz

    do ix=2,nx-2
      do iz=2,nz-2
        snap_Vx(iz,ix,countsnap)= Vx_x(iz,ix)+Vx_z(iz,ix)
        snap_Vz(iz,ix,countsnap)= Vz_x(iz,ix)+Vz_z(iz,ix)
      enddo
    enddo

  endsubroutine

  subroutine snap_op3(countsnap,sigmaxx_x,sigmaxx_z,sigmazz_x,sigmazz_z)
    integer,intent(in)          :: countsnap
    real*8,dimension(nz,nx)       :: sigmaxx_x,sigmaxx_z,sigmazz_x,sigmazz_z
    integer                     :: ix,iz

    do ix=2,nx-2
      do iz=2,nz-2
        snap_P(iz,ix,countsnap) = 0.5*(sigmaxx_x(iz,ix)+sigmaxx_z(iz,ix)+sigmazz_x(iz,ix)+sigmazz_z(iz,ix))
      enddo
    enddo

  endsubroutine

  subroutine snap_op4(countsnap,Vx_x_b,Vx_z_b,Vz_x_b,Vz_z_b)
    integer,intent(in)          :: countsnap
    real*8,dimension(nz,nx)       :: Vx_x_b,Vx_z_b,Vz_x_b,Vz_z_b
    integer                     :: ix,iz

    do ix=2,nx-2
      do iz=2,nz-2
        snap_Vx_b(iz,ix,countsnap)= Vx_x_b(iz,ix)+Vx_z_b(iz,ix)
        snap_Vz_b(iz,ix,countsnap)= Vz_x_b(iz,ix)+Vz_z_b(iz,ix)
      enddo
    enddo

  endsubroutine

  subroutine snap_op5(countsnap,sigmaxx_x_b,sigmaxx_z_b,sigmazz_x_b,sigmazz_z_b)
    integer,intent(in)          :: countsnap
    real*8,dimension(nz,nx)       :: sigmaxx_x_b,sigmaxx_z_b,sigmazz_x_b,sigmazz_z_b
    integer                     :: ix,iz

    do ix=2,nx-2
      do iz=2,nz-2
        snap_P_b(iz,ix,countsnap) = 0.5*(sigmaxx_x_b(iz,ix)+sigmaxx_z_b(iz,ix)+sigmazz_x_b(iz,ix)+sigmazz_z_b(iz,ix))
      enddo
    enddo

  endsubroutine


endmodule
