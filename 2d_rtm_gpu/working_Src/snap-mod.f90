! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! snapshot module
module snap_mod

  implicit none
  integer, private :: nx, nz, abc
  real, dimension(:,:,:), pointer, private :: snap_Vx, snap_Vz, snap_P

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine snap_init(nz_in, nx_in, abc_in, snap_Vx_in, snap_Vz_in, snap_P_in)
    integer             :: nx_in, nz_in, abc_in
    real, dimension(:,:,:), target :: snap_Vx_in, snap_Vz_in, snap_P_in
    nx = nx_in
    nz = nz_in
    abc= abc_in
    snap_Vx => snap_Vx_in
    snap_Vz => snap_Vz_in
    snap_P => snap_P_in
  end subroutine

  subroutine snap_op2(countsnap,Vx,Vz)
    integer,intent(in)          :: countsnap
    !real*16,dimension(nz,nx)    :: Vx,Vz
    real,dimension(nz,nx)       :: Vx,Vz
    integer                     :: ix,iz
    !real*16                     :: aux,dxVx,dzVz,dxVz,dzVx

    do ix=2,nx-2
      do iz=2,nz-2
        snap_Vx(iz,ix,countsnap)= Vx(iz,ix)
        snap_Vz(iz,ix,countsnap)= Vz(iz,ix)
      enddo
    enddo

  endsubroutine

  subroutine snap_op3(countsnap,sigmaxx,sigmazz)
    integer,intent(in)          :: countsnap
    !real*16,dimension(nz,nx)    :: Vx,Vz
    real,dimension(nz,nx)       :: sigmaxx,sigmazz
    integer                     :: ix,iz
    !real*16                     :: aux,dxVx,dzVz,dxVz,dzVx

    do ix=2,nx-2
      do iz=2,nz-2
        snap_P(iz,ix,countsnap) = 0.5*(sigmaxx(iz,ix)+sigmazz(iz,ix))
      enddo
    enddo

  endsubroutine



endmodule
