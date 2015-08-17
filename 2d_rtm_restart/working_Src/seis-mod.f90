! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! seismogram module
module seis_mod

  use sep

  implicit none
  integer, private :: nz,nx,abc,zrec
  real, dimension(:,:), pointer, private :: seis_Vx, seis_Vz, seis_P

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine seis_init(nz_in,nx_in, abc_in,zrec_in, seis_Vx_in, seis_Vz_in, seis_P_in)
    integer,intent(in)  :: nz_in,nx_in, abc_in, zrec_in
    real, dimension(:,:), target :: seis_Vx_in, seis_Vz_in, seis_P_in

    nz = nz_in
    nx = nx_in
    abc = abc_in
    zrec = zrec_in
    seis_Vx => seis_Vx_in
    seis_Vz => seis_Vz_in
    seis_P  => seis_P_in

    if(abc.eq.0) abc=1

  end subroutine

  ! Function that is called by elastic_mod
  function seis_op(countseis,Vx,Vz) result(stat)
    integer,intent(in)          :: countseis
    !real*16,dimension(:)        :: Vx, Vz
    real,dimension(:)        :: Vx, Vz
    integer                     :: stat
    call seis_op2(countseis,Vx,Vz)
    stat=0
  end function

  subroutine seis_op2(countseis,Vx,Vz)
    integer,intent(in)          :: countseis
    !real*16,dimension(nz,nx)    :: Vx, Vz
    real,dimension(nz,nx)       :: Vx, Vz
    integer                     :: ix, iz
    !real*16                     :: dxVx,dxVz,dzVx,dzVz,aux

    do ix=1,nx-2*abc-1
      seis_Vx(countseis,ix) = 0.5*(Vx(zrec,ix+abc)+Vx(zrec-1,ix+abc))
      seis_Vz(countseis,ix) = 0.5*(Vz(zrec,ix+abc)+Vz(zrec,ix+abc-1))
    enddo

  endsubroutine

  subroutine seis_op3(countseis,sigmaxx,sigmazz)
    integer,intent(in)          :: countseis
    !real*16,dimension(nz,nx)    :: Vx, Vz
    real,dimension(nz,nx)       :: sigmaxx, sigmazz
    integer                     :: ix, iz
    !real*16                     :: dxVx,dxVz,dzVx,dzVz,aux

    do ix=1,nx-2*abc-1
      seis_P(countseis,ix) = 0.5*(sigmaxx(zrec,ix+abc)+sigmazz(zrec-1,ix+abc))
    enddo

  endsubroutine

endmodule

