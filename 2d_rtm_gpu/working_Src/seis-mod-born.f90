! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! seismogram module
module seis_mod_b

  use sep

  implicit none
  integer, private :: nz,nx,abc,zrec
  real, dimension(:,:), pointer, private :: seis_Vx, seis_Vz, seis_P
  real, dimension(:,:), pointer, private :: seis_Vx_b, seis_Vz_b, seis_P_b

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine seis_init(nz_in,nx_in, abc_in,zrec_in, seis_Vx_in, seis_Vz_in, seis_P_in, &
                                              seis_Vx_b_in, seis_Vz_b_in, seis_P_b_in)
    integer,intent(in)  :: nz_in,nx_in, abc_in, zrec_in
    real, dimension(:,:), target :: seis_Vx_in, seis_Vz_in, seis_P_in
    real, dimension(:,:), target :: seis_Vx_b_in, seis_Vz_b_in, seis_P_b_in

    nz = nz_in
    nx = nx_in
    abc = abc_in
    zrec = zrec_in
    seis_Vx => seis_Vx_in
    seis_Vz => seis_Vz_in
    seis_P  => seis_P_in
    seis_Vx_b => seis_Vx_b_in
    seis_Vz_b => seis_Vz_b_in
    seis_P_b  => seis_P_b_in

    if(abc.eq.0) abc=1

  end subroutine

  subroutine seis_op2(countseis,Vx_x,Vx_z,Vz_x,Vz_z)
    integer,intent(in)          :: countseis
    !real*16,dimension(nz,nx)    :: Vx, Vz
    real,dimension(nz,nx)       :: Vx_x,Vx_z,Vz_x,Vz_z
    integer                     :: ix, iz
    !real*16                     :: dxVx,dxVz,dzVx,dzVz,aux

    do ix=1,nx-2*abc-1
      seis_Vx(countseis,ix) = 0.5*(Vx_x(zrec,ix+abc)+Vx_z(zrec,ix+abc)+&
                                   Vx_x(zrec-1,ix+abc)+Vx_z(zrec-1,ix+abc))
      seis_Vz(countseis,ix) = 0.5*(Vz_x(zrec,ix+abc)+Vz_z(zrec,ix+abc)+&
                                   Vz_x(zrec,ix+abc-1)+Vz_z(zrec,ix+abc-1))
    enddo

  endsubroutine

  subroutine seis_op3(countseis,sigmaxx_x,sigmaxx_z,sigmazz_x,sigmazz_z)
    integer,intent(in)          :: countseis
    !real*16,dimension(nz,nx)    :: Vx, Vz
    real,dimension(nz,nx)       :: sigmaxx_x,sigmaxx_z,sigmazz_x,sigmazz_z
    integer                     :: ix, iz
    !real*16                     :: dxVx,dxVz,dzVx,dzVz,aux

    do ix=1,nx-2*abc-1
      seis_P(countseis,ix) = 0.5*(sigmaxx_x(zrec,ix+abc)+sigmaxx_z(zrec,ix+abc)+&
                                  sigmazz_x(zrec-1,ix+abc)+sigmazz_z(zrec-1,ix+abc))
    enddo

  endsubroutine

  subroutine seis_op4(countseis,Vx_x,Vx_z,Vz_x,Vz_z)
    integer,intent(in)          :: countseis
    !real*16,dimension(nz,nx)    :: Vx, Vz
    real,dimension(nz,nx)       :: Vx_x,Vx_z,Vz_x,Vz_z
    integer                     :: ix, iz
    !real*16                     :: dxVx,dxVz,dzVx,dzVz,aux

    do ix=1,nx-2*abc-1
      seis_Vx_b(countseis,ix) = 0.5*(Vx_x(zrec,ix+abc)+Vx_z(zrec,ix+abc)+&
                                     Vx_x(zrec-1,ix+abc)+Vx_z(zrec-1,ix+abc))
      seis_Vz_b(countseis,ix) = 0.5*(Vz_x(zrec,ix+abc)+Vz_z(zrec,ix+abc)+&
                                     Vz_x(zrec,ix+abc-1)+Vz_z(zrec,ix+abc-1))
    enddo

  endsubroutine

  subroutine seis_op5(countseis,sigmaxx_x,sigmaxx_z,sigmazz_x,sigmazz_z)
    integer,intent(in)          :: countseis
    !real*16,dimension(nz,nx)    :: Vx, Vz
    real,dimension(nz,nx)       :: sigmaxx_x,sigmaxx_z,sigmazz_x,sigmazz_z
    integer                     :: ix, iz
    !real*16                     :: dxVx,dxVz,dzVx,dzVz,aux

    do ix=1,nx-2*abc-1
      seis_P_b(countseis,ix) = 0.5*(sigmaxx_x(zrec,ix+abc)+sigmaxx_z(zrec,ix+abc)+&
                                    sigmazz_x(zrec-1,ix+abc)+sigmazz_z(zrec-1,ix+abc))
    enddo

  endsubroutine

endmodule

