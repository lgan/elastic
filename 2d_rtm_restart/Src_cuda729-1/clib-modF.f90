! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Jan, 5th, 2015

! finite difference module
module clib_mod

  use, intrinsic :: iso_c_binding

  implicit none

  interface
    subroutine dsigmaxxdx_c(nx,nz,ix,model,dataF,rho_Vx) bind(c,name='dsigmaxxdx_c_')
      import
      integer(c_int),value        :: nz,nx,ix
      real(c_float),dimension(*)  :: model,rho_Vx
      real(c_float),dimension(*)  :: dataF
    endsubroutine

    subroutine dsigmaxzdz_c(nx,nz,ix,model,dataF,rho_Vx) bind(c,name='dsigmaxzdz_c_')
      import
      integer(c_int),value        :: nz,nx,ix
      real(c_float),dimension(*)  :: model,rho_Vx
      real(c_float),dimension(*)  :: dataF
    endsubroutine

    subroutine dsigmaxzdx_c(nx,nz,ix,model,dataF,rho_Vz) bind(c,name='dsigmaxzdx_c_')
      import
      integer(c_int),value        :: nz,nx,ix
      real(c_float),dimension(*)  :: model,rho_Vz
      real(c_float),dimension(*)  :: dataF
    endsubroutine

    subroutine dsigmazzdz_c(nx,nz,ix,model,dataF,rho_Vz) bind(c,name='dsigmazzdz_c_')
      import
      integer(c_int),value        :: nz,nx,ix
      real(c_float),dimension(*)  :: model,rho_Vz
      real(c_float),dimension(*)  :: dataF
    endsubroutine

    subroutine dvzdx_c(nx,nz,ix,model,dataF) bind(c,name='dvzdx_c_')
      import
      integer(c_int),value        :: nz,nx,ix
      real(c_float),dimension(*)  :: model
      real(c_float),dimension(*)  :: dataF
    endsubroutine

    subroutine dvxdz_c(nx,nz,ix,model,dataF) bind(c,name='dvxdz_c_')
      import
      integer(c_int),value        :: nz,nx,ix
      real(c_float),dimension(*)  :: model
      real(c_float),dimension(*)  :: dataF
    endsubroutine

    subroutine dvxdx_c(nx,nz,ix,model,dataF) bind(c,name='dvxdx_c_')
      import
      integer(c_int),value        :: nz,nx,ix
      real(c_float),dimension(*)  :: model
      real(c_float),dimension(*)  :: dataF
    endsubroutine

    subroutine dvzdz_c(nx,nz,ix,model,dataF) bind(c,name='dvzdz_c_')
      import
      integer(c_int),value        :: nz,nx,ix
      real(c_float),dimension(*)  :: model
      real(c_float),dimension(*)  :: dataF
    endsubroutine

  endinterface
end module
