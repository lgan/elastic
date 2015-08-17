! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 12th, 2015

! This module creates the density and Lame
! background models in the staggered grid.

module model_mod

  !use abc_mod

  implicit none

  integer, private :: nz,nx!,abc
  real, private    :: dx,dt
  real, pointer, dimension(:,:),private :: m1,m2,m3,m1_x,m1_z
  real, pointer, dimension(:,:),private :: dm1,dm1_x,dm1_z,dm2,dm3

  contains

  subroutine model_op(nz_in,nx_in,dx_in,dt_in,m1_in,m2_in,m3_in,m1_x_in,m1_z_in)

    integer                     :: iz,ix
    real                        :: dx_in,dt_in,aux1
    integer                     :: nx_in, nz_in!, abc_in
    real, dimension(:,:),target :: m1_in,m2_in,m3_in,m1_x_in,m1_z_in

    nx = nx_in
    nz = nz_in
    dx = dx_in
    dt = dt_in

    m1 => m1_in
    m2 => m2_in
    m3 => m3_in
    m1_x => m1_x_in
    m1_z => m1_z_in

    !Calculate inverse of density
    !These also calculate the values of density at the velocity staggered grid,
    !using the arithmetic average proposed in Moczo, 2002

    !Rho_vx
    do ix=1,nx
      do iz=1,nz-1
        m1_x(iz,ix)=1.0/(0.5*(m1(iz+1,ix) + m1(iz,ix)))
      enddo
    enddo
    !Rho_vz
    do ix=1,nx-1
      do iz=1,nz
        m1_z(iz,ix)=1.0/(0.5*(m1(iz,ix+1) + m1(iz,ix)))
      enddo
    enddo
    !Right side
    do iz=1,nz
      m1_z(iz,nx)=(1.0/m1(iz,nx))
    enddo
    !Bottom
    do ix=1,nx
      m1_x(nz,ix)=(1.0/m1(nz,ix))
    enddo
    !Calculate weigthed average for the lambda variable
    !These put lambda and (lambda+2mi) in their respective staggered grids,
    !using the harmonic average proposed in Moczo, 2002
    do ix=2,nx-1
      do iz=2,nz-1
        m2(iz,ix)=4.0/((1.0/m2(iz,ix))+(1.0/m2(iz+1,ix))+(1.0/m2(iz,ix+1))+(1.0/m2(iz+1,ix+1)))
        if (m3(iz,ix)==0.or.m3(iz+1,ix)==0.or.m3(iz,ix+1)==0.or.m3(iz+1,ix+1)==0) then
          m3(iz,ix)=0.0
        else
          m3(iz,ix)=4.0/((1.0/m3(iz,ix))+(1.0/m3(iz+1,ix))+(1.0/m3(iz,ix+1))+(1.0/m3(iz+1,ix+1)))
        endif
      enddo
    enddo
    do ix=1,nx
      m2(1,ix) = m2(2,ix)
      m2(nz,ix) = m2(nz-1,ix)
      m3(1,ix) = m3(2,ix)
      m3(nz,ix) = m3(nz-1,ix)
    enddo
    do iz=1,nz
      m2(iz,1) = m2(iz,2)
      m2(iz,nx) = m2(iz,nx-1)
      m3(iz,1) = m3(iz,2)
      m3(iz,nx) = m3(iz,nx-1)
    enddo

  endsubroutine

  subroutine model_delta(nz_in,nx_in,dm1_in,dm1_x_in,dm1_z_in,dm2_in,dm3_in)
    integer                     :: iz,ix
    integer                     :: nz_in, nx_in
    real, dimension(:,:),target :: dm1_x_in,dm1_z_in
    real, dimension(:,:),target :: dm1_in,dm2_in,dm3_in

    nx = nx_in
    nz = nz_in
    dm1 => dm1_in
    dm1_x => dm1_x_in
    dm1_z => dm1_z_in
    dm2 => dm2_in
    dm3 => dm3_in

    dm1_x     =0.0
    dm1_z     =0.0

    !dm1_x
    do ix=1,nx
      do iz=1,nz-1
        if (dm1(iz+1,ix).le.1.0e-1.or.dm1(iz,ix).le.1.0e-1) then
          dm1_x(iz,ix)=0.0
        else
          dm1_x(iz,ix)=0.5*(dm1(iz+1,ix)+dm1(iz,ix))
        endif
      enddo
    enddo
    !drhoz
    do ix=1,nx-1
      do iz=1,nz
        if (dm1(iz,ix+1).le.1.0e-1.or.dm1(iz,ix).le.1.0e-1) then
          dm1_z(iz,ix)=0.0
        else
          dm1_z(iz,ix)=0.5*(dm1(iz,ix+1)+dm1(iz,ix))
        endif
      enddo
    enddo
    !Right side
    do iz=1,nz
      if (dm1(iz,nx).le.1.0e-1) then
        dm1_z(iz,nx)=0.0
      else
        dm1_z(iz,nx)=dm1(iz,nx)
      endif
    enddo
    !Bottom
    do ix=1,nx
      if (dm1(nz,ix).le.1.0e-1) then
        dm1_x(nz,ix)=0.0
      else
        dm1_x(nz,ix)=dm1(nz,ix)
      endif
    enddo
    do ix=2,nx-1
      do iz=2,nz-1
        if (dm2(iz,ix)<1e-2.or.dm2(iz+1,ix)<1e-2.or.dm2(iz,ix+1)<1e-2.or.dm2(iz+1,ix+1)<1e-2) then
          dm2(iz,ix)=0.0
        else
          dm2(iz,ix)=4.0/((1.0/dm2(iz,ix))+(1.0/dm2(iz+1,ix))+(1.0/dm2(iz,ix+1))+(1.0/dm2(iz+1,ix+1)))
        endif
        if (dm3(iz,ix)<1e-2.or.dm3(iz+1,ix)<1e-2.or.dm3(iz,ix+1)<1e-2.or.dm3(iz+1,ix+1)<1e-2) then
          dm3(iz,ix)=0.0
        else
          dm3(iz,ix)=4.0/((1.0/dm3(iz,ix))+(1.0/dm3(iz+1,ix))+(1.0/dm3(iz,ix+1))+(1.0/dm3(iz+1,ix+1)))
        endif
      enddo
    enddo
    do ix=1,nx
      dm2(1,ix) = dm2(2,ix)
      dm2(nz,ix) = dm2(nz-1,ix)
      dm3(1,ix) = dm3(2,ix)
      dm3(nz,ix) = dm3(nz-1,ix)
    enddo
    do iz=1,nz
      dm2(iz,1) = dm2(iz,2)
      dm2(iz,nx) = dm2(iz,nx-1)
      dm3(iz,1) = dm3(iz,2)
      dm3(iz,nx) = dm3(iz,nx-1)
    enddo

  endsubroutine

end module
