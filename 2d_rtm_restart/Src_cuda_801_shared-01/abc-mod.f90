! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! November, 10th, 2014

! boundary module
module abc_mod

  use sep

  implicit none
  integer, private   :: nx, nz, iz, ix, abc, abc2
  real, private      :: dt, dx
  real, private, allocatable, dimension(:) :: filter
!  real, private, pointer, dimension(:,:) :: l2mi, rho 
!  logical, private :: free_surface

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine abc_init(nz_in,nx_in, dx_in, dt_in,abc_in)!,l2mi_in,rho_in,free_surface_in)
    integer        :: nx_in,nz_in,abc_in
    real           :: dx_in,dt_in,alfa,aux
!    logical        :: free_surface_in
!    real, dimension(:,:), target :: l2mi_in, rho_in

    nx = nx_in
    nz = nz_in
    dx = dx_in
    dt = dt_in
!    free_surface = free_surface_in
    abc= abc_in
!    l2mi => l2mi_in
!    rho => rho_in

    abc2 = abc+5
    alfa = 0.09
    aux = 0.0

    allocate(filter(abc2))

!    write(0,*)"Max Vp=",maxval(sqrt(l2mi/rho))

    filter=0.0
    do ix=1,abc
      aux=(1.0*ix)/(1.0*abc)
      filter(ix)=exp(-alfa*aux*aux)
    enddo
    do ix=abc+1,abc2
      filter(ix)=filter(abc)
    enddo

    call to_history("n1",abc2, "Dat/filter.H")
    call to_history("o1",0.0, "Dat/filter.H")
    call to_history("d1",1, "Dat/filter.H")
    call sep_write(filter,"Dat/filter.H")

  endsubroutine

  subroutine K(data)
    real, dimension(nz,nx) :: data

    do ix=1,nx
      do iz=1,abc2
        !TOP
        data(iz,ix) = data(iz,ix)*filter(abc2-iz+1)
        !BOTTOM
        data(nz-abc2+iz,ix) = data(nz-abc2+iz,ix)*filter(iz)
      enddo
    enddo
    do ix=1,abc2
      do iz=1,nz
        !LEFT
        data(iz,ix) = data(iz,ix)*filter(abc2-ix+1)
        !RIGHT
        data(iz,nx-abc2+ix) = data(iz,nx-abc2+ix)*filter(ix)
      enddo
    enddo
  endsubroutine

endmodule
