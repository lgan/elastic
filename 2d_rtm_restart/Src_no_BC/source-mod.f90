! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! source module
module source_mod

  use sep

  implicit none
  real, private      :: fp, dt, dx
  integer, private   :: nsource, nx, nz, it, nt

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine source_init(fp_in,dt_in, dx_in, nsource_in, nz_in, nx_in, nt_in)
    real                :: fp_in, dt_in, dx_in
    integer             :: nsource_in, nx_in, nz_in, nt_in

    fp = fp_in
    dt = dt_in
    dx = dx_in
    nsource = nsource_in
    nz = nz_in
    nx = nx_in
    nt = nt_in
  end subroutine

  ! Ricker source
  subroutine source_op2(source)
    real                :: pi, aux, t
    integer             :: it,ft
    real,dimension(nt)  :: source

    source=0.0
    ft = nsource/10 !first nonzero sample
    pi = 3.14159265
    do it=ft,nt
      t = (it*dt)-(nsource*dt/2.0)
      aux = (pi*pi*fp*fp*t*t)
      source(it) = 10.0e7*(1.0 - 2.0*aux)*exp(-aux)
    enddo
  endsubroutine

  ! Distribution of point sources to create an explosive source
  subroutine explosive_source(source,xsource,zsource,data)
    real                     :: source
    integer                  :: xsource,zsource
    real,dimension(nz,nx)    :: data

    data(zsource,xsource) = data(zsource,xsource) + dt*source
  endsubroutine

endmodule
