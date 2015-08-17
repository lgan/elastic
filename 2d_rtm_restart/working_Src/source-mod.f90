! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! source module
module source_mod

  use sep

  implicit none
  real, private      :: fp, dt, dx
  integer, private   :: nsource, nx, nz

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine source_init(fp_in,dt_in, dx_in, nsource_in, nz_in, nx_in)
    real                :: fp_in, dt_in, dx_in
    integer             :: nsource_in, nx_in, nz_in

    fp = fp_in
    dt = dt_in
    dx = dx_in
    nsource = nsource_in
    nz = nz_in
    nx = nx_in
  end subroutine

  ! Function that is called by elastic_mod
  function source_op(it, source) result(stat)
    real,intent(inout)       :: source
    !real,intent(inout)       :: source
    integer,intent(in)          :: it
    integer                     :: stat

    call source_op2(it, source)
    stat=0
  end function

  ! Ricker source
  subroutine source_op2(it, source)
    real                :: pi, aux, t
    integer             :: it
    !real             :: source
    real             :: source

    source=0.0
    pi = 3.14159265
    t = (it*dt)-(nsource*dt/2.0)
    aux = (pi*pi*fp*fp*t*t)
    source = (1.0 - 2.0*aux)*exp(-aux)
    if (it.le.(nsource/10)) source=0

    source = 10.0e7*source

  endsubroutine

  ! Ricker source
  subroutine source_op3(it, source)
    real                :: pi, aux, t
    integer             :: it
    !real             :: source
    real             :: source

    source=0.0
    pi = 3.14159265
    t = (it*dt)-(nsource*dt/2.0)
    aux = (pi*pi*1.3*fp*fp*t*t)
    source = (1.0 - 2.0*aux)*exp(-aux)
    if (it.le.(nsource/10)) source=0

    source = -10.0e7*source

  endsubroutine

  ! Distribution of point sources to create an explosive source
  subroutine explosive_sourcex(source,xsource,zsource,data)!,sourcezz)
    !logical                  :: adj,add
    real                     :: source
    integer                  :: xsource,zsource
    real,dimension(nz,nx)    :: data

    data(zsource  ,xsource) = data(zsource  ,xsource) + source
    !data(zsource  ,xsource-1) = data(zsource  ,xsource-1) + source
    !data(zsource+1,xsource  ) = data(zsource+1,xsource  ) + source
    !data(zsource,xsource-1) = data(zsource,xsource-1) + source
  endsubroutine

  subroutine explosive_sourcez(source,xsource,zsource,data)!,sourcezz)
    !logical                  :: adj,add
    real                     :: source
    integer                  :: xsource,zsource
    real,dimension(nz,nx)    :: data

    data(zsource  ,xsource  ) = data(zsource  ,xsource  ) + source
    !data(zsource  ,xsource  ) = data(zsource  ,xsource  ) + source
    !data(zsource+1,xsource-1) = data(zsource+1,xsource-1) + source
    !data(zsource-1,xsource) = data(zsource-1,xsource) + source
  endsubroutine

endmodule
