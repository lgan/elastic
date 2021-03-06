! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! This main is to be used with the velocity-stress formulation
! Jan, 12th, 2015

! Main program to execute elastic modeling
program elastic_program

  use sep
  use rtm_mod

  use omp_lib

  implicit none

  integer                               :: nz,nx
  integer                               :: stat
  real                                  :: dz,oz,dx,ox,dt
  real, allocatable, dimension(:,:)     :: rho0,lambda0,mi0
  real                                  :: vmax

  !Variables to measure runtime
  real*8                                :: start_time,end_time

  !OMP settings
  integer                               :: node, nodes
  start_time = omp_get_wtime()

  call sep_init()

  ! Read SEP header parameter
  ! Model parameters
  call from_history("n1",nz)
  call from_history("o1",oz)
  call from_history("d1",dz)
  call from_history("n2",nx)
  call from_history("o2",ox)
  call from_history("d2",dx)

  !OMP
  nodes=16

  node = omp_get_num_procs()
  call omp_set_num_threads(nodes)

  ! dt from the Courant stability condition => C=(vx*dt)/dx+(vz*dt)dz<1
  ! For the simple case of dz=dx and using Virieux(1986) relation:
  dt = 0.4*dx/(sqrt(2.0)*2000.0) !the 0.5 is arbitrary, but it was diverging otherwise

  allocate(lambda0(nz,nx))
  allocate(mi0(nz,nx))
  allocate(rho0(nz,nx))

  ! Initialize values
  lambda0=0.0
  mi0=0.0
  rho0=0.0

  call sep_read(lambda0)
  call sep_read(mi0,"mi0")
  call sep_read(rho0,"rho0")

  ! Write outputs to the terminal
  write(0,*)"============================================================="
  write(0,*)"BEGIN OF SUBROUTINE CALLS"
  write(0,*)"============================================================="
  write(0,*)"Values of the stability conditions:"
  write(0,*)"dx=dy=",dx
  write(0,*)"dt=",dt
  write(0,*)"OMP using",nodes,"threads out of",node
  write(0,*)"============================================================="

  ! Initialize the elastic modeling
  call elastic_init(nz,dz,nx,dx,dt)


  stat = elastic_op(lambda0,mi0,rho0)
  if (stat) then
    write(0,*) "ERROR ON ELASTIC_OP FUNCTION CALL ==> EXITING"
    stop
  endif

  end_time = omp_get_wtime()

  ! Write outputs to the terminal
  write(0,*)"============================================================="
  write(0,*)"TOTAL ELAPSED TIME"
  write(0,*)"============================================================="
  write(0,*)"Total=",end_time - start_time, "seconds."
  write(0,*)"============================================================="

end program
