! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! This main is to be used with the velocity-stress formulation
! Jan, 12th, 2015

! Main program to execute elastic modeling
program elastic_program

  use sep
  use elastic_mod
  use snap_mod_rtm
  use seis_mod_rtm

  use omp_lib

  implicit none

  integer                               :: nz,nx,xsource,zsource,nt,abc,zrec,nsnap
  integer                               :: stat
  real                                  :: dz,oz,dx,ox,dt,fp,Tt,ot,dtseis
  real, allocatable, dimension(:,:)     :: rho,lambda,mu
  real, allocatable, dimension(:,:)     :: Vx_rec,Vz_rec,sigmaxx_rec,sigmazz_rec,sigmaxz_rec
  real                                  :: vmax
  logical                               :: free_surface

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

  call from_param("nsnap",nsnap,50)        !# of snapshots, default=50
  call from_param("xsource",xsource)       !position of source in x samples
  call from_param("zsource",zsource)       !position of source in z samples
  call from_param("Ttime",Tt)              !total propagation time
  call from_param("ot",ot,0.0)             !start of seismogram, default=0
  call from_param("fp",fp,25.0)            !Ricker source peak freq.
  call from_param("abc",abc,0)             !# of absorbing border samples
  call from_param("zrec",zrec,1)           !receivers depth
  call from_param("dtseis",dtseis,0.002)   !seismogram sampling, default=2ms
  call from_param("free_surface",free_surface,.false.)

  !OMP
  nodes=16

  node = omp_get_num_procs()
  call omp_set_num_threads(nodes)

  ! dt from the Courant stability condition => C=(vx*dt)/dx+(vz*dt)dz<1
  ! For the simple case of dz=dx and using Virieux(1986) relation:
  dt = 0.4*dx/(sqrt(2.0)*2000.0) !the 0.4 is arbitrary, but it was diverging otherwise
  nt = Tt/dt

  allocate(lambda(nz,nx))
  allocate(mu(nz,nx))
  allocate(rho(nz,nx))
  allocate(Vx_rec(nt,nx))
  allocate(Vz_rec(nt,nx))
  allocate(sigmaxx_rec(nt,nx))
  allocate(sigmazz_rec(nt,nx))
  allocate(sigmaxz_rec(nt,nx))

  ! Initialize values
  lambda = 0.
  mu     = 0.
  rho    = 0.  
  Vx_rec = 0.
  Vz_rec = 0.
  sigmaxx_rec = 0.
  sigmazz_rec = 0.
  sigmaxz_rec = 0.

  call sep_read(lambda)
  call sep_read(mu,"mu")
  call sep_read(rho,"rho")

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
  call elastic_init(nz,dz,nx,dx,dt,xsource,zsource,zrec,fp,nt,abc,free_surface,lambda,mu,rho)

  stat = elastic_op(Vx_rec,Vz_rec,sigmaxx_rec,sigmazz_rec,sigmaxz_rec)

  if (stat) then
    write(0,*) "ERROR ON ELASTIC_OP FUNCTION CALL ==> EXITING"
    stop
  endif

  call seismo(nx,dx,nt,dt,Vx_rec,Vz_rec,sigmaxx_rec,sigmazz_rec,sigmaxz_rec)

  call sep_close()

  end_time = omp_get_wtime()

  ! Write outputs to the terminal
  write(0,*)"============================================================="
  write(0,*)"TOTAL ELAPSED TIME"
  write(0,*)"============================================================="
  write(0,*)"Total=",end_time - start_time, "seconds."
  write(0,*)"============================================================="

end program
