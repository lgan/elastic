! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! This main is to be used with the velocity-stress formulation
! Jan, 12th, 2015

! Main program to execute elastic modeling
program elastic_program

  use sep
  use born_mod

  use omp_lib

  implicit none

  integer                               :: nz,nx,nt,nsource,xsource,zsource,iz,ix
  integer                               :: nsnap, stat, abc, zrec, nseis
  real                                  :: dz,oz,dx,ox,dt,ot,fp, dtseis, Tt
  real, allocatable, dimension(:,:)     :: vp0, vs0, rho0, seis_Vx, seis_Vz, seis_P,lambda0,mi0
  real, allocatable, dimension(:,:)     :: dvp, dvs, drho, dmi, dlambda
  real, allocatable, dimension(:,:,:)   :: snap_Vx,snap_Vz, snap_P
  logical                               :: lame, free_surface
  real                                  :: vmax

  real, allocatable, dimension(:,:)     :: seis_Vx_b,seis_Vz_b,seis_P_b
  real, allocatable, dimension(:,:,:)   :: snap_Vx_b,snap_Vz_b,snap_P_b

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
  !call from_param("nz",nz,1000)
  !call from_param("oz",oz,0.0)
  !call from_param("dz",dz,6.0)
  !call from_param("nx",nx,1000)
  !call from_param("ox",ox,0.0)
  !call from_param("dx",dx,6.0)
  ! Simulation parameters
  call from_param("nsnap",nsnap,50)        !# of snapshots, default=50
  call from_param("xsource",xsource)       !position of source in x samples
  call from_param("zsource",zsource)       !position of source in z samples
  call from_param("Ttime",Tt)              !total propagation time
  call from_param("ot",ot,0.0)             !start of seismogram, default=0
  call from_param("fpeak",fp,25.0)         !Ricker source peak freq.
  call from_param("abc",abc,0)             !# of absorbing border samples
  call from_param("zrec",zrec,1)           !receivers depth
  call from_param("dtseis",dtseis,0.002)   !seismogram sampling, default=2ms
  call from_param("lame",lame,.false.)     !if lame=.true., lame parameters will be
                                           !used instead of vp and vs
  call from_param("free_surface",free_surface,.false.)

  !OMP
  call from_param("nodes",nodes,15)        !number of nodes to use

  node = omp_get_num_procs()
  call omp_set_num_threads(nodes)

  ! dt from the Courant stability condition => C=(vx*dt)/dx+(vz*dt)dz<1
  ! For the simple case of dz=dx and using Virieux(1986) relation:
  dt = 0.2*dx/(sqrt(2.0)*2000.0) !the 0.5 is arbitrary, but it was diverging otherwise

  ! size of seismogram in time
  nseis=Tt/dtseis+1

  ! Write SEP header parameter
  ! seis_Vx.H output parameters
  call to_history("n1",nseis)
  call to_history("o1",ot)
  call to_history("d1",dtseis)
  call to_history("n2",(nx-2*abc))
  call to_history("o2",ox)
  call to_history("d2",dx)
  call to_history("label1","t")
  call to_history("label2","x")
  ! seis_Vz.H output parameters
  call to_history("n1",nseis,"seis_Vz")
  call to_history("o1",ot,"seis_Vz")
  call to_history("d1",dtseis,"seis_Vz")
  call to_history("n2",(nx-2*abc),"seis_Vz")
  call to_history("o2",ox,"seis_Vz")
  call to_history("d2",dx,"seis_Vz")
  call to_history("label1","t","seis_Vz")
  call to_history("label2","x","seis_Vz")
  ! seis_P.H output parameters
  call to_history("n1",nseis,"seis_P")
  call to_history("o1",ot,"seis_P")
  call to_history("d1",dtseis,"seis_P")
  call to_history("n2",(nx-2*abc),"seis_P")
  call to_history("o2",ox,"seis_P")
  call to_history("d2",dx,"seis_P")
  call to_history("label1","t","seis_P")
  call to_history("label2","x","seis_P")
  ! snapshot of x velocities output parameters
  call to_history("n1",nz,"snap_Vx")
  call to_history("o1",oz,"snap_Vx")
  call to_history("d1",dz,"snap_Vx")
  call to_history("n2",nx,"snap_Vx")
  call to_history("o2",ox,"snap_Vx")
  call to_history("d2",dx,"snap_Vx")
  call to_history("n3",nsnap,"snap_Vx")
  call to_history("o3",0,"snap_Vx")
  call to_history("d3",1,"snap_Vx")
  call to_history("label1","z","snap_Vx")
  call to_history("label2","x","snap_Vx")
  call to_history("label3","t","snap_Vx")
  ! snapshot of z velocities output parameters
  call to_history("n1",nz,"snap_Vz")
  call to_history("o1",oz,"snap_Vz")
  call to_history("d1",dz,"snap_Vz")
  call to_history("n2",nx,"snap_Vz")
  call to_history("o2",ox,"snap_Vz")
  call to_history("d2",dx,"snap_Vz")
  call to_history("n3",nsnap,"snap_Vz")
  call to_history("o3",0,"snap_Vz")
  call to_history("d3",1,"snap_Vz")
  call to_history("label1","z","snap_Vz")
  call to_history("label2","x","snap_Vz")
  call to_history("label3","t","snap_Vz")
  ! snapshot of pressure output parameters
  call to_history("n1",nz,"snap_P")
  call to_history("o1",oz,"snap_P")
  call to_history("d1",dz,"snap_P")
  call to_history("n2",nx,"snap_P")
  call to_history("o2",ox,"snap_P")
  call to_history("d2",dx,"snap_P")
  call to_history("n3",nsnap,"snap_P")
  call to_history("o3",0,"snap_P")
  call to_history("d3",1,"snap_P")
  call to_history("label1","z","snap_P")
  call to_history("label2","x","snap_P")
  call to_history("label3","t","snap_P")

  ! seis_Vx_b.H output parameters
  call to_history("n1",nseis,"seis_Vx_b")
  call to_history("o1",ot,"seis_Vx_b")
  call to_history("d1",dtseis,"seis_Vx_b")
  call to_history("n2",(nx-2*abc),"seis_Vx_b")
  call to_history("o2",ox,"seis_Vx_b")
  call to_history("d2",dx,"seis_Vx_b")
  call to_history("label1","t","seis_Vx_b")
  call to_history("label2","x","seis_Vx_b")
  ! seis_Vz_b.H output parameters
  call to_history("n1",nseis,"seis_Vz_b")
  call to_history("o1",ot,"seis_Vz_b")
  call to_history("d1",dtseis,"seis_Vz_b")
  call to_history("n2",(nx-2*abc),"seis_Vz_b")
  call to_history("o2",ox,"seis_Vz_b")
  call to_history("d2",dx,"seis_Vz_b")
  call to_history("label1","t","seis_Vz_b")
  call to_history("label2","x","seis_Vz_b")
  ! seis_P_b.H output parameters
  call to_history("n1",nseis,"seis_P_b")
  call to_history("o1",ot,"seis_P_b")
  call to_history("d1",dtseis,"seis_P_b")
  call to_history("n2",(nx-2*abc),"seis_P_b")
  call to_history("o2",ox,"seis_P_b")
  call to_history("d2",dx,"seis_P_b")
  call to_history("label1","t","seis_P_b")
  call to_history("label2","x","seis_P_b")
  ! snapshot of x velocities output parameters
  call to_history("n1",nz,"snap_Vx_b")
  call to_history("o1",oz,"snap_Vx_b")
  call to_history("d1",dz,"snap_Vx_b")
  call to_history("n2",nx,"snap_Vx_b")
  call to_history("o2",ox,"snap_Vx_b")
  call to_history("d2",dx,"snap_Vx_b")
  call to_history("n3",nsnap,"snap_Vx_b")
  call to_history("o3",0,"snap_Vx_b")
  call to_history("d3",1,"snap_Vx_b")
  call to_history("label1","z","snap_Vx_b")
  call to_history("label2","x","snap_Vx_b")
  call to_history("label3","t","snap_Vx_b")
  ! snapshot of z velocities output parameters
  call to_history("n1",nz,"snap_Vz_b")
  call to_history("o1",oz,"snap_Vz_b")
  call to_history("d1",dz,"snap_Vz_b")
  call to_history("n2",nx,"snap_Vz_b")
  call to_history("o2",ox,"snap_Vz_b")
  call to_history("d2",dx,"snap_Vz_b")
  call to_history("n3",nsnap,"snap_Vz_b")
  call to_history("o3",0,"snap_Vz_b")
  call to_history("d3",1,"snap_Vz_b")
  call to_history("label1","z","snap_Vz_b")
  call to_history("label2","x","snap_Vz_b")
  call to_history("label3","t","snap_Vz_b")
  ! snapshot of pressure output parameters
  call to_history("n1",nz,"snap_P_b")
  call to_history("o1",oz,"snap_P_b")
  call to_history("d1",dz,"snap_P_b")
  call to_history("n2",nx,"snap_P_b")
  call to_history("o2",ox,"snap_P_b")
  call to_history("d2",dx,"snap_P_b")
  call to_history("n3",nsnap,"snap_P_b")
  call to_history("o3",0,"snap_P_b")
  call to_history("d3",1,"snap_P_b")
  call to_history("label1","z","snap_P_b")
  call to_history("label2","x","snap_P_b")
  call to_history("label3","t","snap_P_b")

  ! Allocation of the input/output matrices
  if(lame)then
    allocate(lambda0(nz,nx))
    allocate(mi0(nz,nx))
    allocate(dlambda(nz,nx))
    allocate(dmi(nz,nx))
  else
    allocate(vp0(nz,nx))
    allocate(vs0(nz,nx))
    allocate(dvp(nz,nx))
    allocate(dvs(nz,nx))
  endif
  allocate(rho0(nz,nx))
  allocate(drho(nz,nx))
  allocate(seis_Vx(nseis,(nx-2*abc)))
  allocate(seis_Vz(nseis,(nx-2*abc)))
  allocate(seis_P(nseis,(nx-2*abc)))
  allocate(snap_Vx(nz,nx,nsnap))
  allocate(snap_Vz(nz,nx,nsnap))
  allocate(snap_P(nz,nx,nsnap))
  allocate(seis_Vx_b(nseis,(nx-2*abc)))
  allocate(seis_Vz_b(nseis,(nx-2*abc)))
  allocate(seis_P_b(nseis,(nx-2*abc)))
  allocate(snap_Vx_b(nz,nx,nsnap))
  allocate(snap_Vz_b(nz,nx,nsnap))
  allocate(snap_P_b(nz,nx,nsnap))

  ! Initialize values
  if(lame)then
    lambda0=0.0
    mi0=0.0
    dlambda=0.0
    dmi=0.0
  else
    vp0=0.0
    vs0=0.0
    dvp=0.0
    dvs=0.0
  endif
  rho0=0.0
  drho=0.0
  !Background outputs
  seis_Vx=0.0
  seis_Vz=0.0
  seis_P =0.0
  snap_Vx=0.0
  snap_Vz=0.0
  !Scatterer outputs
  seis_Vx=0.0
  seis_Vz=0.0
  seis_P =0.0
  snap_Vx=0.0
  snap_Vz=0.0

  if(lame)then
    call sep_read(lambda0)
    call sep_read(mi0,"mi0")
    call sep_read(rho0,"rho0")
    call sep_read(dlambda, "dlambda")
    call sep_read(dmi,"dmi")
    call sep_read(drho,"drho")
  else
    call sep_read(vp0)
    call sep_read(vs0,"vs0")
    call sep_read(rho0,"rho0")
    call sep_read(dvp, "dvp")
    call sep_read(dvs,"dvs")
    call sep_read(drho,"drho")
  endif

  nt = int(Tt/dt)

  ! calculate source length in time steps (making it 4 times the average
  ! distance between side lobes in a Ricker wavelet)
  nsource=4*nint(sqrt(6.0)/3.14159265/fp/dt)

  ! Write outputs to the terminal
  write(0,*)"============================================================="
  write(0,*)"BEGIN OF SUBROUTINE CALLS"
  write(0,*)"============================================================="
  write(0,*)"Values of the stability conditions:"
  write(0,*)"dx=dy=",dx
  write(0,*)"dt=",dt
  write(0,*)"source length=",nsource
  write(0,*)"seismogram length=",nseis
  write(0,*)"OMP using",nodes,"threads out of",node
  write(0,*)"============================================================="

  ! Initialize the elastic modeling
  call elastic_init(nz,dz,nx,dx,nt,dt,xsource,zsource,nsource,fp,nsnap,nseis,abc,zrec,free_surface)


  if(lame)then
    stat = elastic_op(lambda0,mi0,rho0,dlambda,dmi,drho,seis_Vx,seis_Vz,seis_P, &
                      seis_Vx_b, seis_Vz_b, seis_P_b, snap_Vx, snap_Vz, snap_P, &
                      snap_Vx_b, snap_Vz_b, snap_P_b, lame)
  else
    stat = elastic_op(vp0,vs0,rho0,dvp,dvs,drho,seis_Vx,seis_Vz,seis_P, &
                      seis_Vx_b, seis_Vz_b, seis_P_b, snap_Vx, snap_Vz, snap_P, &
                      snap_Vx_b, snap_Vz_b, snap_P_b, lame)
  endif
  if (stat) then
    write(0,*) "ERROR ON ELASTIC_OP FUNCTION CALL ==> EXITING"
    stop
  endif

  ! Write outputs
  call sep_write(seis_Vx)
  call sep_write(seis_Vz,"seis_Vz")
  call sep_write(seis_P,"seis_P")
  call sep_write(snap_Vx,"snap_Vx")
  call sep_write(snap_Vz,"snap_Vz")
  call sep_write(snap_P,"snap_P")
  call sep_write(seis_Vx_b,"seis_Vx_b")
  call sep_write(seis_Vz_b,"seis_Vz_b")
  call sep_write(seis_P_b,"seis_P_b")
  call sep_write(snap_Vx_b,"snap_Vx_b")
  call sep_write(snap_Vz_b,"snap_Vz_b")
  call sep_write(snap_P_b,"snap_P_b")

  call sep_close

  end_time = omp_get_wtime()

  ! Write outputs to the terminal
  write(0,*)"============================================================="
  write(0,*)"TOTAL ELAPSED TIME"
  write(0,*)"============================================================="
  write(0,*)"Total=",end_time - start_time, "seconds."
  write(0,*)"============================================================="

end program
