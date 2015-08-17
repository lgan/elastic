! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! June, 25th, 2014

! Main program to pick maximun amplitudes
program pickmax_program

  use sep

  implicit none

  integer                               :: nx,nt,it,ix
  integer                               :: minx,maxx,itmin,itmax
  real                                  :: mint,maxt
  integer                               :: size_x
  real                                  :: ot,dt,ox,dx
  real, allocatable, dimension(:,:)     :: seis
  real, allocatable, dimension(:)       :: val

  call sep_init()

  ! Read SEP header parameter
  ! Model parameters
  call from_history("n1",nt)
  call from_history("o1",ot)
  call from_history("d1",dt)
  call from_history("n2",nx)
  call from_history("o2",ox)
  call from_history("d2",dx)
  ! Picking parameters
  call from_param("minx",minx,1)           !1st  X trace to read, default=  1
  call from_param("maxx",maxx,nx)          !last X trace to read, default= nx
  call from_param("mint",mint,0.0)         !1st  t value to read, default=  0
  call from_param("maxt",maxt,nt*dt)       !last t value to read, default= nt*dt

  size_x = maxx-minx+1
  itmin = 1 + mint/dt
  itmax = maxt/dt

  ! Write SEP header parameter
  ! valp.H, vals.H output parameters
  call to_history("n1",size_x)
  call to_history("o1",ox)
  call to_history("d1",dx)

  ! Allocation of the input/output matrices
  allocate(seis(nt,nx))
  allocate(val(size_x))

  ! Initialize values
  seis=0.0
  val =0.0

  ! Write outputs to the terminal
  write(0,*)"============================================================="
  write(0,*)"AUTOMATIC PICKING OF MAX AMPLITUDES PER TRACE"
  write(0,*)"============================================================="
  write(0,*)"Total number of traces=",nx
  write(0,*)"Number of traces to read=",size_x
  write(0,*)"Total seismogram length=",nt
  write(0,*)"Seismogram window length=",(maxt-mint)
  write(0,*)"============================================================="

  ! Read seismograms
  call sep_read(seis)

  write(0,*)"SEISMOGRAM READ => OK"

  write(0,*)"itmin=",itmin
  write(0,*)"itmax=",itmax

  do ix=minx,maxx 
     val(ix)=0.0
     do it=itmin,itmax
       !if(abs(seis(it,ix)).gt.val(ix)) then
         val(ix)=val(ix)+abs(seis(it,ix))**2
       !endif
     enddo
     write(0,*)"val=",val(ix),"ix=",ix
   enddo
       
  ! Write outputs
  call sep_write(val)

  call sep_close

end program
