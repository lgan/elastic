! Apply Linear interpolation operator
program Linear_interpolation

  use sep
  use lint_mod_template

  implicit none

  integer                           :: it,nrec,nseis,nout
  real                              :: dt, dtseis, oseis, orec, drec
  logical                           :: adj, add
  real, dimension(:,:), allocatable :: seisinput,seisoutput

  call sep_init()

  call from_history("n1", nseis)
  call from_history("o1", oseis)
  call from_history("d1", dtseis)
  call from_history("n2", nrec)
  call from_history("o2", orec)
  call from_history("d2", drec)

  call from_param("adj", adj)
  call from_param("add", add)
  call from_param("dt",dt)

  nout=0.5+(nseis*dtseis)/dt

  call to_history("n1", nout)
  call to_history("o1", oseis)
  call to_history("d1", dt)
  call to_history("n2", nrec)
  call to_history("o2", orec)
  call to_history("d2", drec)

  allocate(seisinput(nseis,nrec))
  allocate(seisoutput(nout,nrec))
  seisinput = 0.
  seisoutput = 0.

  call sep_read(seisinput)
  
  call lint_init(dt,dtseis,seisinput,nrec,nseis,nout)

  do it=1,nout
    call lint_op(adj,add,seisoutput,it)
  enddo

  call sep_write(seisoutput)
  call sep_write(seisinput,"seisinput.H")
  call sep_close()

end program
