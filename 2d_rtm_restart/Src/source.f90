! Source output
program source_program

  use sep
  use source_mod

  implicit none

  integer                           :: nx,nz,nsource,it,nt
  real                              :: source
  real                              :: dt, fp, ot,dx
  real, allocatable, dimension(:)   :: data

  call sep_init()

  call from_param("dt", dt)
  call from_param("nt", nt)
  call from_param("fp", fp)
  call from_param("ot", ot)

  call to_history("n1", nt)
  call to_history("o1", ot)
  call to_history("d1", dt)
  call to_history("label1","t")

  allocate(data(nt))
  data = 0.
  source=0.0
  nx=1
  nz=1
  nsource=nt
  call source_init(fp,dt,dx,nsource,nz,nx)

  do it=1,nt
    call source_op2(it,source)
    data(it) = data(it) + source
  enddo
  call sep_write(data)

  call sep_close()

end program
