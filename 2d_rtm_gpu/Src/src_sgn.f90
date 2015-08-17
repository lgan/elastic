program src_sgn
  use sep

  implicit none

 
  integer :: tt,nt
  real, dimension(:), allocatable :: src
  real::PI, a, src_fpeak, dt, src_t



  write(0,*) "src_sgn.x < dt= nt= src_fpeak= > output.H"

  call sep_init()



  call from_param ("dt", dt)
  call from_param ("nt", nt)
  call from_param ("src_fpeak", src_fpeak)


  allocate(src(nt))


  
  PI=3.1415926 
  a=(PI*src_fpeak)**2 
  do tt=1,nt
     src_t=((tt-1.)*dt-1./src_fpeak)**2.
     src(tt)=-(1.-2.*a*src_t)*exp(-a*src_t)
  end do



  call to_history('n1',nt)
  call to_history('d1',dt)
  call to_history('o1',0.)



  call sep_write(src)

  deallocate(src)

  call sep_close()
end program

