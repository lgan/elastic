! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! Dot product tests
program Dot_product_test_pml

  use sep
  use pml_mod_double
  use dptest_mod

implicit none

  integer            :: n1, n2, n_model, n_data, abc
  real               :: dx,dt
  real               :: vel
  real, dimension(2) :: dot1, dot2
  real*16, allocatable, dimension(:,:) :: vp

  call sep_init()

  call from_param("n1", n1)
  call from_param("n2", n2)
  call from_param("dx", dx)
  call from_param("dt", dt)
  call from_param("abc", abc)
  call from_param("vel", vel)

  allocate(vp(n1,n2))

  vp = vel

  n_model = n1*n2
  n_data = n1*n2

  call pml_init(n1, n2, dx, dt, abc, vp)

  call dptest(boundary_op0_x, n_model, n_data, dot1, dot2)

  write(0,*) "Derivative for pml_0_x"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(boundary_op0_x2, n_model, n_data, dot1, dot2)

  write(0,*) "Derivative for pml_0_x2"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(boundary_op0_z, n_model, n_data, dot1, dot2)

  write(0,*) "Derivative for pml_0_z"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(boundary_op0_z2, n_model, n_data, dot1, dot2)

  write(0,*) "Derivative for pml_1_z2"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call sep_close()

end program
