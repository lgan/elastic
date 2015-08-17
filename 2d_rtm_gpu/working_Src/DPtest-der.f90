! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! Dot product tests
program Dot_product_test

  use sep
  use rtmlib_mod_double
  use dptest_mod

implicit none

  integer            :: n1, n2, n_model, n_data, n_model2
  real, dimension(2) :: dot1, dot2
  real               :: dx
  real*8, dimension(:,:), allocatable :: rho_Vx,rho_Vz,l,mi,l2mi
  real*8, dimension(:,:), allocatable :: drho_Vx,drho_Vz,dl,dmi,dl2mi

  call sep_init()

  call from_param("n1", n1)
  call from_param("n2", n2)

  n_model = n1*n2
  n_model2 = n_model*2
  n_data = n1*n2

  allocate(rho_Vx(n1,n2))
  allocate(rho_Vz(n1,n2))
  allocate(l(n1,n2))
  allocate(mi(n1,n2))
  allocate(l2mi(n1,n2))
  allocate(drho_Vx(n1,n2))
  allocate(drho_Vz(n1,n2))
  allocate(dl(n1,n2))
  allocate(dmi(n1,n2))
  allocate(dl2mi(n1,n2))

  call random_number(rho_Vx)
  call random_number(rho_Vz)
  call random_number(l)
  call random_number(mi)
  call random_number(l2mi)
  call random_number(drho_Vx)
  call random_number(drho_Vz)
  call random_number(dl)
  call random_number(dmi)
  call random_number(dl2mi)

  dx=4.0

!  rho_Vx = 1.
!  rho_Vz = 1.
!  l = 1.
!  mi = 1.
!  l2mi = 1.
!  drho_Vx = 1.
!  drho_Vz = 1.
!  dl = 1.
!  dmi = 1.
!  dl2mi = 1.
  
  call rtmmod_init(n1,n2,dx,rho_Vx,rho_Vz,l,mi,l2mi,drho_Vx,drho_Vz,dl,dmi,dl2mi)

  call dptest(funAx, n_model, n_data, dot1, dot2)

  write(0,*) "Ax"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funBx, n_model, n_data, dot1, dot2)

  write(0,*) "Bx"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funAz, n_model, n_data, dot1, dot2)

  write(0,*) "Az"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funBz, n_model, n_data, dot1, dot2)

  write(0,*) "Bz"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funCx, n_model, n_data, dot1, dot2)

  write(0,*) "Cx"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funCz, n_model, n_data, dot1, dot2)

  write(0,*) "Cz"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funD, n_model, n_data, dot1, dot2)

  write(0,*) "D"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funE, n_model, n_data, dot1, dot2)

  write(0,*) "E"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funF, n_model, n_data, dot1, dot2)

  write(0,*) "F"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funG, n_model, n_data, dot1, dot2)

  write(0,*) "G"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funH, n_model, n_data, dot1, dot2)

  write(0,*) "H"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funJ, n_model, n_data, dot1, dot2)

  write(0,*) "J"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funT, n_model2, n_data, dot1, dot2)

  write(0,*) "T"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funS1, n_model2, n_data, dot1, dot2)

  write(0,*) "S1"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funS2, n_model2, n_data, dot1, dot2)

  write(0,*) "S2"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funS3, n_model2, n_data, dot1, dot2)

  write(0,*) "S3"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funS4, n_model2, n_data, dot1, dot2)

  write(0,*) "S4"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call dptest(funS5, n_model2, n_data, dot1, dot2)

  write(0,*) "S5"
  write(0,*) dot1
  write(0,*) dot2

  call to_history("dot1", dot1)
  call to_history("dot2", dot2)

  call sep_close()

end program
