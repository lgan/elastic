module fourier

  implicit none
  include "fftw3.f"

  complex, dimension(:), allocatable, private  :: in1, out1
  integer, private                             :: n1
  integer*8, private                           :: plan1f, plan1b
  contains

  subroutine ft_init(n1_in)
    integer               :: n1_in
    n1 = n1_in
    allocate(in1(n1), out1(n1))
    in1 = 0.
    out1 = 0.
    call sfftw_plan_dft_1d(plan1f, n1, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE)
    call sfftw_plan_dft_1d(plan1b, n1, out1, in1, FFTW_BACKWARD, FFTW_ESTIMATE)
  end subroutine

  subroutine ft1axis(adj, sign1, cdata)
  ! Forward and Ajoint/Inverse Fourier Transform on axis 1
    logical, intent(in)        :: adj
    integer                    :: i2
    complex, dimension(:)      :: cdata
    real                       :: sign1
    if(adj) then
      out1(1:(n1+1)/2) = cdata(n1-(n1+1)/2+1:n1)
      out1((n1+1)/2+1:n1) = cdata(1:n1-(n1+1)/2)
      call sfftw_execute_dft(plan1b, out1, in1)
      cdata(:) = in1
    else
      in1 = cdata(:)
      call sfftw_execute_dft(plan1f, in1, out1)
      cdata(:) = sign1*out1/sqrt(1.*n1)
      cdata(1:n1-(n1+1)/2) = out1((n1+1)/2+1:n1)
      cdata(n1-(n1+1)/2+1:n1) = out1(1:(n1+1)/2)
    endif
    cdata = sign1*cdata/sqrt(1.*n1)
  end subroutine

end module
