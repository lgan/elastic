! Linear interpolation operator module
module lint_mod_template

  implicit none
  real, private                         :: dt, dtseis
  integer, private                      :: nrec, nseis, nout
  real, private, dimension(:,:), pointer:: seis

  contains

  subroutine lint_init(dt_in,dtseis_in,seis_in,nrec_in,nseis_in,nout_in)
    integer                           :: nrec_in, nseis_in, nout_in
    real                              :: dt_in, dtseis_in
    real, dimension(:,:), target      :: seis_in

    nrec = nrec_in
    nseis = nseis_in
    dt = dt_in
    nout = nout_in
    dtseis = dtseis_in
    seis => seis_in

  end subroutine

  subroutine lint_op(adj, add, data, it)
    logical,intent(in)          :: adj, add
    integer                     :: it, ix, iseis
    real, dimension(:,:)        :: data
    real                        :: temp, f, g

    if (adj) then
      if (add==0) data=0.
      do ix=1,nrec
        temp = it*dt/dtseis
        iseis = temp
        f = temp - iseis
        g = 1.0 - f
        data(iseis+1,ix) =  data(iseis+1,ix) + g*seis(nout-it+1,ix) 
        data(iseis+2,ix) =  data(iseis+2,ix) + f*seis(nout-it+1,ix) 
      enddo
    else
      do ix=1,nrec
        temp = (it*dt)/dtseis
        iseis = temp
        f = temp-iseis
        g = 1.0 - f
        data(nout-it+1,ix) = data(nout-it+1,ix) + (g*seis(iseis+1,ix) + f*seis(iseis+2,ix))
      enddo
    endif

  endsubroutine

end module
