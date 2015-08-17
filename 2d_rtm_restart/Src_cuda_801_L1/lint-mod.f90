! Linear interpolation operator module
module lint_mod

  implicit none
  real, private                         :: dt, dtseis
  integer, private                      :: nrec, nseis
  real, private, dimension(:,:), pointer:: seis_Vx, seis_Vz

  contains

  subroutine lint_init(dt_in,dtseis_in,seis_Vx_in, seis_Vz_in,nrec_in,nseis_in)
    integer                           :: nrec_in, nseis_in, nout_in
    real                              :: dt_in, dtseis_in
    real, dimension(:,:), target      :: seis_Vx_in, seis_Vz_in

    nrec = nrec_in
    nseis = nseis_in
    dt = dt_in
    dtseis = dtseis_in
    seis_Vx => seis_Vx_in
    seis_Vz => seis_Vz_in

  end subroutine

  subroutine lint_op(adj, add, data1, data2, it)
    logical,intent(in)          :: adj, add
    integer                     :: it, ix, iseis
    real, dimension(:,:)        :: data1, data2
    real                        :: temp, f, g

    if (adj) then
      if (add==0) data1=0.
      do ix=1,nrec
        temp = it*dt/dtseis
        iseis = temp
        f = temp - iseis
        g = 1.0 - f
        data1(iseis+1,ix) =  data1(iseis+1,ix) + g*seis_Vx(it,ix) 
        data1(iseis+2,ix) =  data1(iseis+2,ix) + f*seis_Vx(it,ix) 
        data2(iseis+1,ix) =  data2(iseis+1,ix) + g*seis_Vz(it,ix) 
        data2(iseis+2,ix) =  data2(iseis+2,ix) + f*seis_Vz(it,ix) 
      enddo
    else
      do ix=1,nrec
        temp = (it*dt)/dtseis
        iseis = temp
        f = temp-iseis
        g = 1.0 - f
        data1(it,ix) = data1(it,ix) + (g*seis_Vx(iseis+1,ix) + f*seis_Vx(iseis+2,ix))
        data2(it,ix) = data2(it,ix) + (g*seis_Vz(iseis+1,ix) + f*seis_Vz(iseis+2,ix))
      enddo
    endif

  endsubroutine

end module
