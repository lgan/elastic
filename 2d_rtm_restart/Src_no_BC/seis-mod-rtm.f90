! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! seismogram module
module seis_mod_rtm

  use sep

  implicit none
  integer, private :: nx,nt,nz
  real, private    :: dx,dt,dz
  real, dimension(:,:), pointer, private :: Vx_rec,Vz_rec
  real, dimension(:,:), pointer, private :: sigmaxx_rec,sigmazz_rec,sigmaxz_rec
  real, dimension(:,:), pointer, private :: dm1,dm2,dm3
  logical,private                        :: lame

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine seismo(nx_in,dx_in,nt_in,dt_in,Vx_rec_in,Vz_rec_in,sigmaxx_rec_in,sigmazz_rec_in,sigmaxz_rec_in)
    integer,intent(in)           :: nx_in,nt_in
    real,intent(in)              :: dx_in,dt_in
    real, dimension(:,:), target :: Vx_rec_in,Vz_rec_in,sigmaxx_rec_in,sigmazz_rec_in,sigmaxz_rec_in

    nx = nx_in
    dx = dx_in
    nt = nt_in
    dt = dt_in
    Vx_rec => Vx_rec_in
    Vz_rec => Vz_rec_in
    sigmaxx_rec => sigmaxx_rec_in
    sigmazz_rec => sigmazz_rec_in
    sigmaxz_rec => sigmaxz_rec_in

    call to_history("n1", nt, "Vx_rec")
    call to_history("o1", 0., "Vx_rec")
    call to_history("d1", dt, "Vx_rec")
    call to_history("n2", nx, "Vx_rec")
    call to_history("o2", 0., "Vx_rec")
    call to_history("d2", dx, "Vx_rec")
         
    call sep_write(Vx_rec, "Vx_rec")

    call to_history("n1", nt, "Vz_rec")
    call to_history("o1", 0., "Vz_rec")
    call to_history("d1", dt, "Vz_rec")
    call to_history("n2", nx, "Vz_rec")
    call to_history("o2", 0., "Vz_rec")
    call to_history("d2", dx, "Vz_rec")
         
    call sep_write(Vz_rec, "Vz_rec")

    call to_history("n1", nt, "sigmaxx_rec")
    call to_history("o1", 0., "sigmaxx_rec")
    call to_history("d1", dt, "sigmaxx_rec")
    call to_history("n2", nx, "sigmaxx_rec")
    call to_history("o2", 0., "sigmaxx_rec")
    call to_history("d2", dx, "sigmaxx_rec")
         
    call sep_write(sigmaxx_rec, "sigmaxx_rec")

    call to_history("n1", nt, "sigmazz_rec")
    call to_history("o1", 0., "sigmazz_rec")
    call to_history("d1", dt, "sigmazz_rec")
    call to_history("n2", nx, "sigmazz_rec")
    call to_history("o2", 0., "sigmazz_rec")
    call to_history("d2", dx, "sigmazz_rec")
         
    call sep_write(sigmazz_rec, "sigmazz_rec")

    call to_history("n1", nt, "sigmaxz_rec")
    call to_history("o1", 0., "sigmaxz_rec")
    call to_history("d1", dt, "sigmaxz_rec")
    call to_history("n2", nx, "sigmaxz_rec")
    call to_history("o2", 0., "sigmaxz_rec")
    call to_history("d2", dx, "sigmaxz_rec")
         
    call sep_write(sigmaxz_rec, "sigmaxz_rec")

  endsubroutine

  subroutine img(nx_in,dx_in,nz_in,dz_in,dm1_in,dm2_in,dm3_in,lame_in)
    integer,intent(in)           :: nx_in,nz_in
    real,intent(in)              :: dx_in,dz_in
    real, dimension(:,:), target :: dm1_in,dm2_in,dm3_in
    logical                      :: lame_in

    nx = nx_in
    dx = dx_in
    nz = nz_in
    dz = dz_in
    dm1 => dm1_in
    dm2 => dm2_in
    dm3 => dm3_in
    lame = lame_in

    call to_history("n1", nz, "drho")
    call to_history("o1", 0., "drho")
    call to_history("d1", dz, "drho")
    call to_history("n2", nx, "drho")
    call to_history("o2", 0., "drho")
    call to_history("d2", dx, "drho")
         
    call sep_write(dm1, "drho")

    if (lame) then
      call to_history("n1", nz, "dlambda")
      call to_history("o1", 0., "dlambda")
      call to_history("d1", dz, "dlambda")
      call to_history("n2", nx, "dlambda")
      call to_history("o2", 0., "dlambda")
      call to_history("d2", dx, "dlambda")
      call sep_write(dm2, "dlambda")

      call to_history("n1", nz, "dmu")
      call to_history("o1", 0., "dmu")
      call to_history("d1", dz, "dmu")
      call to_history("n2", nx, "dmu")
      call to_history("o2", 0., "dmu")
      call to_history("d2", dx, "dmu")    
      call sep_write(dm3, "dmu")
    else
      call to_history("n1", nz, "dvp")
      call to_history("o1", 0., "dvp")
      call to_history("d1", dz, "dvp")
      call to_history("n2", nx, "dvp")
      call to_history("o2", 0., "dvp")
      call to_history("d2", dx, "dvp")
      call sep_write(dm2, "dvp")
         
      call to_history("n1", nz, "dvs")
      call to_history("o1", 0., "dvs")
      call to_history("d1", dz, "dvs")
      call to_history("n2", nx, "dvs")
      call to_history("o2", 0., "dvs")
      call to_history("d2", dx, "dvs")    
      call sep_write(dm3, "dvs")
    endif


  endsubroutine

endmodule
