! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 20th, 2015

! seismogram module
module snap_mod_rtm

  use sep

  implicit none
  integer, private :: nz,nx,nt,abc
  real, private    :: dz,dx,dt
  real, dimension(:,:,:), pointer, private :: Vx0t,Vz0t,Vx0dx,Vz0dz,Vx0dz,Vz0dx 
  real, dimension(:,:,:), pointer, private :: Vx0,Vz0,sigmaxx0,sigmazz0,sigmaxz0

  integer, private :: it,countsnap,nsnap,dsnap

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine snap_movie(nz_in,dz_in,nx_in,dx_in,nt_in,dt_in,abc_in,Vx0t_in,Vz0t_in,Vx0dx_in,Vz0dz_in,Vx0dz_in,Vz0dx_in,Vx0_in,Vz0_in,sigmaxx0_in,sigmazz0_in,sigmaxz0_in)
    integer,intent(in)           :: nz_in,nx_in,nt_in,abc_in
    real,intent(in)              :: dz_in,dx_in,dt_in
    real,dimension(:,:,:),target :: Vx0t_in,Vz0t_in,Vx0dx_in,Vz0dz_in,Vx0dz_in,Vz0dx_in
    real,dimension(:,:,:),target :: Vx0_in,Vz0_in,sigmaxx0_in,sigmazz0_in,sigmaxz0_in

    real,allocatable,dimension(:,:,:) :: data0

    nz = nz_in
    dz = dz_in
    nx = nx_in
    dx = dx_in
    nt = nt_in
    dt = dt_in
    abc = abc_in
    Vx0t => Vx0t_in
    Vz0t => Vz0t_in
    Vx0dx => Vx0dx_in
    Vz0dz => Vz0dz_in
    Vx0dz => Vx0dz_in
    Vz0dx => Vz0dx_in
    Vx0 => Vx0_in
    Vz0 => Vz0_in
    sigmaxx0 => sigmaxx0_in
    sigmazz0 => sigmazz0_in
    sigmaxz0 => sigmaxz0_in
    nsnap=100
    dsnap=nt/nsnap

    allocate(data0(nz-(2*abc+6),5*(nx-2*abc),nsnap))

    data0 = 0.
    countsnap=2
    do it=1,nt
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        data0(:,1:(nx-2*abc),countsnap) = Vx0t(abc+6:nz-abc,abc+1:nx-abc,it)
        data0(:,(nx-2*abc)+1:2*(nx-2*abc),countsnap) = Vx0(abc+6:nz-abc,abc+1:nx-abc,it)
        data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),countsnap) = data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),countsnap-1) + (Vx0t(abc+6:nz-abc,abc+1:nx-abc,it)*Vx0(abc+6:nz-abc,abc+1:nx-abc,it)) +&
                                                                                                            (Vz0t(abc+6:nz-abc,abc+1:nx-abc,it)*Vz0(abc+6:nz-abc,abc+1:nx-abc,it))
        data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),countsnap) = Vz0t(abc+6:nz-abc,abc+1:nx-abc,it)
        data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),countsnap) = Vz0(abc+6:nz-abc,abc+1:nx-abc,it)
       countsnap = countsnap+1
      endif
    enddo
    data0(:,1:(nx-2*abc),:) = data0(:,1:(nx-2*abc),:)/maxval(data0(:,1:(nx-2*abc),:))
    data0(:,(nx-2*abc)+1:2*(nx-2*abc),:) = data0(:,(nx-2*abc)+1:2*(nx-2*abc),:)/maxval(data0(:,(nx-2*abc)+1:2*(nx-2*abc),:))
    data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:) = data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:)/maxval(data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:))
    data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:) = data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:)/maxval(data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:))
    data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:) = data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:)/maxval(data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:))

    call to_history("n1", nz-2*abc-6, "Dat/r_rho_time.H")
    call to_history("o1", 0., "Dat/r_rho_time.H")
    call to_history("d1", dz, "Dat/r_rho_time.H")
    call to_history("n2", 5*(nx-2*abc), "Dat/r_rho_time.H")
    call to_history("o2", 0., "Dat/r_rho_time.H")
    call to_history("d2", dx, "Dat/r_rho_time.H")
    call to_history("n3", nsnap, "Dat/r_rho_time.H")
    call to_history("o3", 0., "Dat/r_rho_time.H")
    call to_history("d3", dt*dsnap, "Dat/r_rho_time.H")
         
    call sep_write(data0,"Dat/r_rho_time.H")

    data0 = 0.
    countsnap=2
    do it=1,nt
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        data0(:,1:(nx-2*abc),countsnap) = -(Vx0dx(abc+6:nz-abc,abc+1:nx-abc,it)+Vz0dz(abc+6:nz-abc,abc+1:nx-abc,it))
        data0(:,(nx-2*abc)+1:2*(nx-2*abc),countsnap) = sigmaxx0(abc+6:nz-abc,abc+1:nx-abc,it)
        data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),countsnap) = data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),countsnap-1) + (-(Vx0dx(abc+6:nz-abc,abc+1:nx-abc,it)+Vz0dz(abc+6:nz-abc,abc+1:nx-abc,it))*sigmaxx0(abc+6:nz-abc,abc+1:nx-abc,it)) +&
                                                                                (-(Vx0dx(abc+6:nz-abc,abc+1:nx-abc,it)+Vz0dz(abc+6:nz-abc,abc+1:nx-abc,it))*sigmazz0(abc+6:nz-abc,abc+1:nx-abc,it))
        data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),countsnap) = -(Vx0dx(abc+6:nz-abc,abc+1:nx-abc,it)+Vz0dz(abc+6:nz-abc,abc+1:nx-abc,it))
        data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),countsnap) = sigmazz0(abc+6:nz-abc,abc+1:nx-abc,it)
       countsnap = countsnap+1
      endif
    enddo
    data0(:,1:(nx-2*abc),:) = data0(:,1:(nx-2*abc),:)/maxval(data0(:,1:(nx-2*abc),:))
    data0(:,(nx-2*abc)+1:2*(nx-2*abc),:) = data0(:,(nx-2*abc)+1:2*(nx-2*abc),:)/maxval(data0(:,(nx-2*abc)+1:2*(nx-2*abc),:))
    data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:) = data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:)/maxval(data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:))
    data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:) = data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:)/maxval(data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:))
    data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:) = data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:)/maxval(data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:))

    call to_history("n1", nz-2*abc-6, "Dat/r_l_time.H")
    call to_history("o1", 0., "Dat/r_l_time.H")
    call to_history("d1", dz, "Dat/r_l_time.H")
    call to_history("n2", 5*(nx-2*abc), "Dat/r_l_time.H")
    call to_history("o2", 0., "Dat/r_l_time.H")
    call to_history("d2", dx, "Dat/r_l_time.H")
    call to_history("n3", nsnap, "Dat/r_l_time.H")
    call to_history("o3", 0., "Dat/r_l_time.H")
    call to_history("d3", dt*dsnap, "Dat/r_l_time.H")
         
    call sep_write(data0,"Dat/r_l_time.H")

    deallocate(data0)
    allocate(data0(nz-2*abc-6,7*(nx-2*abc),nsnap))
    data0 = 0.
    countsnap=2
    do it=1,nt
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        data0(:,1:(nx-2*abc),countsnap) = -2.0*Vx0dx(abc+6:nz-abc,abc+1:nx-abc,it)
        data0(:,(nx-2*abc)+1:2*(nx-2*abc),countsnap) = sigmaxx0(abc+6:nz-abc,abc+1:nx-abc,it)
        data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),countsnap) = data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),countsnap-1) + (-2.0*Vx0dx(abc+6:nz-abc,abc+1:nx-abc,it)*sigmaxx0(abc+6:nz-abc,abc+1:nx-abc,it)+&
                                                                                                             -2.0*Vz0dz(abc+6:nz-abc,abc+1:nx-abc,it)*sigmazz0(abc+6:nz-abc,abc+1:nx-abc,it)+&
                                                                           -(Vz0dx(abc+6:nz-abc,abc+1:nx-abc,it)+Vx0dz(abc+6:nz-abc,abc+1:nx-abc,it))*sigmaxz0(abc+6:nz-abc,abc+1:nx-abc,it))
        data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),countsnap) = -2.0*Vz0dz(abc+6:nz-abc,abc+1:nx-abc,it)
        data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),countsnap) = sigmazz0(abc+6:nz-abc,abc+1:nx-abc,it)
        data0(:,5*(nx-2*abc)+1:6*(nx-2*abc),countsnap) = -(Vz0dx(abc+6:nz-abc,abc+1:nx-abc,it)+Vx0dz(abc+6:nz-abc,abc+1:nx-abc,it))
        data0(:,6*(nx-2*abc)+1:7*(nx-2*abc),countsnap) = sigmaxz0(abc+6:nz-abc,abc+1:nx-abc,it)
       countsnap = countsnap+1
      endif
    enddo
    data0(:,1:(nx-2*abc),:) = data0(:,1:(nx-2*abc),:)/maxval(data0(:,1:(nx-2*abc),:))
    data0(:,(nx-2*abc)+1:2*(nx-2*abc),:) = data0(:,(nx-2*abc)+1:2*(nx-2*abc),:)/maxval(data0(:,(nx-2*abc)+1:2*(nx-2*abc),:))
    data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:) = data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:)/maxval(data0(:,2*(nx-2*abc)+1:3*(nx-2*abc),:))
    data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:) = data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:)/maxval(data0(:,3*(nx-2*abc)+1:4*(nx-2*abc),:))
    data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:) = data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:)/maxval(data0(:,4*(nx-2*abc)+1:5*(nx-2*abc),:))
    data0(:,5*(nx-2*abc)+1:6*(nx-2*abc),:) = data0(:,5*(nx-2*abc)+1:6*(nx-2*abc),:)/maxval(data0(:,5*(nx-2*abc)+1:6*(nx-2*abc),:))
    data0(:,6*(nx-2*abc)+1:7*(nx-2*abc),:) = data0(:,6*(nx-2*abc)+1:7*(nx-2*abc),:)/maxval(data0(:,6*(nx-2*abc)+1:7*(nx-2*abc),:))

    call to_history("n1", nz-2*abc-6, "Dat/r_mi_time.H")
    call to_history("o1", 0., "Dat/r_mi_time.H")
    call to_history("d1", dz, "Dat/r_mi_time.H")
    call to_history("n2", 7*(nx-2*abc), "Dat/r_mi_time.H")
    call to_history("o2", 0., "Dat/r_mi_time.H")
    call to_history("d2", dx, "Dat/r_mi_time.H")
    call to_history("n3", nsnap, "Dat/r_mi_time.H")
    call to_history("o3", 0., "Dat/r_mi_time.H")
    call to_history("d3", dt*dsnap, "Dat/r_mi_time.H")
         
    call sep_write(data0,"Dat/r_mi_time.H")

  endsubroutine

  subroutine snap_wfld(nz_in,dz_in,nx_in,dx_in,nt_in,dt_in,abc_in,Vx0_in,Vz0_in,sigmaxx0_in,sigmazz0_in,sigmaxz0_in)
    integer,intent(in)           :: nz_in,nx_in,nt_in,abc_in
    real,intent(in)              :: dz_in,dx_in,dt_in
    real,dimension(:,:,:),target :: Vx0_in,Vz0_in,sigmaxx0_in,sigmazz0_in,sigmaxz0_in

    real,allocatable,dimension(:,:,:) :: data0

    nz = nz_in
    dz = dz_in
    nx = nx_in
    dx = dx_in
    nt = nt_in
    dt = dt_in
    abc = abc_in
    Vx0 => Vx0_in
    Vz0 => Vz0_in
    sigmaxx0 => sigmaxx0_in
    sigmazz0 => sigmazz0_in
    sigmaxz0 => sigmaxz0_in
    nsnap=100
    dsnap=nt/nsnap

    allocate(data0(nz-(2*abc+6),(nx-2*abc),nsnap))

    data0 = 0.
    countsnap=1
    do it=1,nt
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        data0(:,:,countsnap) = Vx0(abc+6:nz-abc,abc+1:nx-abc,it)
        countsnap = countsnap+1
      endif
    enddo

    call to_history("n1", nz-2*abc-6, "Dat/Vx0.H")
    call to_history("o1", 0., "Dat/Vx0.H")
    call to_history("d1", dz, "Dat/Vx0.H")
    call to_history("n2", nx-2*abc, "Dat/Vx0.H")
    call to_history("o2", 0., "Dat/Vx0.H")
    call to_history("d2", dx, "Dat/Vx0.H")
    call to_history("n3", nsnap, "Dat/Vx0.H")
    call to_history("o3", 0., "Dat/Vx0.H")
    call to_history("d3", dt*dsnap, "Dat/Vx0.H")
         
    call sep_write(data0,"Dat/Vx0.H")

    data0 = 0.
    countsnap=1
    do it=1,nt
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        data0(:,:,countsnap) = Vz0(abc+6:nz-abc,abc+1:nx-abc,it)
        countsnap = countsnap+1
      endif
    enddo

    call to_history("n1", nz-2*abc-6, "Dat/Vz0.H")
    call to_history("o1", 0., "Dat/Vz0.H")
    call to_history("d1", dz, "Dat/Vz0.H")
    call to_history("n2", nx-2*abc, "Dat/Vz0.H")
    call to_history("o2", 0., "Dat/Vz0.H")
    call to_history("d2", dx, "Dat/Vz0.H")
    call to_history("n3", nsnap, "Dat/Vz0.H")
    call to_history("o3", 0., "Dat/Vz0.H")
    call to_history("d3", dt*dsnap, "Dat/Vz0.H")
         
    call sep_write(data0,"Dat/Vz0.H")

    data0 = 0.
    countsnap=1
    do it=1,nt
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        data0(:,:,countsnap) = sigmaxx0(abc+6:nz-abc,abc+1:nx-abc,it)
        countsnap = countsnap+1
      endif
    enddo

    call to_history("n1", nz-2*abc-6, "Dat/sigmaxx0.H")
    call to_history("o1", 0., "Dat/sigmaxx0.H")
    call to_history("d1", dz, "Dat/sigmaxx0.H")
    call to_history("n2", nx-2*abc, "Dat/sigmaxx0.H")
    call to_history("o2", 0., "Dat/sigmaxx0.H")
    call to_history("d2", dx, "Dat/sigmaxx0.H")
    call to_history("n3", nsnap, "Dat/sigmaxx0.H")
    call to_history("o3", 0., "Dat/sigmaxx0.H")
    call to_history("d3", dt*dsnap, "Dat/sigmaxx0.H")
         
    call sep_write(data0,"Dat/sigmaxx0.H")

    data0 = 0.
    countsnap=1
    do it=1,nt
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        data0(:,:,countsnap) = sigmazz0(abc+6:nz-abc,abc+1:nx-abc,it)
        countsnap = countsnap+1
      endif
    enddo

    call to_history("n1", nz-2*abc-6, "Dat/sigmazz0.H")
    call to_history("o1", 0., "Dat/sigmazz0.H")
    call to_history("d1", dz, "Dat/sigmazz0.H")
    call to_history("n2", nx-2*abc, "Dat/sigmazz0.H")
    call to_history("o2", 0., "Dat/sigmazz0.H")
    call to_history("d2", dx, "Dat/sigmazz0.H")
    call to_history("n3", nsnap, "Dat/sigmazz0.H")
    call to_history("o3", 0., "Dat/sigmazz0.H")
    call to_history("d3", dt*dsnap, "Dat/sigmazz0.H")
         
    call sep_write(data0,"Dat/sigmazz0.H")

    data0 = 0.
    countsnap=1
    do it=1,nt
      if (mod(it,dsnap)==0.and.countsnap.le.nsnap) then
        data0(:,:,countsnap) = sigmaxz0(abc+6:nz-abc,abc+1:nx-abc,it)
        countsnap = countsnap+1
      endif
    enddo

    call to_history("n1", nz-2*abc-6, "Dat/sigmaxz0.H")
    call to_history("o1", 0., "Dat/sigmaxz0.H")
    call to_history("d1", dz, "Dat/sigmaxz0.H")
    call to_history("n2", nx-2*abc, "Dat/sigmaxz0.H")
    call to_history("o2", 0., "Dat/sigmaxz0.H")
    call to_history("d2", dx, "Dat/sigmaxz0.H")
    call to_history("n3", nsnap, "Dat/sigmaxz0.H")
    call to_history("o3", 0., "Dat/sigmaxz0.H")
    call to_history("d3", dt*dsnap, "Dat/sigmaxz0.H")
         
    call sep_write(data0,"Dat/sigmaxz0.H")

  endsubroutine

endmodule

