! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! November, 10th, 2014

! boundary module
module pml_mod

  use sep

  implicit none
  integer, private   :: nx, nz, iz, ix, abc
  real, private      :: dt, dx
  real, private, allocatable, dimension(:) :: filter
  real, private, allocatable, dimension(:,:) :: d0, d0_xl, d0_xr, d0_zt, d0_zb
  real, private, allocatable, dimension(:,:) :: pml0_xl, pml0_xr, pml0_zt, pml0_zb
  real, private, allocatable, dimension(:,:) :: pml1_xl, pml1_xr, pml1_zt, pml1_zb
  real, pointer, private, dimension(:,:) :: vp
  logical, private :: free_surface

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine pml_init(nz_in,nx_in, dx_in, dt_in,abc_in, vp_in, free_surface_in)
    integer        :: nx_in, nz_in, abc_in
    real           :: R,dx_in,dt_in
    real, dimension(:,:),target :: vp_in
    logical        :: free_surface_in

    nx = nx_in
    nz = nz_in
    dx = dx_in
    dt = dt_in
    free_surface = free_surface_in
    abc= abc_in
    vp => vp_in

    allocate(filter(abc))
    allocate(d0(nz,nx))
    allocate(d0_xl(nz,abc))
    allocate(d0_xr(nz,abc))
    allocate(d0_zt(abc,nx))
    allocate(d0_zb(abc,nx))
    allocate(pml0_xl(nz,abc))
    allocate(pml0_xr(nz,abc))
    allocate(pml0_zt(abc,nx))
    allocate(pml0_zb(abc,nx))
    allocate(pml1_xl(nz,abc))
    allocate(pml1_xr(nz,abc))
    allocate(pml1_zt(abc,nx))
    allocate(pml1_zb(abc,nx))

    R=0.001 !Theoretical reflection parameter - free variable (~0.001)

    d0=(log(1.0/R))*(1.5*vp)/(abc*dx) !PML absorption parameter

    filter=0.0

    do ix=1,abc
      filter(ix)=(ix/(1.0*abc))*(ix/(1.0*abc))
    enddo

    !Filter in the x direction
    do ix=1,abc
      do iz=1,nz
      d0_xl(iz,ix) = d0(iz,ix)*filter(abc-ix+1)
      d0_xr(iz,ix) = d0(iz,ix)*filter(ix)
      enddo
    enddo
    pml0_xl = 1.0-(dt*d0_xl/2.0)
    pml0_xr = 1.0-(dt*d0_xr/2.0)
    pml1_xl = 1.0/(1.0+(dt*d0_xl/2.0))
    pml1_xr = 1.0/(1.0+(dt*d0_xr/2.0))
    !Filter in the z direction
    do ix=1,nx
      do iz=1,abc
      d0_zt(iz,ix) = d0(iz,ix)*filter(abc-iz+1)
      d0_zb(iz,ix) = d0(iz,ix)*filter(iz)
      enddo
    enddo
    pml0_zt = 1.0-(dt*d0_zt/2.0)
    pml0_zb = 1.0-(dt*d0_zb/2.0)
    pml1_zt = 1.0/(1.0+(dt*d0_zt/2.0))
    pml1_zb = 1.0/(1.0+(dt*d0_zb/2.0))
    if (free_surface) then
      pml0_zt = 1.0
      pml1_zt = 1.0
    endif

    !pml0_x = 1.0
    !pml0_z = 1.0
    !pml1_x = 1.0
    !pml1_z = 1.0

!    call to_history("n1",nz, "Dat/d0_xl.H")
!    call to_history("o1",0.0, "Dat/d0_xl.H")
!    call to_history("d1",1, "Dat/d0_xl.H")
!    call to_history("n2",nx, "Dat/d0_xl.H")
!    call to_history("o2",0.0, "Dat/d0_xl.H")
!    call to_history("d2",1, "Dat/d0_xl.H")
!    call sep_write(d0_xl,"Dat/d0_xl.H")

!    call to_history("n1",nz, "Dat/d0_xr.H")
!    call to_history("o1",0.0, "Dat/d0_xr.H")
!    call to_history("d1",1, "Dat/d0_xr.H")
!    call to_history("n2",nx, "Dat/d0_xr.H")
!    call to_history("o2",0.0, "Dat/d0_xr.H")
!    call to_history("d2",1, "Dat/d0_xr.H")
!    call sep_write(d0_xr,"Dat/d0_xr.H")

!    call to_history("n1",nz, "Dat/d0_zt.H")
!    call to_history("o1",0.0, "Dat/d0_zt.H")
!    call to_history("d1",1, "Dat/d0_zt.H")
!    call to_history("n2",nx, "Dat/d0_zt.H")
!    call to_history("o2",0.0, "Dat/d0_zt.H")
!    call to_history("d2",1, "Dat/d0_zt.H")
!    call sep_write(d0_zt,"Dat/d0_zt.H")

!    call to_history("n1",nz, "Dat/d0_zb.H")
!    call to_history("o1",0.0, "Dat/d0_zb.H")
!    call to_history("d1",1, "Dat/d0_zb.H")
!    call to_history("n2",nx, "Dat/d0_zb.H")
!    call to_history("o2",0.0, "Dat/d0_zb.H")
!    call to_history("d2",1, "Dat/d0_zb.H")
!    call sep_write(d0_zb,"Dat/d0_zb.H")

  endsubroutine

!============================================================
!FUNCTIONS
!============================================================
  function funK(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real*16,dimension(:)         :: model,data
    integer                      :: stat
    call K(adj,add,model,data)
    stat=0
  endfunction


!============================================================
!SUBROUTINES
!============================================================
  subroutine K(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real*16,dimension(nz,nx),intent(inout)    :: model,data

    if (adj)
    else
      if (.not.add) model=0.
      do ix=1,abc
        do iz=1,nz      
          model(iz,ix+5) = model(iz,ix+5) + pml1_xl(iz,ix)*data(iz,ix+5)
          model(iz,nx-abc-5+ix) = model(iz,nx-abc-5+ix) + pml1_xr(iz,ix)*data(iz,nx-abc-5+ix)
        enddo
      enddo
      do ix=abc+1,nx-abc
        do iz=1,nz
        enddo
      enddo
  endif
    do ix=1,abc
      do iz=1,nz      
      enddo
    enddo
  endsubroutine

  subroutine L(data)
    real,dimension(nz,nx)     :: data

    do ix=1,abc
      do iz=1,nz
        data(iz,ix+5) = pml0_xl(iz,ix)*data(iz,ix+5)
      enddo
    enddo
    do ix=1,abc
      do iz=1,nz
        data(iz,nx-abc-5+ix) = pml0_xr(iz,ix)*data(iz,nx-abc-5+ix)
      enddo
    enddo
  endsubroutine

  subroutine M(data)
    real,dimension(nz,nx)     :: data

    do ix=1,nx
      do iz=1,abc
        data(iz+5,ix) = pml1_zt(iz,ix)*data(iz+5,ix)
      enddo
    enddo
    do ix=1,nx
      do iz=1,abc
        data(nz-abc-5+iz,ix) = pml1_zb(iz,ix)*data(nz-abc-5+iz,ix)
      enddo
    enddo
  endsubroutine

  subroutine N(data)
    real,dimension(nz,nx)     :: data

    do ix=1,nx
      do iz=1,abc
        data(iz+5,ix) = pml0_zt(iz,ix)*data(iz+5,ix)
      enddo
    enddo
    do ix=1,nx
      do iz=1,abc
        data(nz-abc-5+iz,ix) = pml0_zb(iz,ix)*data(nz-abc-5+iz,ix)
      enddo
    enddo
  endsubroutine

  subroutine O(data)
    real,dimension(nz,nx)     :: data

    do ix=1,abc-1
      do iz=1,nz      
        data(iz,ix+5) = 0.5*(pml1_xl(iz,ix+1)+ &
                             pml1_xl(iz,ix  ))*data(iz,ix+5)
      enddo
    enddo
    do ix=1,abc-1
      do iz=1,nz      
        data(iz,nx-abc-5+ix) = 0.5*(pml1_xr(iz,ix+1)+ &
                                    pml1_xr(iz,ix  ))*data(iz,nx-abc-5+ix)
      enddo
    enddo
  endsubroutine

  subroutine P(data)
    real,dimension(nz,nx)     :: data

    do ix=1,abc-1
      do iz=1,nz      
        data(iz,ix+5) = 0.5*(pml0_xl(iz,ix+1)+ &
                             pml0_xl(iz,ix  ))*data(iz,ix+5)
      enddo
    enddo
    do ix=1,abc-1
      do iz=1,nz      
        data(iz,nx-abc-5+ix) = 0.5*(pml0_xr(iz,ix+1)+ &
                                    pml0_xr(iz,ix  ))*data(iz,nx-abc-5+ix)
      enddo
    enddo
  endsubroutine

  subroutine Q(data)
    real,dimension(nz,nx)     :: data

    do ix=1,nx
      do iz=2,abc      
        data(iz+5,ix) = 0.5*(pml1_zt(iz  ,ix)+ &
                             pml1_zt(iz-1,ix))*data(iz+5,ix)
      enddo
    enddo
    do ix=1,nx
      do iz=2,abc      
        data(nz-abc-5+iz,ix) = 0.5*(pml1_zb(iz  ,ix)+ &
                                    pml1_zb(iz-1,ix))*data(nz-abc-5+iz,ix)
      enddo
    enddo
  endsubroutine

  subroutine R(data)
    real,dimension(nz,nx)     :: data

    do ix=1,nx
      do iz=2,abc      
        data(iz+5,ix) = 0.5*(pml0_zt(iz  ,ix)+ &
                             pml0_zt(iz-1,ix))*data(iz+5,ix)
      enddo
    enddo
    do ix=1,nx
      do iz=2,abc      
        data(nz-abc-5+iz,ix) = 0.5*(pml0_zb(iz  ,ix)+ &
                                    pml0_zb(iz-1,ix))*data(nz-abc-5+iz,ix)
      enddo
    enddo
  endsubroutine

endmodule
