! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Mar, 24th, 2015

! finite difference module
module rtmlib_mod

  use sep
  use omp_lib

  implicit none
  integer, private   :: nx, nz, iz, ix
  real, private      :: dx, dt, aux
  real, private      :: c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0
  real, pointer, dimension(:,:),private :: m1,m1_x,m1_z,m2,m3
  real, allocatable, dimension(:,:), private :: aux_m2,aux_m3,aux_m2m3
  logical,private    :: lame

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine rtmmod_init(nz_in,nx_in,dx_in, dt_in, m1_in,m1_x_in,m1_z_in,m2_in,m3_in,lame_in)
    integer        :: nx_in, nz_in
    real           :: dx_in, dt_in
    real, dimension(:,:),target :: m1_in,m1_x_in, m1_z_in, m2_in, m3_in
    logical        :: lame_in
    !OMP settings
!    integer                               :: node, nodes

    !OMP
    !nodes=16
    !node = omp_get_num_procs()
    !call omp_set_num_threads(nodes)

    nx = nx_in
    nz = nz_in
    dx = dx_in
    dt = dt_in
    m1 => m1_in
    m1_x => m1_x_in
    m1_z => m1_z_in
    m2 => m2_in
    m3 => m3_in
    lame = lame_in

    allocate(aux_m2(nz,nx))
    allocate(aux_m3(nz,nx))
    allocate(aux_m2m3(nz,nx))

    aux = dt/dx

    m1_x = aux*m1_x
    m1_z = aux*m1_z
    if (lame) then
      aux_m2 = aux*m2
      aux_m3 = aux*m3
      aux_m2m3 = aux*(m2+2.0*m3)
    else
      aux_m2 = aux*m1*(m2*m2-2.0*m3*m3)
      aux_m3 = aux*m1*m3*m3
      aux_m2m3 = aux*m1*m2*m2
    endif

  endsubroutine

!=============================================================================
!=============================================================================
! Functions (for the adjoints, dot product test, etc)
  function funAx(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call Ax(adj,add,model,data)
    stat=0
  endfunction

  function funBx(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call Bx(adj,add,model,data)
    stat=0
  endfunction

  function funAz(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call Az(adj,add,model,data)
    stat=0
  endfunction

  function funBz(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call Bz(adj,add,model,data)
    stat=0
  endfunction

  function funCx(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call Cx(adj,add,model,data)
    stat=0
  endfunction

  function funCz(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call Cz(adj,add,model,data)
    stat=0
  endfunction


  function funD(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call D(adj,add,model,data)
    stat=0
  endfunction

  function funE(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call E(adj,add,model,data)
    stat=0
  endfunction

  function funF(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call F(adj,add,model,data)
    stat=0
  endfunction

  function funG(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call G(adj,add,model,data)
    stat=0
  endfunction

  function funH(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call H(adj,add,model,data)
    stat=0
  endfunction

  function funJ(adj,add,model,data) result(stat)
    logical,intent(in)           :: adj,add
    real,dimension(:)          :: model,data
    integer                      :: stat
    call J(adj,add,model,data)
    stat=0
  endfunction

!=============================================================================
!!!! 10TH ORDER OPERATOR
  subroutine Ax(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data
    integer :: ith

    if(adj)then
      if(.not.add) model = 0.
      do ix=6,nx-5
!!$OMP PARALLEL DO SHARED(model,data,m1_x,ix) PRIVATE(ix,iz)
!        write(0,*) omp_get_thread_num()
        do iz=5,nz-5
          model(iz,ix+4) = model(iz,ix+4) + m1_x(iz,ix)*c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + m1_x(iz,ix)*c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + m1_x(iz,ix)*c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + m1_x(iz,ix)*c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + m1_x(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - m1_x(iz,ix)*c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - m1_x(iz,ix)*c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - m1_x(iz,ix)*c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - m1_x(iz,ix)*c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - m1_x(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data,m1_x,c1,c2,c3,c4,c5,nx,nz) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
        data(iz,ix) = data(iz,ix) + m1_x(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
                                                 c2*(model(iz,ix+3) - model(iz,ix-4))+&
                                                 c3*(model(iz,ix+2) - model(iz,ix-3))+&
                                                 c4*(model(iz,ix+1) - model(iz,ix-2))+&
                                                 c5*(model(iz,ix  ) - model(iz,ix-1)))        
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Bx(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data,m1_x) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + m1_x(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + m1_x(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + m1_x(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + m1_x(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + m1_x(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - m1_x(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - m1_x(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - m1_x(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - m1_x(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - m1_x(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data,m1_x) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + m1_x(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                   c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                   c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                   c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                   c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine Az(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          model(iz,ix+5) = model(iz,ix+5) + m1_z(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + m1_z(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + m1_z(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + m1_z(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + m1_z(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - m1_z(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - m1_z(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - m1_z(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - m1_z(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - m1_z(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + m1_z(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                   c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                   c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                   c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                   c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Bz(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          model(iz+4,ix) = model(iz+4,ix) + m1_z(iz,ix)*c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + m1_z(iz,ix)*c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + m1_z(iz,ix)*c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + m1_z(iz,ix)*c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + m1_z(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - m1_z(iz,ix)*c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - m1_z(iz,ix)*c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - m1_z(iz,ix)*c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - m1_z(iz,ix)*c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - m1_z(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + m1_z(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                                   c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                                   c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                                   c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                                   c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine Cx(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          model(iz,ix) = model(iz,ix) + dt*m1_x(iz,ix)*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          data(iz,ix) = data(iz,ix) + dt*m1_x(iz,ix)*model(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Cz(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          model(iz,ix) = model(iz,ix) + dt*m1_z(iz,ix)*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          data(iz,ix) = data(iz,ix) + dt*m1_z(iz,ix)*model(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine D(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz,ix+5) = model(iz,ix+5) + aux_m2m3(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + aux_m2m3(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m2m3(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m2m3(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m2m3(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - aux_m2m3(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m2m3(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m2m3(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m2m3(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m2m3(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + aux_m2m3(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                       c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                       c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                       c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                       c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine E(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + aux_m2(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + aux_m2(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m2(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m2(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m2(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - aux_m2(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m2(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m2(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m2(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m2(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + aux_m2(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                     c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                     c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                     c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                     c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine F(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz,ix+5) = model(iz,ix+5) + aux_m2(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + aux_m2(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m2(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m2(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m2(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - aux_m2(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m2(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m2(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m2(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m2(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + aux_m2(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                     c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                     c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                     c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                     c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine G(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + aux_m2m3(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + aux_m2m3(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m2m3(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m2m3(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m2m3(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - aux_m2m3(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m2m3(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m2m3(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m2m3(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m2m3(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + aux_m2m3(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                       c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                       c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                       c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                       c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine H(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          model(iz+4,ix) = model(iz+4,ix) + aux_m3(iz,ix)*c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + aux_m3(iz,ix)*c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + aux_m3(iz,ix)*c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + aux_m3(iz,ix)*c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + aux_m3(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - aux_m3(iz,ix)*c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - aux_m3(iz,ix)*c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - aux_m3(iz,ix)*c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - aux_m3(iz,ix)*c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - aux_m3(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + aux_m3(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                                     c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                                     c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                                     c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                                     c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine J(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj)then
      if(.not.add) model = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          model(iz,ix+4) = model(iz,ix+4) + aux_m3(iz,ix)*c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + aux_m3(iz,ix)*c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + aux_m3(iz,ix)*c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + aux_m3(iz,ix)*c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + aux_m3(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - aux_m3(iz,ix)*c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - aux_m3(iz,ix)*c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - aux_m3(iz,ix)*c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - aux_m3(iz,ix)*c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - aux_m3(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
        data(iz,ix) = data(iz,ix) + aux_m3(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
                                                   c2*(model(iz,ix+3) - model(iz,ix-4))+&
                                                   c3*(model(iz,ix+2) - model(iz,ix-3))+&
                                                   c4*(model(iz,ix+1) - model(iz,ix-2))+&
                                                   c5*(model(iz,ix  ) - model(iz,ix-1)))
        enddo
      enddo
!!!$OMP END PARALLEL DO
    endif
  endsubroutine

end module
