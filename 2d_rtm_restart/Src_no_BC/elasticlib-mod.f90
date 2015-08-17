! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Mar, 24th, 2015

! finite difference module
module elasticlib_mod

  use sep
  use omp_lib

  implicit none
  integer, private   :: nx, nz, iz, ix
  real, private      :: dx, dt
  real, private    :: c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0
  real, pointer, dimension(:,:),private :: rho_Vx,rho_Vz,l,mi,l2mi

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine elasticmod_init(nz_in,nx_in,dx_in, dt_in, rho_Vx_in,rho_Vz_in,l_in,mi_in,l2mi_in)
    integer        :: nx_in, nz_in
    real           :: dx_in, dt_in
    real, dimension(:,:),target :: rho_Vx_in, rho_Vz_in, l_in, mi_in, l2mi_in

    nx = nx_in
    nz = nz_in
    dx = dx_in
    dt = dt_in
    rho_Vx => rho_Vx_in
    rho_Vz => rho_vz_in
    l => l_in
    mi => mi_in
    l2mi => l2mi_in

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

    if(adj)then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
          model(iz,ix+4) = model(iz,ix+4) + rho_Vx(iz,ix)*c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + rho_Vx(iz,ix)*c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + rho_Vx(iz,ix)*c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + rho_Vx(iz,ix)*c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + rho_Vx(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - rho_Vx(iz,ix)*c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - rho_Vx(iz,ix)*c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - rho_Vx(iz,ix)*c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - rho_Vx(iz,ix)*c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - rho_Vx(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
        data(iz,ix) = data(iz,ix) + rho_Vx(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
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
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + rho_Vx(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + rho_Vx(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + rho_Vx(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + rho_Vx(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + rho_Vx(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - rho_Vx(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - rho_Vx(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - rho_Vx(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - rho_Vx(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - rho_Vx(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + rho_Vx(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                     c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                     c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                     c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                     c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine Az(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          model(iz,ix+5) = model(iz,ix+5) + rho_Vz(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + rho_Vz(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + rho_Vz(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + rho_Vz(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + rho_Vz(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - rho_Vz(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - rho_Vz(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - rho_Vz(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - rho_Vz(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - rho_Vz(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + rho_Vz(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                     c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                     c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                     c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                     c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Bz(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          model(iz+4,ix) = model(iz+4,ix) + rho_Vz(iz,ix)*c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + rho_Vz(iz,ix)*c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + rho_Vz(iz,ix)*c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + rho_Vz(iz,ix)*c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + rho_Vz(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - rho_Vz(iz,ix)*c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - rho_Vz(iz,ix)*c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - rho_Vz(iz,ix)*c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - rho_Vz(iz,ix)*c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - rho_Vz(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + rho_Vz(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                                     c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                                     c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                                     c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                                     c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine Cx(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          model(iz,ix) = model(iz,ix) + dt*rho_Vx(iz,ix)*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          data(iz,ix) = data(iz,ix) + dt*rho_Vx(iz,ix)*model(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Cz(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          model(iz,ix) = model(iz,ix) + dt*rho_Vz(iz,ix)*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=1,nz
          data(iz,ix) = data(iz,ix) + dt*rho_Vz(iz,ix)*model(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine D(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz,ix+5) = model(iz,ix+5) + l2mi(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + l2mi(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + l2mi(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + l2mi(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + l2mi(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - l2mi(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - l2mi(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - l2mi(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - l2mi(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - l2mi(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + l2mi(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                   c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                   c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                   c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                   c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine E(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + l(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + l(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + l(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + l(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + l(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - l(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - l(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - l(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - l(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - l(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + l(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine F(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz,ix+5) = model(iz,ix+5) + l(iz,ix)*c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + l(iz,ix)*c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + l(iz,ix)*c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + l(iz,ix)*c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + l(iz,ix)*c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - l(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - l(iz,ix)*c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - l(iz,ix)*c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - l(iz,ix)*c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - l(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + l(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine G(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + l2mi(iz,ix)*c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + l2mi(iz,ix)*c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + l2mi(iz,ix)*c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + l2mi(iz,ix)*c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + l2mi(iz,ix)*c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - l2mi(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - l2mi(iz,ix)*c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - l2mi(iz,ix)*c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - l2mi(iz,ix)*c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - l2mi(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=5,nz-5
          data(iz,ix) = data(iz,ix) + l2mi(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                   c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                   c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                   c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                   c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine H(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          model(iz+4,ix) = model(iz+4,ix) + mi(iz,ix)*c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + mi(iz,ix)*c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + mi(iz,ix)*c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + mi(iz,ix)*c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + mi(iz,ix)*c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - mi(iz,ix)*c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - mi(iz,ix)*c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - mi(iz,ix)*c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - mi(iz,ix)*c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - mi(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + mi(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                                 c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                                 c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                                 c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                                 c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine J(adj,add,model,data)
    logical,intent(in)          :: adj,add
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj)then
      if(.not.add) model = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
          model(iz,ix+4) = model(iz,ix+4) + mi(iz,ix)*c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + mi(iz,ix)*c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + mi(iz,ix)*c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + mi(iz,ix)*c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + mi(iz,ix)*c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - mi(iz,ix)*c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - mi(iz,ix)*c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - mi(iz,ix)*c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - mi(iz,ix)*c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - mi(iz,ix)*c1*data(iz,ix)
        enddo
      enddo
!!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-5
        do iz=6,nz-5
        data(iz,ix) = data(iz,ix) + mi(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
                                               c2*(model(iz,ix+3) - model(iz,ix-4))+&
                                               c3*(model(iz,ix+2) - model(iz,ix-3))+&
                                               c4*(model(iz,ix+1) - model(iz,ix-2))+&
                                               c5*(model(iz,ix  ) - model(iz,ix-1)))
        enddo
      enddo
!!$OMP END PARALLEL DO
    endif
  endsubroutine

end module
