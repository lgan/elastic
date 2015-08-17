! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! Mar, 24th, 2015

! finite difference module
module bornlib_mod

  use sep

  implicit none
  integer, private   :: nx, nz, iz, ix
  real, private      :: c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0
  real, pointer, dimension(:,:),private :: rho_Vx,rho_Vz,l,mi,l2mi
  real, pointer, dimension(:,:),private :: drho_Vx,drho_vz,dl,dmi,dl2mi

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine bornmod_init(nz_in,nx_in,rho_Vx_in,rho_Vz_in,l_in,mi_in,l2mi_in,drho_Vx_in,drho_Vz_in,dl_in,dmi_in,dl2mi_in)
    integer        :: nx_in, nz_in
    real, dimension(:,:),target :: rho_Vx_in, rho_Vz_in, l_in, mi_in, l2mi_in
    real, dimension(:,:),target :: drho_Vx_in, drho_Vz_in, dl_in, dmi_in, dl2mi_in

    nx = nx_in
    nz = nz_in
    rho_Vx => rho_Vx_in
    rho_Vz => rho_vz_in
    l => l_in
    mi => mi_in
    l2mi => l2mi_in
    drho_Vx => drho_Vx_in
    drho_Vz => drho_vz_in
    dl => dl_in
    dmi => dmi_in
    dl2mi => dl2mi_in

  endsubroutine

!=============================================================================
!=============================================================================
!=============================================================================
!!!! 10TH ORDER OPERATOR
  subroutine Ax(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj)then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-4
        do iz=1,nz
          model(iz,ix+4) = model(iz,ix+4) + c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-4
        do iz=6,nz
        data(iz,ix) = data(iz,ix) + rho_Vx(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
                                                   c2*(model(iz,ix+3) - model(iz,ix-4))+&
                                                   c3*(model(iz,ix+2) - model(iz,ix-3))+&
                                                   c4*(model(iz,ix+1) - model(iz,ix-2))+&
                                                   c5*(model(iz,ix  ) - model(iz,ix-1)))        
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Bx(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + rho_Vx(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                     c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                     c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                     c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                     c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine Az(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=1,nz
          model(iz,ix+5) = model(iz,ix+5) + c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz
          data(iz,ix) = data(iz,ix) + rho_Vz(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                                     c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                                     c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                                     c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                                     c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine Bz(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-4
          model(iz+4,ix) = model(iz+4,ix) + c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-4
          data(iz,ix) = data(iz,ix) + rho_Vz(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                                     c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                                     c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                                     c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                                     c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine D(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=1,nz
          model(iz,ix+5) = model(iz,ix+5) + c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz
          data(iz,ix) = data(iz,ix) + l2mi(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                      c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                      c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                      c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                      c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine E(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + l(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                      c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                      c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                      c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                      c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine F(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=1,nz
          model(iz,ix+5) = model(iz,ix+5) + c1*data(iz,ix)
          model(iz,ix+4) = model(iz,ix+4) + c2*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + c3*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + c4*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + c5*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) - c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - c4*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - c3*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - c2*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz
          data(iz,ix) = data(iz,ix) + l(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                      c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                      c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                      c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                      c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine G(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=5,nz-5
          model(iz+5,ix) = model(iz+5,ix) + c1*data(iz,ix)
          model(iz+4,ix) = model(iz+4,ix) + c2*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + c3*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + c4*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + c5*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) - c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - c4*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - c3*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - c2*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + l2mi(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                      c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                      c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                      c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                      c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine H(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-4
          model(iz+4,ix) = model(iz+4,ix) + c1*data(iz,ix)
          model(iz+3,ix) = model(iz+3,ix) + c2*data(iz,ix)
          model(iz+2,ix) = model(iz+2,ix) + c3*data(iz,ix)
          model(iz+1,ix) = model(iz+1,ix) + c4*data(iz,ix)
          model(iz  ,ix) = model(iz  ,ix) + c5*data(iz,ix)
          model(iz-1,ix) = model(iz-1,ix) - c5*data(iz,ix)
          model(iz-2,ix) = model(iz-2,ix) - c4*data(iz,ix)
          model(iz-3,ix) = model(iz-3,ix) - c3*data(iz,ix)
          model(iz-4,ix) = model(iz-4,ix) - c2*data(iz,ix)
          model(iz-5,ix) = model(iz-5,ix) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-4
          data(iz,ix) = data(iz,ix) + mi(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                      c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                      c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                      c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                      c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine J(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj)then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-4
        do iz=1,nz
          model(iz,ix+4) = model(iz,ix+4) + c1*data(iz,ix)
          model(iz,ix+3) = model(iz,ix+3) + c2*data(iz,ix)
          model(iz,ix+2) = model(iz,ix+2) + c3*data(iz,ix)
          model(iz,ix+1) = model(iz,ix+1) + c4*data(iz,ix)
          model(iz,ix  ) = model(iz,ix  ) + c5*data(iz,ix)
          model(iz,ix-1) = model(iz,ix-1) - c5*data(iz,ix)
          model(iz,ix-2) = model(iz,ix-2) - c4*data(iz,ix)
          model(iz,ix-3) = model(iz,ix-3) - c3*data(iz,ix)
          model(iz,ix-4) = model(iz,ix-4) - c2*data(iz,ix)
          model(iz,ix-5) = model(iz,ix-5) - c1*data(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-4
        do iz=6,nz
        data(iz,ix) = data(iz,ix) + mi(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
                                               c2*(model(iz,ix+3) - model(iz,ix-4))+&
                                               c3*(model(iz,ix+2) - model(iz,ix-3))+&
                                               c4*(model(iz,ix+1) - model(iz,ix-2))+&
                                               c5*(model(iz,ix  ) - model(iz,ix-1)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

!=======================================================
! Born source scaling terms
!=======================================================

  subroutine S1(adj,add,model1,model2,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model1,model2,data

    if(adj) then
      if(.not.add) then
        model1 = 0.
        model2 = 0.
      endif
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model1,model2,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + (model1(iz,ix)-model2(iz,ix))*drho_Vx(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine S2(adj,add,model1,model2,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model1,model2,data

    if(adj) then
      if(.not.add) then
        model1 = 0.
        model2 = 0.
      endif
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model1,model2,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + (model1(iz,ix)-model2(iz,ix))*drho_Vz(iz,ix)
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine S3(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout) :: model,data

    if(adj) then
      if(.not.add) model = 0.
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz
          data(iz,ix) = data(iz,ix) + dl2mi(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                      c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                      c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                      c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                      c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine S4(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + dl(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                                 c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                                 c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                                 c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                                 c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine S5(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=5,nx-5
        do iz=6,nz
          data(iz,ix) = data(iz,ix) + dl(iz,ix)*(c1*(model(iz,ix+5) - model(iz,ix-4))+&
                                      c2*(model(iz,ix+4) - model(iz,ix-3))+&
                                      c3*(model(iz,ix+3) - model(iz,ix-2))+&
                                      c4*(model(iz,ix+2) - model(iz,ix-1))+&
                                      c5*(model(iz,ix+1) - model(iz,ix  )))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine S6(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-5
          data(iz,ix) = data(iz,ix) + dl2mi(iz,ix)*(c1*(model(iz+5,ix) - model(iz-4,ix))+&
                                      c2*(model(iz+4,ix) - model(iz-3,ix))+&
                                      c3*(model(iz+3,ix) - model(iz-2,ix))+&
                                      c4*(model(iz+2,ix) - model(iz-1,ix))+&
                                      c5*(model(iz+1,ix) - model(iz  ,ix)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine S7(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj) then
      if(.not.add) model = 0.
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=1,nx
        do iz=6,nz-4
          data(iz,ix) = data(iz,ix) + dmi(iz,ix)*(c1*(model(iz+4,ix) - model(iz-5,ix))+&
                                      c2*(model(iz+3,ix) - model(iz-4,ix))+&
                                      c3*(model(iz+2,ix) - model(iz-3,ix))+&
                                      c4*(model(iz+1,ix) - model(iz-2,ix))+&
                                      c5*(model(iz  ,ix) - model(iz-1,ix)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 
  endsubroutine

  subroutine S8(adj,add,model,data)
    logical,intent(in)          :: adj,add
    !real*16,dimension(nz,nx)    :: model,data
    real,dimension(nz,nx),intent(inout)    :: model,data

    if(adj)then
      if(.not.add) model = 0.
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data) PRIVATE(ix,iz)
      do ix=6,nx-4
        do iz=6,nz
        data(iz,ix) = data(iz,ix) + dmi(iz,ix)*(c1*(model(iz,ix+4) - model(iz,ix-5))+&
                                               c2*(model(iz,ix+3) - model(iz,ix-4))+&
                                               c3*(model(iz,ix+2) - model(iz,ix-3))+&
                                               c4*(model(iz,ix+1) - model(iz,ix-2))+&
                                               c5*(model(iz,ix  ) - model(iz,ix-1)))
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

end module
