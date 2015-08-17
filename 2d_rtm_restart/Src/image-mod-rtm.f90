! Program created by Gustavo Catao Alves at
! the Stanford Exploration Project
! May, 6th, 2014

! seismogram module
module image_mod_rtm

  use sep
  use omp_lib

  implicit none
  integer, private :: nx,nz,nt,ix,iz,it
  real, private    :: dx,dz,dt
  real, private    :: c1=35.0/294912.0,c2=-405.0/229376.0,c3=567.0/40960.0,c4=-735.0/8192.0,c5=19845.0/16384.0
  real, dimension(:,:,:), pointer, private :: Vx0,Vz0

  !Added later
  real, dimension(:,:), pointer, private   :: m1,m2,m3
  logical, private  :: lame

  real, allocatable, dimension(:,:,:), private :: S1,S2,S3,S4,S5,S6

  contains

  ! Initialization of all variables needed by this subroutine
  subroutine image_init(nx_in,dx_in,nz_in,dz_in,nt_in,dt_in,Vx0_in,Vz0_in,lame_in)
    integer,intent(in)              :: nx_in,nz_in,nt_in
    real,intent(in)                 :: dx_in,dz_in,dt_in
    real,dimension(:,:,:),target    :: Vx0_in,Vz0_in
    logical                         :: lame_in

    nx = nx_in
    dx = dx_in
    nz = nz_in
    dz = dz_in
    nt = nt_in
    dt = dt_in
    lame = lame_in

    Vx0 => Vx0_in
    Vz0 => Vz0_in

    allocate(S1(nz,nx,nt))
    allocate(S2(nz,nx,nt))
    allocate(S3(nz,nx,nt))
    allocate(S4(nz,nx,nt))
    allocate(S5(nz,nx,nt))
    allocate(S6(nz,nx,nt))

    S1 = 0.
    S2 = 0.
    S3 = 0.
    S4 = 0.
    S5 = 0.
    S6 = 0.

    call T(.false.,.false.,Vx0,S1)
    call T(.false.,.false.,Vz0,S2)
    call vxdx(.false.,.false.,Vx0,S4)
    call vzdz(.false.,.false.,Vz0,S5)
    S3=S4+S5
    !S4=2.0*S4
    !S5=2.0*S5
    call vxdz(.false.,.false.,Vx0,S6)
    call vzdx(.false.,.true.,Vz0,S6)
    
  endsubroutine

  subroutine image_init2(m1_in,m2_in,m3_in)
    real,dimension(:,:),target  :: m1_in,m2_in,m3_in

    m1 => m1_in
    m2 => m2_in
    m3 => m3_in

  endsubroutine

!==================================
! FUNCTION CALLS
!==================================

!  function image_op(adj,add,drho,dl,dmu,Vx,Vz,sigmaxx,sigmazz,sigmaxz) result(stat)
  function image_op(adj,add,dm1,dm2,dm3,Vx,Vz,sigmaxx,sigmazz,sigmaxz) result(stat)
    logical,intent(in)        :: adj,add
    !real,dimension(:)         :: drho,dl,dmu,Vx,Vz,sigmaxx,sigmazz,sigmaxz
    real,dimension(:)         :: dm1,dm2,dm3,Vx,Vz,sigmaxx,sigmazz,sigmaxz
    integer                   :: stat

    call image_op2(adj,add,dm1,dm2,dm3,Vx,Vz,sigmaxx,sigmazz,sigmaxz)
!    call image_op2(adj,add,drho,dl,dmu,Vx,Vz,sigmaxx,sigmazz,sigmaxz)
    stat=0
  endfunction

  function funT(adj,add,model,data) result(stat)
    logical,intent(in)        :: adj,add
    real,dimension(:)         :: model,data
    integer                   :: stat
    call T(adj,add,model,data)
    stat=0
  endfunction

  function funvxdx(adj,add,model,data) result(stat)
    logical,intent(in)        :: adj,add
    real,dimension(:)         :: model,data
    integer                   :: stat
    call vxdx(adj,add,model,data)
    stat=0
  endfunction

  function funvzdz(adj,add,model,data) result(stat)
    logical,intent(in)        :: adj,add
    real,dimension(:)         :: model,data
    integer                   :: stat
    call vzdz(adj,add,model,data)
    stat=0
  endfunction

  function funvxdz(adj,add,model,data) result(stat)
    logical,intent(in)        :: adj,add
    real,dimension(:)         :: model,data
    integer                   :: stat
    call vxdz(adj,add,model,data)
    stat=0
  endfunction

  function funvzdx(adj,add,model,data) result(stat)
    logical,intent(in)        :: adj,add
    real,dimension(:)         :: model,data
    integer                   :: stat
    call vzdx(adj,add,model,data)
    stat=0
  endfunction

!==================================
! SUBROUTINES
!==================================
!   subroutine image_op2(adj,add,drho,dl,dmu,Vx,Vz,sigmaxx,sigmazz,sigmaxz)
   subroutine image_op2(adj,add,dm1,dm2,dm3,Vx,Vz,sigmaxx,sigmazz,sigmaxz)
    logical,intent(in)                     :: adj,add
!    real,dimension(nz,nx),intent(inout)    :: drho,dl,dmu
    real,dimension(nz,nx),intent(inout)    :: dm1,dm2,dm3
    real,dimension(nz,nx,nt),intent(inout) :: Vx,Vz,sigmaxx,sigmazz,sigmaxz

! IMAGING CONDITION FOR LAME PARAMETERS
    if (lame) then
      if(adj)then
        write(0,*)maxval(S1),maxval(S2),maxval(S3),maxval(S4),maxval(S5),maxval(S6)
        write(0,*)maxval(Vx),maxval(Vz),maxval(sigmaxx),maxval(sigmazz),maxval(sigmaxz)
        if(.not.add) then
          dm1 = 0.
          dm2 = 0.
          dm3 = 0.
        endif
        do it=1,nt
          dm1(:,:) = dm1(:,:) + (Vx(:,:,it)*S1(:,:,it)) &
                              + (Vz(:,:,it)*S2(:,:,it))
          dm2(:,:) = dm2(:,:) - (sigmaxx(:,:,it)*S3(:,:,it)) &
                              - (sigmazz(:,:,it)*S3(:,:,it))
          dm3(:,:) = dm3(:,:) - (sigmaxx(:,:,it)*2.0*S4(:,:,it)) &
                              - (sigmazz(:,:,it)*2.0*S5(:,:,it)) &
                              - (sigmaxz(:,:,it)*S6(:,:,it))
        enddo
      else
        if(.not.add) then
          Vx = 0.
          Vz = 0.
          sigmaxx = 0.
          sigmazz = 0.
          sigmaxz = 0.
        endif
        do it=1,nt
          Vx(:,:,it) = Vx(:,:,it) + (dm1(:,:)*S1(:,:,it))
          Vz(:,:,it) = Vz(:,:,it) + (dm1(:,:)*S2(:,:,it))
          sigmaxx(:,:,it) = sigmaxx(:,:,it) - (dm2(:,:)*S3(:,:,it)) &
                                            - (dm3(:,:)*2.0*S4(:,:,it))
          sigmazz(:,:,it) = sigmazz(:,:,it) - (dm2(:,:)*S3(:,:,it)) &
                                            - (dm3(:,:)*2.0*S5(:,:,it))
          sigmaxz(:,:,it) = sigmaxz(:,:,it) - (dm3(:,:)*S6(:,:,it))
        enddo
      endif
! IMAGING CONDITION FOR VELOCITY PARAMETERS
    else
      if(adj)then
        if(.not.add) then
          dm1 = 0.
          dm2 = 0.
          dm3 = 0.
        endif
        do it=1,nt
          dm1(:,:) = dm1(:,:) + (Vx(:,:,it)*S1(:,:,it)) &
                              + (Vz(:,:,it)*S2(:,:,it)) &
                              - (sigmaxx(:,:,it)*(m2(:,:)*m2(:,:)*S4(:,:,it) &
                                                 +(m2(:,:)*m2(:,:)-2.0*m3(:,:)*m3(:,:))*S5(:,:,it))) &
                              - (sigmazz(:,:,it)*(m2(:,:)*m2(:,:)*S5(:,:,it) &
                                                 +(m2(:,:)*m2(:,:)-2.0*m3(:,:)*m3(:,:))*S4(:,:,it))) &
                              - (sigmaxz(:,:,it)*(m3(:,:)*m3(:,:)*S6(:,:,it)))
          dm2(:,:) = dm2(:,:) - (sigmaxx(:,:,it)*(2.0*m1(:,:)*m2(:,:)*S3(:,:,it))) &
                              - (sigmazz(:,:,it)*(2.0*m1(:,:)*m2(:,:)*S3(:,:,it)))
          dm3(:,:) = dm3(:,:) + (sigmaxx(:,:,it)*(4.0*m1(:,:)*m3(:,:)*S5(:,:,it))) &
                              + (sigmazz(:,:,it)*(4.0*m1(:,:)*m3(:,:)*S4(:,:,it))) &
                              - (sigmaxz(:,:,it)*(2.0*m1(:,:)*m3(:,:)*S6(:,:,it)))
        enddo
      else
        if(.not.add) then
          Vx = 0.
          Vz = 0.
          sigmaxx = 0.
          sigmazz = 0.
          sigmaxz = 0.
        endif
        do it=1,nt
          Vx(:,:,it) = Vx(:,:,it) + (dm1(:,:)*S1(:,:,it))
          Vz(:,:,it) = Vz(:,:,it) + (dm1(:,:)*S2(:,:,it))
          sigmaxx(:,:,it) = sigmaxx(:,:,it) - (dm1(:,:)*(m2(:,:)*m2(:,:)*S4(:,:,it) &
                                                        +(m2(:,:)*m2(:,:)-2.0*m3(:,:)*m3(:,:))*S5(:,:,it))) &
                                            - (dm2(:,:)*(2.0*m1(:,:)*m2(:,:)*S3(:,:,it))) &
                                            + (dm3(:,:)*(4.0*m1(:,:)*m3(:,:)*S5(:,:,it)))
          sigmazz(:,:,it) = sigmazz(:,:,it) - (dm1(:,:)*(m2(:,:)*m2(:,:)*S5(:,:,it) &
                                                        +(m2(:,:)*m2(:,:)-2.0*m3(:,:)*m3(:,:))*S4(:,:,it))) &
                                            - (dm2(:,:)*(2.0*m1(:,:)*m2(:,:)*S3(:,:,it))) &
                                            + (dm3(:,:)*(4.0*m1(:,:)*m3(:,:)*S4(:,:,it)))
          sigmaxz(:,:,it) = sigmaxz(:,:,it) - (dm1(:,:)*(m3(:,:)*m3(:,:)*S6(:,:,it))) &
                                            - (dm3(:,:)*(2.0*m1(:,:)*m3(:,:)*S6(:,:,it)))
        enddo
      endif
    endif


!    if(adj)then
!      if(.not.add) then
!        drho = 0.
!        dl = 0.
!        dmu = 0.
!      endif
!      do it=1,nt
!        drho(:,:)= drho(:,:) + S1(:,:,it)*Vx(:,:,it) + S2(:,:,it)*Vz(:,:,it)
!        dl(:,:)  = dl(:,:)   + S3(:,:,it)*(sigmaxx(:,:,it)+sigmazz(:,:,it))
!        dmu(:,:) = dmu(:,:)  + S4(:,:,it)*sigmaxx(:,:,it) +&
!                               S5(:,:,it)*sigmazz(:,:,it) +&
!                               S6(:,:,it)*sigmaxz(:,:,it)
!      enddo
!    else
!      if(.not.add) then
!        Vx = 0.
!        Vz = 0.
!        sigmaxx = 0.
!        sigmazz = 0.
!        sigmaxz = 0.
!      endif
!      do it=1,nt
!        Vx(:,:,it)=Vx(:,:,it)+S1(:,:,it)*drho(:,:)
!        Vz(:,:,it)=Vz(:,:,it)+S2(:,:,it)*drho(:,:)
!        sigmaxx(:,:,it)=sigmaxx(:,:,it)+S3(:,:,it)*dl(:,:)+S4(:,:,it)*dmu(:,:)
!        sigmazz(:,:,it)=sigmazz(:,:,it)+S3(:,:,it)*dl(:,:)+S5(:,:,it)*dmu(:,:)
!        sigmaxz(:,:,it)=sigmaxz(:,:,it)+S6(:,:,it)*dmu(:,:)
!      enddo
!    endif

!    deallocate(S1)
!    deallocate(S2)
!    deallocate(S3)
!    deallocate(S4)
!    deallocate(S5)
!    deallocate(S6)    

  endsubroutine

  subroutine T(adj,add,model,data)
    real                                   :: aux
    logical,intent(in)                     :: adj,add
    real,dimension(nz,nx,nt),intent(inout) :: model,data

    aux=1.0/dt

    if(adj) then
      if(.not.add) model = 0.
      do it=1,nt-1
!$OMP PARALLEL DO SHARED(model,data,aux,it) PRIVATE(ix,iz)
        do ix=1,nx
          do iz=1,nz
            model(iz,ix,(it+1)) = model(iz,ix,(it+1)) + aux*data(iz,ix,it)
            model(iz,ix,it) = model(iz,ix,it) - aux*data(iz,ix,it)
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
    else
      if(.not.add) data = 0.
      do it=1,nt-1
!$OMP PARALLEL DO SHARED(model,data,aux,it) PRIVATE(ix,iz)
        do ix=1,nx
          do iz=1,nz
            data(iz,ix,it) = data(iz,ix,it) + aux*(model(iz,ix,(it+1))-model(iz,ix,it))
          enddo
        enddo
!$OMP END PARALLEL DO
      enddo
    endif 
  endsubroutine

  subroutine vxdx(adj,add,model,data)
    real                                   :: aux
    logical,intent(in)                     :: adj,add
    real,dimension(nz,nx,nt),intent(inout) :: model,data

    aux=1.0/dx

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data,aux) PRIVATE(ix,iz,it)
      do it=1,nt
        do ix=5,nx-5
          do iz=5,nz-5
            model(iz,ix+5,it) = model(iz,ix+5,it) + c1*data(iz,ix,it)*aux
            model(iz,ix+4,it) = model(iz,ix+4,it) + c2*data(iz,ix,it)*aux
            model(iz,ix+3,it) = model(iz,ix+3,it) + c3*data(iz,ix,it)*aux
            model(iz,ix+2,it) = model(iz,ix+2,it) + c4*data(iz,ix,it)*aux
            model(iz,ix+1,it) = model(iz,ix+1,it) + c5*data(iz,ix,it)*aux
            model(iz,ix  ,it) = model(iz,ix  ,it) - c5*data(iz,ix,it)*aux
            model(iz,ix-1,it) = model(iz,ix-1,it) - c4*data(iz,ix,it)*aux
            model(iz,ix-2,it) = model(iz,ix-2,it) - c3*data(iz,ix,it)*aux
            model(iz,ix-3,it) = model(iz,ix-3,it) - c2*data(iz,ix,it)*aux
            model(iz,ix-4,it) = model(iz,ix-4,it) - c1*data(iz,ix,it)*aux
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data,aux) PRIVATE(ix,iz,it)
      do it=1,nt
        do ix=5,nx-5
          do iz=5,nz-5
            data(iz,ix,it) = data(iz,ix,it) + (c1*(model(iz,ix+5,it) - model(iz,ix-4,it))+&
                                               c2*(model(iz,ix+4,it) - model(iz,ix-3,it))+&
                                               c3*(model(iz,ix+3,it) - model(iz,ix-2,it))+&
                                               c4*(model(iz,ix+2,it) - model(iz,ix-1,it))+&
                                               c5*(model(iz,ix+1,it) - model(iz,ix  ,it)))*aux
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine vzdz(adj,add,model,data)
    real                                   :: aux
    logical,intent(in)                     :: adj,add
    real,dimension(nz,nx,nt),intent(inout) :: model,data

    aux=1.0/dx

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data,aux) PRIVATE(ix,iz,it)
      do it=1,nt
        do ix=5,nx-5
          do iz=5,nz-5
            model(iz+5,ix,it) = model(iz+5,ix,it) + c1*data(iz,ix,it)*aux
            model(iz+4,ix,it) = model(iz+4,ix,it) + c2*data(iz,ix,it)*aux
            model(iz+3,ix,it) = model(iz+3,ix,it) + c3*data(iz,ix,it)*aux
            model(iz+2,ix,it) = model(iz+2,ix,it) + c4*data(iz,ix,it)*aux
            model(iz+1,ix,it) = model(iz+1,ix,it) + c5*data(iz,ix,it)*aux
            model(iz  ,ix,it) = model(iz  ,ix,it) - c5*data(iz,ix,it)*aux
            model(iz-1,ix,it) = model(iz-1,ix,it) - c4*data(iz,ix,it)*aux
            model(iz-2,ix,it) = model(iz-2,ix,it) - c3*data(iz,ix,it)*aux
            model(iz-3,ix,it) = model(iz-3,ix,it) - c2*data(iz,ix,it)*aux
            model(iz-4,ix,it) = model(iz-4,ix,it) - c1*data(iz,ix,it)*aux
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data,aux) PRIVATE(ix,iz,it)
      do it=1,nt
        do ix=5,nx-5
          do iz=5,nz-5
            data(iz,ix,it) = data(iz,ix,it) + (c1*(model(iz+5,ix,it) - model(iz-4,ix,it))+&
                                               c2*(model(iz+4,ix,it) - model(iz-3,ix,it))+&
                                               c3*(model(iz+3,ix,it) - model(iz-2,ix,it))+&
                                               c4*(model(iz+2,ix,it) - model(iz-1,ix,it))+&
                                               c5*(model(iz+1,ix,it) - model(iz  ,ix,it)))*aux
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
    endif 

  endsubroutine

  subroutine vxdz(adj,add,model,data)
    real                                   :: aux
    logical,intent(in)                     :: adj,add
    real,dimension(nz,nx,nt),intent(inout) :: model,data

    aux=1.0/dx

    if(adj) then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data,aux) PRIVATE(ix,iz,it)
      do it=1,nt
        do ix=6,nx-5
          do iz=6,nz-5
            model(iz+4,ix,it) = model(iz+4,ix,it) + c1*data(iz,ix,it)*aux
            model(iz+3,ix,it) = model(iz+3,ix,it) + c2*data(iz,ix,it)*aux
            model(iz+2,ix,it) = model(iz+2,ix,it) + c3*data(iz,ix,it)*aux
            model(iz+1,ix,it) = model(iz+1,ix,it) + c4*data(iz,ix,it)*aux
            model(iz  ,ix,it) = model(iz  ,ix,it) + c5*data(iz,ix,it)*aux
            model(iz-1,ix,it) = model(iz-1,ix,it) - c5*data(iz,ix,it)*aux
            model(iz-2,ix,it) = model(iz-2,ix,it) - c4*data(iz,ix,it)*aux
            model(iz-3,ix,it) = model(iz-3,ix,it) - c3*data(iz,ix,it)*aux
            model(iz-4,ix,it) = model(iz-4,ix,it) - c2*data(iz,ix,it)*aux
            model(iz-5,ix,it) = model(iz-5,ix,it) - c1*data(iz,ix,it)*aux
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data,aux) PRIVATE(ix,iz,it)
      do it=1,nt
        do ix=6,nx-5
          do iz=6,nz-5
            data(iz,ix,it) = data(iz,ix,it) + (c1*(model(iz+4,ix,it) - model(iz-5,ix,it))+&
                                               c2*(model(iz+3,ix,it) - model(iz-4,ix,it))+&
                                               c3*(model(iz+2,ix,it) - model(iz-3,ix,it))+&
                                               c4*(model(iz+1,ix,it) - model(iz-2,ix,it))+&
                                               c5*(model(iz  ,ix,it) - model(iz-1,ix,it)))*aux
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

  subroutine vzdx(adj,add,model,data)
    real                                   :: aux
    logical,intent(in)                     :: adj,add
    real,dimension(nz,nx,nt),intent(inout) :: model,data

    aux=1.0/dx

    if(adj)then
      if(.not.add) model = 0.
!$OMP PARALLEL DO SHARED(model,data,aux) PRIVATE(ix,iz,it)
      do it=1,nt
        do ix=6,nx-5
          do iz=6,nz-5
            model(iz,ix+4,nt) = model(iz,ix+4,it) + c1*data(iz,ix,it)*aux
            model(iz,ix+3,it) = model(iz,ix+3,it) + c2*data(iz,ix,it)*aux
            model(iz,ix+2,it) = model(iz,ix+2,it) + c3*data(iz,ix,it)*aux
            model(iz,ix+1,it) = model(iz,ix+1,it) + c4*data(iz,ix,it)*aux
            model(iz,ix  ,it) = model(iz,ix  ,it) + c5*data(iz,ix,it)*aux
            model(iz,ix-1,it) = model(iz,ix-1,it) - c5*data(iz,ix,it)*aux
            model(iz,ix-2,it) = model(iz,ix-2,it) - c4*data(iz,ix,it)*aux
            model(iz,ix-3,it) = model(iz,ix-3,it) - c3*data(iz,ix,it)*aux
            model(iz,ix-4,it) = model(iz,ix-4,it) - c2*data(iz,ix,it)*aux
            model(iz,ix-5,it) = model(iz,ix-5,it) - c1*data(iz,ix,it)*aux
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
    else
      if(.not.add) data = 0.
!$OMP PARALLEL DO SHARED(model,data,aux) PRIVATE(ix,iz,it)
      do it=1,nt
        do ix=6,nx-5
          do iz=6,nz-5
            data(iz,ix,it) = data(iz,ix,it) + (c1*(model(iz,ix+4,it) - model(iz,ix-5,it))+&
                                               c2*(model(iz,ix+3,it) - model(iz,ix-4,it))+&
                                               c3*(model(iz,ix+2,it) - model(iz,ix-3,it))+&
                                               c4*(model(iz,ix+1,it) - model(iz,ix-2,it))+&
                                               c5*(model(iz,ix  ,it) - model(iz,ix-1,it)))*aux
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
    endif
  endsubroutine

endmodule
