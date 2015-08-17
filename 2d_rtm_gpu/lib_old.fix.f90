# 1 "<stdin>"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "<stdin>"
! finite difference module
module elasticlib_mod

  use sep

  implicit none
  integer, private :: nx, nz, ik
  real, private :: c1,c2,c3,c4,c5,r1,r2,r3

  contains

  subroutine elasticmod_init(nx_in,nz_in)
    integer :: nx_in, nz_in

    nx = nx_in
    nz = nz_in
  endsubroutine

  function dx4(model,iz,ix) result(data)
    integer,intent(in) :: iz,ix
    real,dimension(nz,nx) :: model, data

    c1=-1.0/24.0
    c2= 9.0/8.0

    data=c1*(model(iz,ix+1)-model(iz,ix-2))+c2*(model(iz,ix)-model(iz,ix-1))

  endfunction

  function dz4(model,iz,ix) result(data)
    integer,intent(in) :: iz,ix
    real,dimension(nz,nx) :: model, data

    c1=-1.0/24.0
    c2= 9.0/8.0

    data=c1*(model(iz+1,ix)-model(iz-2,ix))+c2*(model(iz,ix)-model(iz-1,ix))

  endfunction

    if (sp.eq.4) then
    elseif (sp.eq.10) then
    !Coefficients from Liu et al. 2009, Geo.J.Int.
      c1= 19845.0/16384.0
      c2= -735.0/8192.0
      c3= 567.0/40960.0
      c4= -405.0/229376.0
      c5= 35.0/294912.0
    endif

    ! Coefficients from Mathematica
    r1= 21.0/23.0
    r2= 3.0/23.0
    r3= -1.0/23.0
    ! Coefficients from Ghrist, 2000
    !r1= 99.0/96.0
    !r2= -3.0/96.0
    !r3= 1.0/96.0
    !4th order time derivative by hand, backward coefficients
    !r1= 18.0/11.0
    !r2=-18.0/22.0
    !r3= 6.0/33.0

  end subroutine

  function elasticmod_op(it,source,Vx,Vz,sigmaxx,sigmazz,sigmaxz,flag) result(stat)

    integer,intent(in):: it
    real,intent(in) :: source
    real,dimension(:) :: Vx, Vz
    real,dimension(:) :: sigmaxx, sigmazz, sigmaxz
    integer,dimension(:) :: flag
    integer :: stat

    if (sp.eq.4.and.ti.eq.2) then
      call fd10s2t_op2(it,source,Vx,Vz,sigmaxx,sigmazz,sigmaxz)
    elseif (sp.eq.10.and.ti.eq.4) then
    call fd10s4t_op2(it,source,Vx,Vz,sigmaxx,sigmazz,sigmaxz,flag)
    endif
    stat=0
  end function
  subroutine fd10s2t_op2(it,source,Vx,Vz,sigmaxx,sigmazz,sigmaxz)
    integer :: it
    real :: source
    real,dimension(nz,nx) :: Vx,Vz
    real :: Exx,Ezz,Exz, rhox, rhoz
    real,dimension(nz,nx) :: sigmaxx,sigmazz,sigmaxz
    real :: aux2,dVxx,dVxz,dVzx,dVzz

    ! Boundary conditions for stress
    ! Free surface (Stress = 0) at surface
    !do ik=1,3
    ! do ix=1,nx
    ! do iz=1,6
    ! sigmazz(iz,ix,ik) = 0.0
    ! sigmaxz(iz,ix,ik) = 0.0
    ! enddo
    ! enddo
    !enddo


    ! Equations of momentum conservation
    do ix=6,nx-5
      do iz=6,nz-5

         rhox = 0.5*(rho(iz,ix) + rho(iz,ix-1))
         rhoz = 0.5*(rho(iz,ix) + rho(iz-1,ix))

         Vx(iz,ix) = Vx(iz,ix) + aux1*rhox &
                   *(c1*(sigmaxx(iz ,ix )-sigmaxx(iz ,ix-1) &
                        +sigmaxz(iz ,ix )-sigmaxz(iz-1,ix ))&
                   + c2*(sigmaxx(iz ,ix+1)-sigmaxx(iz ,ix-2) &
                        +sigmaxz(iz+1,ix )-sigmaxz(iz-2,ix ))&
                   + c3*(sigmaxx(iz ,ix+2)-sigmaxx(iz ,ix-3) &
                        +sigmaxz(iz+2,ix )-sigmaxz(iz-3,ix ))&
                   + c4*(sigmaxx(iz ,ix+3)-sigmaxx(iz ,ix-4) &
                        +sigmaxz(iz+3,ix )-sigmaxz(iz-4,ix ))&
                   + c5*(sigmaxx(iz ,ix+4)-sigmaxx(iz ,ix-5) &
                        +sigmaxz(iz+4,ix )-sigmaxz(iz-5,ix )))

         Vz(iz,ix) = Vz(iz,ix) + aux1*rhoz &
                   *(c1*(sigmazz(iz ,ix )-sigmazz(iz-1,ix ) &
                        +sigmaxz(iz ,ix )-sigmaxz(iz ,ix-1))&
                   + c2*(sigmazz(iz+1,ix )-sigmazz(iz-2,ix ) &
                        +sigmaxz(iz ,ix+1)-sigmaxz(iz ,ix-2))&
                   + c3*(sigmazz(iz+2,ix )-sigmazz(iz-3,ix ) &
                        +sigmaxz(iz ,ix+2)-sigmaxz(iz ,ix-3))&
                   + c4*(sigmazz(iz+3,ix )-sigmazz(iz-4,ix ) &
                        +sigmaxz(iz ,ix+3)-sigmaxz(iz ,ix-4))&
                   + c5*(sigmazz(iz+4,ix )-sigmazz(iz-5,ix ) &
                        +sigmaxz(iz ,ix+4)-sigmaxz(iz ,ix-5)))
      enddo
    enddo

    ! Boundary conditions for displacements
    ! Displacements at top
    !do ik=1,3
    ! do ix=1,nx
    ! do iz=1,5
    ! Vx(iz,ix,ik)=0.0
    ! Vz(iz,ix,ik)=0.0
    ! enddo
    ! enddo
    !enddo

    ! Stress strain relations in the X and Z directions
    do ix=6,nx-5
      do iz=6,nz-5

        !Derivatives for the calculating the strains
        dVxx = c1*(Vx(iz,ix+1)-Vx(iz,ix )) &
             + c2*(Vx(iz,ix+2)-Vx(iz,ix-1)) &
             + c3*(Vx(iz,ix+3)-Vx(iz,ix-2)) &
             + c4*(Vx(iz,ix+4)-Vx(iz,ix-3)) &
             + c5*(Vx(iz,ix+5)-Vx(iz,ix-4))
        dVxz = c1*(Vx(iz+1,ix)-Vx(iz ,ix)) &
             + c2*(Vx(iz+2,ix)-Vx(iz-1,ix)) &
             + c3*(Vx(iz+3,ix)-Vx(iz-2,ix)) &
             + c4*(Vx(iz+4,ix)-Vx(iz-3,ix)) &
             + c5*(Vx(iz+5,ix)-Vx(iz-4,ix))
        dVzx = c1*(Vz(iz,ix )-Vz(iz,ix-1)) &
             + c2*(Vz(iz,ix+1)-Vz(iz,ix-2)) &
             + c3*(Vz(iz,ix+2)-Vz(iz,ix-3)) &
             + c4*(Vz(iz,ix+3)-Vz(iz,ix-4)) &
             + c5*(Vz(iz,ix+4)-Vz(iz,ix-5))
        dVzz = c1*(Vz(iz ,ix)-Vz(iz-1,ix)) &
             + c2*(Vz(iz+1,ix)-Vz(iz-2,ix)) &
             + c3*(Vz(iz+2,ix)-Vz(iz-3,ix)) &
             + c4*(Vz(iz+3,ix)-Vz(iz-4,ix)) &
             + c5*(Vz(iz+4,ix)-Vz(iz-5,ix))

        !Finite strain approximation (includes second order terms)
        Exx = dVxx !+ 0.5*(dVxx*dVxx+dVxz*dVxz)
        Ezz = dVzz !+ 0.5*(dVzz*dVzz+dVxz*dVxz)
        Exz = 0.5*(dVxz+dVzx)!+0.5*(dVxx*dVxz+dVzx*dVzz)

        !Stress
        sigmaxx(iz,ix) = sigmaxx(iz,ix) &
                       + dt*(aux2*Exx+lambda(iz,ix)*Ezz)
        sigmazz(iz,ix) = sigmazz(iz,ix) &
                       + dt*(aux2*Ezz+lambda(iz,ix)*Exx)
        sigmaxz(iz,ix) = sigmaxx(iz,ix) &
                       + dt*(2*mi(iz,ix)*Exz)
      enddo
    enddo

    ! Placement of source
    source = source*(dt/(dx*dx))*0.25
    sigmaxx(zsource ,xsource-1) = sigmaxx(zsource ,xsource-1) + source
    sigmaxx(zsource ,xsource ) = sigmaxx(zsource ,xsource ) + source
    sigmaxx(zsource+1,xsource-1) = sigmaxx(zsource+1,xsource-1) + source
    sigmaxx(zsource+1,xsource ) = sigmaxx(zsource+1,xsource ) + source
    sigmazz(zsource ,xsource-1) = sigmazz(zsource ,xsource-1) + source
    sigmazz(zsource ,xsource ) = sigmazz(zsource ,xsource ) + source
    sigmazz(zsource+1,xsource-1) = sigmazz(zsource+1,xsource-1) + source
    sigmazz(zsource+1,xsource ) = sigmazz(zsource+1,xsource ) + source

  endsubroutine

  subroutine fd10s4t_op2(it,source,Vx,Vz,sigmaxx,sigmazz,sigmaxz,flag)
    integer :: it
    real :: source
    real,dimension(nz,nx,3) :: Vx,Vz
    real :: Exx,Ezz,Exz, rhox, rhoz
    integer,dimension(3) :: flag
    real,dimension(nz,nx,3) :: sigmaxx,sigmazz,sigmaxz
    real :: aux2,dVxx,dVxz,dVzx,dVzz

    ! Boundary conditions for stress
    ! Free surface (Stress = 0) at surface
    !do ik=1,3
    ! do ix=1,nx
    ! do iz=1,6
    ! sigmazz(iz,ix,ik) = 0.0
    ! sigmaxz(iz,ix,ik) = 0.0
    ! enddo
    ! enddo
    !enddo

    ! Equations of momentum conservation
    do ix=6,nx-5
      do iz=6,nz-5

         rhox = 0.5*(rho(iz,ix) + rho(iz,ix-1))
         rhoz = 0.5*(rho(iz,ix) + rho(iz-1,ix))

         Vx(iz,ix,flag(3)) = r1*Vx(iz,ix,flag(2)) &
                           + r2*Vx(iz,ix,flag(1)) &
                           + r3*Vx(iz,ix,flag(3)) &
                           + aux1*rhox*(48.0/23.0) &
                   *(c1*(sigmaxx(iz ,ix ,flag(2))-sigmaxx(iz ,ix-1,flag(2)) &
                        +sigmaxz(iz ,ix ,flag(2))-sigmaxz(iz-1,ix ,flag(2)))&
                   + c2*(sigmaxx(iz ,ix+1,flag(2))-sigmaxx(iz ,ix-2,flag(2)) &
                        +sigmaxz(iz+1,ix ,flag(2))-sigmaxz(iz-2,ix ,flag(2)))&
                   + c3*(sigmaxx(iz ,ix+2,flag(2))-sigmaxx(iz ,ix-3,flag(2)) &
                        +sigmaxz(iz+2,ix ,flag(2))-sigmaxz(iz-3,ix ,flag(2)))&
                   + c4*(sigmaxx(iz ,ix+3,flag(2))-sigmaxx(iz ,ix-4,flag(2)) &
                        +sigmaxz(iz+3,ix ,flag(2))-sigmaxz(iz-4,ix ,flag(2)))&
                   + c5*(sigmaxx(iz ,ix+4,flag(2))-sigmaxx(iz ,ix-5,flag(2)) &
                        +sigmaxz(iz+4,ix ,flag(2))-sigmaxz(iz-5,ix ,flag(2))))

         Vz(iz,ix,flag(3)) = r1*Vz(iz,ix,flag(2)) &
                           + r2*Vz(iz,ix,flag(1)) &
                           + r3*Vz(iz,ix,flag(3)) &
                           + aux1*rhoz*(48.0/23.0) &
                   *(c1*(sigmazz(iz ,ix ,flag(2))-sigmazz(iz-1,ix ,flag(2)) &
                        +sigmaxz(iz ,ix ,flag(2))-sigmaxz(iz ,ix-1,flag(2)))&
                   + c2*(sigmazz(iz+1,ix ,flag(2))-sigmazz(iz-2,ix ,flag(2)) &
                        +sigmaxz(iz ,ix+1,flag(2))-sigmaxz(iz ,ix-2,flag(2)))&
                   + c3*(sigmazz(iz+2,ix ,flag(2))-sigmazz(iz-3,ix ,flag(2)) &
                        +sigmaxz(iz ,ix+2,flag(2))-sigmaxz(iz ,ix-3,flag(2)))&
                   + c4*(sigmazz(iz+3,ix ,flag(2))-sigmazz(iz-4,ix ,flag(2)) &
                        +sigmaxz(iz ,ix+3,flag(2))-sigmaxz(iz ,ix-4,flag(2)))&
                   + c5*(sigmazz(iz+4,ix ,flag(2))-sigmazz(iz-5,ix ,flag(2)) &
                        +sigmaxz(iz ,ix+4,flag(2))-sigmaxz(iz ,ix-5,flag(2))))
      enddo
    enddo

    ! Boundary conditions for displacements
    ! Displacements at top
    !do ik=1,3
    ! do ix=1,nx
    ! do iz=1,5
    ! Vx(iz,ix,ik)=0.0
    ! Vz(iz,ix,ik)=0.0
    ! enddo
    ! enddo
    !enddo

    ! Stress strain relations in the X and Z directions
    do ix=6,nx-5
      do iz=6,nz-5
        aux2 = lambda(iz,ix)+2.0*mi(iz,ix)

        !Derivatives for calculating the strains
        dVxx = c1*(Vx(iz,ix+1,flag(3))-Vx(iz,ix ,flag(3))) &
             + c2*(Vx(iz,ix+2,flag(3))-Vx(iz,ix-1,flag(3))) &
             + c3*(Vx(iz,ix+3,flag(3))-Vx(iz,ix-2,flag(3))) &
             + c4*(Vx(iz,ix+4,flag(3))-Vx(iz,ix-3,flag(3))) &
             + c5*(Vx(iz,ix+5,flag(3))-Vx(iz,ix-4,flag(3)))
        dVxz = c1*(Vx(iz+1,ix,flag(3))-Vx(iz ,ix,flag(3))) &
             + c2*(Vx(iz+2,ix,flag(3))-Vx(iz-1,ix,flag(3))) &
             + c3*(Vx(iz+3,ix,flag(3))-Vx(iz-2,ix,flag(3))) &
             + c4*(Vx(iz+4,ix,flag(3))-Vx(iz-3,ix,flag(3))) &
             + c5*(Vx(iz+5,ix,flag(3))-Vx(iz-4,ix,flag(3)))
        dVzx = c1*(Vz(iz,ix ,flag(3))-Vz(iz,ix-1,flag(3))) &
             + c2*(Vz(iz,ix+1,flag(3))-Vz(iz,ix-2,flag(3))) &
             + c3*(Vz(iz,ix+2,flag(3))-Vz(iz,ix-3,flag(3))) &
             + c4*(Vz(iz,ix+3,flag(3))-Vz(iz,ix-4,flag(3))) &
             + c5*(Vz(iz,ix+4,flag(3))-Vz(iz,ix-5,flag(3)))
        dVzz = c1*(Vz(iz ,ix,flag(3))-Vz(iz-1,ix,flag(3))) &
             + c2*(Vz(iz+1,ix,flag(3))-Vz(iz-2,ix,flag(3))) &
             + c3*(Vz(iz+2,ix,flag(3))-Vz(iz-3,ix,flag(3))) &
             + c4*(Vz(iz+3,ix,flag(3))-Vz(iz-4,ix,flag(3))) &
             + c5*(Vz(iz+4,ix,flag(3))-Vz(iz-5,ix,flag(3)))

        !Finite strain approximation (includes second order terms)
        Exx = dVxx !+ 0.5*(dVxx*dVxx+dVxz*dVxz)
        Ezz = dVzz !+ 0.5*(dVzz*dVzz+dVxz*dVxz)
        Exz = 0.5*(dVxz+dVzx)!+0.5*(dVxx*dVxz+dVzx*dVzz)

        !Stress
        sigmaxx(iz,ix,flag(3)) = r1*sigmaxx(iz,ix,flag(2)) &
                               + r2*sigmaxx(iz,ix,flag(1)) &
                               + r3*sigmaxx(iz,ix,flag(3)) &
                               + aux1*(48.0/23.0)*(aux2*Exx+lambda(iz,ix)*Ezz)
        sigmazz(iz,ix,flag(3)) = r1*sigmazz(iz,ix,flag(2)) &
                               + r2*sigmazz(iz,ix,flag(1)) &
                               + r3*sigmazz(iz,ix,flag(3)) &
                               + aux1*(48.0/23.0)*(aux2*Ezz+lambda(iz,ix)*Exx)
        sigmaxz(iz,ix,flag(3)) = r1*sigmaxx(iz,ix,flag(2)) &
                               + r2*sigmaxx(iz,ix,flag(1)) &
                               + r3*sigmaxx(iz,ix,flag(3)) &
                               + aux1*(48.0/23.0)*(2*mi(iz,ix)*Exz)
        dVxx = 0.0
        dVxz = 0.0
        dVzx = 0.0
        dVzz = 0.0
        Exx = 0.0
        Ezz = 0.0
        Exz = 0.0
      enddo
    enddo

    ! Placement of source (explosive source)
    source = 1000000*source*(dt/(dx*dx))*0.25
    sigmaxx(zsource ,xsource-1,flag(3)) = sigmaxx(zsource ,xsource-1,flag(3)) + source
    sigmaxx(zsource ,xsource ,flag(3)) = sigmaxx(zsource ,xsource ,flag(3)) + source
    sigmaxx(zsource+1,xsource-1,flag(3)) = sigmaxx(zsource+1,xsource-1,flag(3)) + source
    sigmaxx(zsource+1,xsource ,flag(3)) = sigmaxx(zsource+1,xsource ,flag(3)) + source
    sigmazz(zsource ,xsource-1,flag(3)) = sigmazz(zsource ,xsource-1,flag(3)) + source
    sigmazz(zsource ,xsource ,flag(3)) = sigmazz(zsource ,xsource ,flag(3)) + source
    sigmazz(zsource+1,xsource-1,flag(3)) = sigmazz(zsource+1,xsource-1,flag(3)) + source
    sigmazz(zsource+1,xsource ,flag(3)) = sigmazz(zsource+1,xsource ,flag(3)) + source

  endsubroutine

endmodule
