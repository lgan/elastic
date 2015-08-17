! Main program to create 2D models
program create_model_lame

  use sep

  implicit none

  logical                               :: lvl, water, anom, ref, circle
  real                                  :: mu1,mu2
  integer                               :: nz,nx, iz, ix, thick, z1, nlvl, waterz
  integer                               :: anomz, anomx, circlex, circlez, circleradius
  real                                  :: dz,oz,dx,ox
  real                                  :: lambda1, lambda2, rho2, vz, rho1, lvlmu, lvllambda, lvlrho, lvlvz1, lvlvz2
  real                                  :: waterlambda, watermu, waterrho
  real                                  :: anomstr, anomsigma, aux
  real                                  :: circlelambda,circlemu,circlerho
  real                                  :: anomlambda,anommu,anomrho
  real, allocatable, dimension(:,:)     :: lambda, mu, rho

  call sep_init()

  ! Model parameters
  ! Dimensions
  call from_param("nz",nz)                    !# of samples in depth
  call from_param("oz",oz,0.0)                !origin for depth coord
  call from_param("dz",dz)                    !grid spacing in depth
  call from_param("nx",nx)                    !# of samples in x
  call from_param("ox",ox,0.0)                !origin for x coord
  call from_param("dx",dx)                    !grid spacing in x

  ! Water layer
  call from_param("water",water,.false.)      !if there is a water layer
  call from_param("waterz",waterz,0)          !depth of water layer
  call from_param("waterlambda",waterlambda,2.25E9)   !lambda of water (2.25E9 for default)
  call from_param("watermu",watermu,0.0)      !mu of water (0.0 for default)
  call from_param("waterrho",waterrho,1000.0) !rho of water (1000.0 for default)

  ! Low velocity layer
  call from_param("lvl",lvl,.false.)          !if there is a low velocity layer
  call from_param("nlvl",nlvl,0)              !thickness of lvl
  call from_param("lvlvz1",lvlvz1,0.0)        !vertical increase in mu in lvl
  call from_param("lvlvz2",lvlvz2,0.0)        !vertical increase in lambda in lvl
  call from_param("lvllambda",lvllambda,-1.0) !by default, it is a poisson solid
  call from_param("lvlmu",lvlmu,0.0)          !mu of lvl (lambda is calculated)
  call from_param("lvlrho",lvlrho,1000.0)     !rho of lvl

  ! Background properties
  call from_param("lambda1",lambda1,0.0)      !lambda of background
  call from_param("mu1",mu1,-1.0)             !mu of background (-1 for mu=lambda/2)
  call from_param("rho1",rho1,1000.0)         !density of background
  call from_param("vz",vz,0.0)                !vertical increase in velocity

  ! Reflector
  call from_param("ref",ref,.false.)          !if there is a reflector
  call from_param("z1",z1,0)                  !depth of reflector
  call from_param("thick",thick,0)            !thickness of reflector    
  call from_param("lambda2",lambda2,-1.0)             !lambda of reflector
  call from_param("mu2",mu2,-1.0)             !mu of reflector
  call from_param("rho2",rho2,-1.0)           !rho of reflector

  ! Gaussian anomaly
  call from_param("anom",anom,.false.)        !if there is a gaussian anomaly
  call from_param("anomz",anomz,0)            !depth of center of gaussian
  call from_param("anomx",anomx,0)            !x coord of center of gaussian
  call from_param("anomsigma",anomsigma,0.0)  !standard deviation of gaussian
  call from_param("anomstr",anomstr,0.0)      !strength of anomaly to background (>0 for posit, <0 for neg)
  call from_param("anomrho",anomrho,0.0)      !density anomaly max. value
  call from_param("anomlambda",anomlambda,0.0)!lambda anomaly max. value
  call from_param("anommu",anommu,0.0)        !mu anomaly max. value

  ! Circle
  call from_param("circle",circle,.false.)    !if there is a disc-shaped reflector
  call from_param("circlez",circlez,0)          !depth of center of circle
  call from_param("circlex",circlex,0)          !x coord of center of circle
  call from_param("circleradius",circleradius,0)  !radius of circle
  call from_param("circlelambda",circlelambda,-1.0) !lambda of reflector
  call from_param("circlemu",circlemu,-1.0)         !mu of reflector
  call from_param("circlerho",circlerho,-1.0)       !rho of reflector

  ! Write SEP header parameter
  call to_history("n1",nz)
  call to_history("o1",oz)
  call to_history("d1",dz)
  call to_history("n2",nx)
  call to_history("o2",ox)
  call to_history("d2",dx)
  call to_history("label1","z")
  call to_history("label2","x")
  call to_history("n1",nz,"rho")
  call to_history("o1",oz,"rho")
  call to_history("d1",dz,"rho")
  call to_history("n2",nx,"rho")
  call to_history("o2",ox,"rho")
  call to_history("d2",dx,"rho")
  call to_history("label1","z","rho")
  call to_history("label2","x","rho")
  call to_history("n1",nz,"mu")
  call to_history("o1",oz,"mu")
  call to_history("d1",dz,"mu")
  call to_history("n2",nx,"mu")
  call to_history("o2",ox,"mu")
  call to_history("d2",dx,"mu")
  call to_history("label1","z","mu")
  call to_history("label2","x","mu")

  ! Allocation of the input/output matrices
  allocate(lambda(nz,nx))
  allocate(mu(nz,nx))
  allocate(rho(nz,nx))

  ! Initialize values
  lambda=0.
  mu=0.
  rho=0.

  ! Create water layer
  if(water) then
    do ix=1,nx
      do iz=1,waterz
        lambda(iz,ix)  = waterlambda
        mu(iz,ix)  = watermu
        rho(iz,ix) = waterrho
      enddo
    enddo
  endif

  ! Create low velocity layer
  ! Poisson solid 
  if(lvl) then
    do ix=1,nx
      do iz=1+waterz,waterz+nlvl
        mu(iz,ix)=lvlmu+(iz-waterz-1)*lvlvz1
        if(lvllambda.ne.-1.0)then
          lambda(iz,ix)=lvllambda+(iz-waterz-1)*lvlvz2
        else
          lambda(iz,ix)=mu(iz,ix)
        endif
        rho(iz,ix)=lvlrho
      enddo
    enddo
  endif

  ! Create background lame parameters
  do ix=1,nx
    do iz=1+waterz+nlvl,nz
      lambda(iz,ix)=lambda1+vz*(iz-waterz-nlvl)
      if (mu1.ne.-1.0) then
        mu(iz,ix)=mu1+0.5*vz*(iz-waterz-nlvl)
      else
        mu(iz,ix)=0.50*lambda(iz,ix)
      endif
      rho(iz,ix)=rho1
    enddo
  enddo

  ! Insert fluid saturated horizontal layer
  if(ref) then
    do ix=1,nx
      do iz=z1,z1+thick
        lambda(iz,ix)=lambda2
        if(mu2.ne.-1.0)then
          mu(iz,ix)=mu2
        else
          mu(iz,ix)=0.50*mu(iz,ix)
        endif
        rho(iz,ix)=rho2
      enddo
    enddo
  endif

  ! Insert gaussian anomaly
  if(anom) then
    do ix=1,nx
      do iz=1,nz
        aux = ((ix-anomx)*(ix-anomx)+(iz-anomz)*(iz-anomz))/(2*anomsigma*anomsigma)
        if (anomstr.eq.0.0) then
          lambda(iz,ix)= lambda(iz,ix)+(anomlambda-lambda(iz,ix))*exp(-aux)
          mu(iz,ix)=mu(iz,ix)+(anommu-mu(iz,ix))*exp(-aux)
          rho(iz,ix)=rho(iz,ix)+(anomrho-rho(iz,ix))*exp(-aux)
        else
          lambda(iz,ix)=lambda(iz,ix)*(1.0+anomstr*exp(-aux))
          mu(iz,ix)=mu(iz,ix)*(1.0+anomstr*exp(-aux))
          rho(iz,ix)=rho(iz,ix)*(1.0+anomstr*exp(-aux))
        endif
      enddo
    enddo
  endif

  ! Insert circle
  if(circle) then
    do ix=1,nx
      do iz=1,nz
        aux = ((ix-circlex)*(ix-circlex)+(iz-circlez)*(iz-circlez))
        if (aux.le.circleradius) then
          lambda(iz,ix)=circlelambda
          mu(iz,ix)=circlemu
          rho(iz,ix)=circlerho
        endif
      enddo
    enddo
  endif

  ! Write outputs
  call sep_write(lambda)
  call sep_write(mu,"mu")
  call sep_write(rho,"rho")

  call sep_close

end program
