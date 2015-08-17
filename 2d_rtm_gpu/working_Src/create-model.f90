! Main program to create 2D models
program create_model

  use sep

  implicit none

  logical                               :: lvl, water, anom, ref, scatter
  real                                  :: vs1,vs2
  integer                               :: nz,nx, iz, ix, thick, z1, nlvl, waterz, lvlvz1, lvlvz2
  integer                               :: anomz, anomx
  integer                               :: scatterz, scatterx
  real                                  :: dz,oz,dx,ox
  real                                  :: vp1, vp2, rho2, vz, rho1, lvlvs, lvlvp, lvlrho
  real                                  :: watervp, watervs, waterrho
  real                                  :: anomstr, anomsigma, aux
  real                                  :: scattervp, scattervs, scatterrho
  real, allocatable, dimension(:,:)     :: vp, vs, rho

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
  call from_param("watervp",watervp,1500.0)   !vp of water (1500.0 for default)
  call from_param("watervs",watervs,0.0)      !vp of water (0.0 for default)
  call from_param("waterrho",waterrho,1000.0) !vp of water (1000.0 for default)

  ! Low velocity layer
  call from_param("lvl",lvl,.false.)          !if there is a low velocity layer
  call from_param("nlvl",nlvl,0)              !thickness of lvl
  call from_param("lvlvz1",lvlvz1,0)          !vertical increase in vs in lvl
  call from_param("lvlvz2",lvlvz2,0)          !vertical increase in vp in lvl
  call from_param("lvlvp",lvlvp,-1.0)         !by default, it is a poisson solid
  call from_param("lvlvs",lvlvs,0.0)          !vs of lvl (vp is calculated)
  call from_param("lvlrho",lvlrho,1000.0)     !vs of lvl (vp is calculated)

  ! Background properties
  call from_param("vp1",vp1)                  !vp of background
  call from_param("vs1",vs1,-1.0)             !vs of background (-1 for vs=vp/2)
  call from_param("rho1",rho1)                !density of background
  call from_param("vz",vz,0.0)                !vertical increase in velocity

  ! Reflector
  call from_param("ref",ref,.false.)          !if there is a reflector
  call from_param("z1",z1,0)                  !depth of reflector
  call from_param("thick",thick,0)            !thickness of reflector    
  call from_param("vp2",vp2,-1.0)             !vp of reflector
  call from_param("vs2",vs2,-1.0)             !vs of reflector
  call from_param("rho2",rho2,-1.0)           !rho of reflector

  ! Gaussian anomaly
  call from_param("anom",anom,.false.)        !if there is a gaussian anomaly
  call from_param("anomz",anomz,0)            !depth of center of gaussian
  call from_param("anomx",anomx,0)            !x coord of center of gaussian
  call from_param("anomsigma",anomsigma,0.0)  !standard deviation of gaussian
  call from_param("anomstr",anomstr,0.0)      !strength of anomaly to background (>0 for posit, <0 for neg)

  ! Point scatterers
  call from_param("scatter",scatter,.false.)  !if there is a point scatterer
  call from_param("scatterz",scatterz,0)      !depth of scatterer
  call from_param("scatterx",scatterx,0)      !x coord of scatterer
  call from_param("scatterrho",scatterrho,0.0)!rho of scatterer
  call from_param("scattervp",scattervp,0.0)  !vp of scatterer
  call from_param("scattervs",scattervs,0.0)  !vs of scatterer

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
  call to_history("n1",nz,"vs")
  call to_history("o1",oz,"vs")
  call to_history("d1",dz,"vs")
  call to_history("n2",nx,"vs")
  call to_history("o2",ox,"vs")
  call to_history("d2",dx,"vs")
  call to_history("label1","z","vs")
  call to_history("label2","x","vs")

  ! Allocation of the input/output matrices
  allocate(vp(nz,nx))
  allocate(vs(nz,nx))
  allocate(rho(nz,nx))

  ! Initialize values
  vp=0.
  vs=0.
  rho=0.

  ! Create water layer
  if(water) then
    do ix=1,nx
      do iz=1,waterz
        vp(iz,ix)  = watervp
        vs(iz,ix)  = watervs
        rho(iz,ix) = waterrho
      enddo
    enddo
  endif

  ! Create low velocity layer
  ! Poisson solid vp=sqrt(2)*vs
  if(lvl) then
    do ix=1,nx
      do iz=1+waterz,waterz+nlvl
        vs(iz,ix)=lvlvs+(iz-waterz-1)*lvlvz1
        if(lvlvp.ne.-1.0)then
          vp(iz,ix)=lvlvp+(iz-water-1)*lvlvz2
        else
          vp(iz,ix)=sqrt(3.0)*vs(iz,ix)
        endif
        rho(iz,ix)=lvlrho
      enddo
    enddo
  endif

  ! Create background velocity
  do ix=1,nx
    do iz=1+waterz+nlvl,nz
      vp(iz,ix)=vp1+vz*(iz-waterz-nlvl)
      if (vs1.ne.-1.0) then
        vs(iz,ix)=vs1+0.5*vz*(iz-waterz-nlvl)
      else
        vs(iz,ix)=0.50*vp(iz,ix)
      endif
      rho(iz,ix)=rho1
    enddo
  enddo

  ! Insert fluid saturated horizontal layer
  if(ref) then
    do ix=1,nx
      do iz=z1,z1+thick
        vp(iz,ix)=vp2
        if(vs2.ne.-1.0)then
          vs(iz,ix)=vs2
        else
          vs(iz,ix)=0.50*vs(iz,ix)
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
        vp(iz,ix)=vp(iz,ix)*(1.0+anomstr*exp(-aux))
        vs(iz,ix)=vs(iz,ix)*(1.0+anomstr*exp(-aux))
        rho(iz,ix)=rho(iz,ix)*(1.0+anomstr*exp(-aux))
      enddo
    enddo
  endif

  ! Insert scatterer
  if(scatter) then
    vp(scatterz,scatterx)=scattervp
    vs(scatterz,scatterx)=scattervs
    rho(scatterz,scatterx)=scatterrho
  endif
 

  ! Write outputs
  call sep_write(vp)
  call sep_write(vs,"vs")
  call sep_write(rho,"rho")

  call sep_close

end program
