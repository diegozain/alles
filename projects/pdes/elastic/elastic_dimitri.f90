program seismic_CPML_2D_iso_second

! 2D elastic finite-difference code in velocity and stress formulation
! with Convolutional-PML (C-PML) absorbing conditions for an isotropic medium

! Dimitri Komatitsch, University of Pau, France, April 2007.

! The second-order staggered-grid formulation of Madariaga (1976) and Virieux (1986) is used:
!
!            ^ y
!            |
!            |
!
!            +-------------------+
!            |                   |
!            |                   |
!            |                   |
!            |                   |
!            |        v_y        |
!   sigma_xy +---------+         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            |         |         |
!            +---------+---------+  ---> x
!           v_x    sigma_xx
!                  sigma_yy
!


implicit none

! total number of grid points in each direction of the grid
integer, parameter :: NX = 101
integer, parameter :: NY = 641

! size of a grid cell
double precision, parameter :: DELTAX = 10.d0
double precision, parameter :: DELTAY = DELTAX

! flags to add PML layers to the edges of the grid
logical, parameter :: USE_PML_XMIN = .true.
logical, parameter :: USE_PML_XMAX = .true.
logical, parameter :: USE_PML_YMIN = .true.
logical, parameter :: USE_PML_YMAX = .true.

! thickness of the PML layer in grid points
integer, parameter :: NPOINTS_PML = 10

! P-velocity, S-velocity and density
double precision, parameter :: cp = 3300.d0
double precision, parameter :: cs = cp / 1.732d0
double precision, parameter :: density = 2800.d0

! total number of time steps
integer, parameter :: NSTEP = 2000

! time step in seconds
double precision, parameter :: DELTAT = 2.d-3

! parameters for the source
double precision, parameter :: f0 = 7.d0
double precision, parameter :: t0 = 1.20d0 / f0
double precision, parameter :: factor = 1.d7

! source
integer, parameter :: ISOURCE = NX - 2*NPOINTS_PML - 1
integer, parameter :: JSOURCE = 2 * NY / 3 + 1
double precision, parameter :: xsource = (ISOURCE - 1) * DELTAX
double precision, parameter :: ysource = (JSOURCE - 1) * DELTAY
! angle of source force in degrees and clockwise, with respect to the vertical (Y) axis
double precision, parameter :: ANGLE_FORCE = 135.d0

! value of PI
double precision, parameter :: PI = 3.141592653589793238462643d0

! conversion from degrees to radians
double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0

! zero
double precision, parameter :: ZERO = 0.d0

! large value for maximum
double precision, parameter :: HUGEVAL = 1.d+30

! velocity threshold above which we consider that the code became unstable
double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
double precision, dimension(NX,NY) :: vx,vy,sigma_xx,sigma_yy,sigma_xy,lambda,mu,rho

! to interpolate material parameters at the right location in the staggered grid cell
double precision lambda_half_x,mu_half_x,lambda_plus_two_mu_half_x,mu_half_y,rho_half_x_half_y

! for evolution of total energy in the medium
double precision :: epsilon_xx,epsilon_yy,epsilon_xy
double precision, dimension(NSTEP) :: total_energy_kinetic,total_energy_potential

! power to compute d0 profile
double precision, parameter :: NPOWER = 2.d0

! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
double precision, parameter :: K_MAX_PML = 1.d0
double precision, parameter :: ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte

! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here
double precision, dimension(NX,NY) :: &
    memory_dvx_dx, &
    memory_dvx_dy, &
    memory_dvy_dx, &
    memory_dvy_dy, &
    memory_dsigma_xx_dx, &
    memory_dsigma_yy_dy, &
    memory_dsigma_xy_dx, &
    memory_dsigma_xy_dy

double precision :: &
    value_dvx_dx, &
    value_dvx_dy, &
    value_dvy_dx, &
    value_dvy_dy, &
    value_dsigma_xx_dx, &
    value_dsigma_yy_dy, &
    value_dsigma_xy_dx, &
    value_dsigma_xy_dy

! 1D arrays for the damping profiles
double precision, dimension(NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half
double precision, dimension(NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half

double precision :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
double precision :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized

! for the source
double precision :: a,t,force_x,force_y,source_term

integer :: i,j,it,irec

double precision :: Courant_number,velocnorm

!---
!--- program starts here
!---

print *
print *,'2D elastic finite-difference code in velocity and stress formulation with C-PML'
print *

! display size of the model
print *
print *,'NX = ',NX
print *,'NY = ',NY
print *
print *,'size of the model along X = ',(NX - 1) * DELTAX
print *,'size of the model along Y = ',(NY - 1) * DELTAY
print *
print *,'Total number of grid points = ',NX * NY
print *

!--- define profile of absorption in PML region

! thickness of the PML layer in meters
thickness_PML_x = NPOINTS_PML * DELTAX
thickness_PML_y = NPOINTS_PML * DELTAY

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
Rcoef = 0.001d0

! check that NPOWER is okay
if (NPOWER < 1) stop 'NPOWER must be greater than 1'

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
d0_x = - (NPOWER + 1) * cp * log(Rcoef) / (2.d0 * thickness_PML_x)
d0_y = - (NPOWER + 1) * cp * log(Rcoef) / (2.d0 * thickness_PML_y)

print *,'d0_x = ',d0_x
print *,'d0_y = ',d0_y
print *

d_x(:) = ZERO
d_x_half(:) = ZERO
K_x(:) = 1.d0
K_x_half(:) = 1.d0
alpha_x(:) = ZERO
alpha_x_half(:) = ZERO
a_x(:) = ZERO
a_x_half(:) = ZERO

d_y(:) = ZERO
d_y_half(:) = ZERO
K_y(:) = 1.d0
K_y_half(:) = 1.d0
alpha_y(:) = ZERO
alpha_y_half(:) = ZERO
a_y(:) = ZERO
a_y_half(:) = ZERO

! damping in the X direction

! origin of the PML layer (position of right edge minus thickness, in meters)
xoriginleft = thickness_PML_x
xoriginright = (NX-1)*DELTAX - thickness_PML_x

do i = 1,NX

! abscissa of current grid point along the damping profile
  xval = DELTAX * dble(i-1)

!---------- left edge
  if (USE_PML_XMIN) then

! define damping profile at the grid points
    abscissa_in_PML = xoriginleft - xval
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_x
      d_x(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

! define damping profile at half the grid points
    abscissa_in_PML = xoriginleft - (xval + DELTAX/2.d0)
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_x
      d_x_half(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

  endif

!---------- right edge
  if (USE_PML_XMAX) then

! define damping profile at the grid points
    abscissa_in_PML = xval - xoriginright
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_x
      d_x(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

! define damping profile at half the grid points
    abscissa_in_PML = xval + DELTAX/2.d0 - xoriginright
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_x
      d_x_half(i) = d0_x * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

  endif

! just in case, for -5 at the end
  if (alpha_x(i) < ZERO) alpha_x(i) = ZERO
  if (alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

  b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT)
  b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DELTAT)

! this to avoid division by zero outside the PML
  if (abs(d_x(i)) > 1.d-6) a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
  if (abs(d_x_half(i)) > 1.d-6) a_x_half(i) = d_x_half(i) * &
    (b_x_half(i) - 1.d0) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)))

enddo

! damping in the Y direction

! origin of the PML layer (position of right edge minus thickness, in meters)
yoriginbottom = thickness_PML_y
yorigintop = (NY-1)*DELTAY - thickness_PML_y

do j = 1,NY

! abscissa of current grid point along the damping profile
  yval = DELTAY * dble(j-1)

!---------- bottom edge
  if (USE_PML_YMIN) then

! define damping profile at the grid points
    abscissa_in_PML = yoriginbottom - yval
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

! define damping profile at half the grid points
    abscissa_in_PML = yoriginbottom - (yval + DELTAY/2.d0)
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y_half(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

  endif

!---------- top edge
  if (USE_PML_YMAX) then

! define damping profile at the grid points
    abscissa_in_PML = yval - yorigintop
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

! define damping profile at half the grid points
    abscissa_in_PML = yval + DELTAY/2.d0 - yorigintop
    if (abscissa_in_PML >= ZERO) then
      abscissa_normalized = abscissa_in_PML / thickness_PML_y
      d_y_half(j) = d0_y * abscissa_normalized**NPOWER
! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-2
      K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
      alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
    endif

  endif

  b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DELTAT)
  b_y_half(j) = exp(- (d_y_half(j) / K_y_half(j) + alpha_y_half(j)) * DELTAT)

! this to avoid division by zero outside the PML
  if (abs(d_y(j)) > 1.d-6) a_y(j) = d_y(j) * (b_y(j) - 1.d0) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)))
  if (abs(d_y_half(j)) > 1.d-6) a_y_half(j) = d_y_half(j) * &
    (b_y_half(j) - 1.d0) / (K_y_half(j) * (d_y_half(j) + K_y_half(j) * alpha_y_half(j)))

enddo

! compute the Lame parameters and density
do j = 1,NY
  do i = 1,NX
      rho(i,j) = density
      mu(i,j) = density*cs*cs
      lambda(i,j) = density*(cp*cp - 2.d0*cs*cs)
  enddo
enddo

! print position of the source
print *,'Position of the source:'
print *
print *,'x = ',xsource
print *,'y = ',ysource
print *

! check the Courant stability condition for the explicit time scheme
! R. Courant et K. O. Friedrichs et H. Lewy (1928)
Courant_number = cp * DELTAT * sqrt(1.d0/DELTAX**2 + 1.d0/DELTAY**2)
print *,'Courant number is ',Courant_number
print *
if (Courant_number > 1.d0) stop 'time step is too large, simulation will be unstable'

! initialize arrays
vx(:,:) = ZERO
vy(:,:) = ZERO
sigma_xx(:,:) = ZERO
sigma_yy(:,:) = ZERO
sigma_xy(:,:) = ZERO

! PML
memory_dvx_dx(:,:) = ZERO
memory_dvx_dy(:,:) = ZERO
memory_dvy_dx(:,:) = ZERO
memory_dvy_dy(:,:) = ZERO
memory_dsigma_xx_dx(:,:) = ZERO
memory_dsigma_yy_dy(:,:) = ZERO
memory_dsigma_xy_dx(:,:) = ZERO
memory_dsigma_xy_dy(:,:) = ZERO

! initialize total energy
total_energy_kinetic(:) = ZERO
total_energy_potential(:) = ZERO

!---
!---  beginning of time loop
!---

do it = 1,NSTEP

!------------------------------------------------------------
! compute stress sigma and update memory variables for C-PML
!------------------------------------------------------------

do j = 2,NY
  do i = 1,NX-1

! interpolate material parameters at the right location in the staggered grid cell
    lambda_half_x = 0.5d0 * (lambda(i+1,j) + lambda(i,j))
    mu_half_x = 0.5d0 * (mu(i+1,j) + mu(i,j))
    lambda_plus_two_mu_half_x = lambda_half_x + 2.d0 * mu_half_x

    value_dvx_dx = (vx(i+1,j) - vx(i,j)) / DELTAX
    value_dvy_dy = (vy(i,j) - vy(i,j-1)) / DELTAY

    memory_dvx_dx(i,j) = b_x_half(i) * memory_dvx_dx(i,j) + a_x_half(i) * value_dvx_dx
    memory_dvy_dy(i,j) = b_y(j) * memory_dvy_dy(i,j) + a_y(j) * value_dvy_dy

    value_dvx_dx = value_dvx_dx / K_x_half(i) + memory_dvx_dx(i,j)
    value_dvy_dy = value_dvy_dy / K_y(j) + memory_dvy_dy(i,j)

    sigma_xx(i,j) = sigma_xx(i,j) + &
       (lambda_plus_two_mu_half_x * value_dvx_dx + lambda_half_x * value_dvy_dy) * DELTAT

    sigma_yy(i,j) = sigma_yy(i,j) + &
       (lambda_half_x * value_dvx_dx + lambda_plus_two_mu_half_x * value_dvy_dy) * DELTAT

  enddo
enddo

do j = 1,NY-1
  do i = 2,NX

! interpolate material parameters at the right location in the staggered grid cell
    mu_half_y = 0.5d0 * (mu(i,j+1) + mu(i,j))

    value_dvy_dx = (vy(i,j) - vy(i-1,j)) / DELTAX
    value_dvx_dy = (vx(i,j+1) - vx(i,j)) / DELTAY

    memory_dvy_dx(i,j) = b_x(i) * memory_dvy_dx(i,j) + a_x(i) * value_dvy_dx
    memory_dvx_dy(i,j) = b_y_half(j) * memory_dvx_dy(i,j) + a_y_half(j) * value_dvx_dy

    value_dvy_dx = value_dvy_dx / K_x(i) + memory_dvy_dx(i,j)
    value_dvx_dy = value_dvx_dy / K_y_half(j) + memory_dvx_dy(i,j)

    sigma_xy(i,j) = sigma_xy(i,j) + mu_half_y * (value_dvy_dx + value_dvx_dy) * DELTAT

  enddo
enddo

!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------

do j = 2,NY
  do i = 2,NX

    value_dsigma_xx_dx = (sigma_xx(i,j) - sigma_xx(i-1,j)) / DELTAX
    value_dsigma_xy_dy = (sigma_xy(i,j) - sigma_xy(i,j-1)) / DELTAY

    memory_dsigma_xx_dx(i,j) = b_x(i) * memory_dsigma_xx_dx(i,j) + a_x(i) * value_dsigma_xx_dx
    memory_dsigma_xy_dy(i,j) = b_y(j) * memory_dsigma_xy_dy(i,j) + a_y(j) * value_dsigma_xy_dy

    value_dsigma_xx_dx = value_dsigma_xx_dx / K_x(i) + memory_dsigma_xx_dx(i,j)
    value_dsigma_xy_dy = value_dsigma_xy_dy / K_y(j) + memory_dsigma_xy_dy(i,j)

    vx(i,j) = vx(i,j) + (value_dsigma_xx_dx + value_dsigma_xy_dy) * DELTAT / rho(i,j)

  enddo
enddo

do j = 1,NY-1
  do i = 1,NX-1

! interpolate density at the right location in the staggered grid cell
    rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))

    value_dsigma_xy_dx = (sigma_xy(i+1,j) - sigma_xy(i,j)) / DELTAX
    value_dsigma_yy_dy = (sigma_yy(i,j+1) - sigma_yy(i,j)) / DELTAY

    memory_dsigma_xy_dx(i,j) = b_x_half(i) * memory_dsigma_xy_dx(i,j) + a_x_half(i) * value_dsigma_xy_dx
    memory_dsigma_yy_dy(i,j) = b_y_half(j) * memory_dsigma_yy_dy(i,j) + a_y_half(j) * value_dsigma_yy_dy

    value_dsigma_xy_dx = value_dsigma_xy_dx / K_x_half(i) + memory_dsigma_xy_dx(i,j)
    value_dsigma_yy_dy = value_dsigma_yy_dy / K_y_half(j) + memory_dsigma_yy_dy(i,j)

    vy(i,j) = vy(i,j) + (value_dsigma_xy_dx + value_dsigma_yy_dy) * DELTAT / rho_half_x_half_y

  enddo
enddo

! add the source (force vector located at a given grid point)
a = pi*pi*f0*f0
t = dble(it-1)*DELTAT

! Gaussian
! source_term = factor * exp(-a*(t-t0)**2)

! first derivative of a Gaussian
source_term = - factor * 2.d0*a*(t-t0)*exp(-a*(t-t0)**2)

! Ricker source time function (second derivative of a Gaussian)
! source_term = factor * (1.d0 - 2.d0*a*(t-t0)**2)*exp(-a*(t-t0)**2)

force_x = sin(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term
force_y = cos(ANGLE_FORCE * DEGREES_TO_RADIANS) * source_term

! define location of the source
i = ISOURCE
j = JSOURCE

! interpolate density at the right location in the staggered grid cell
rho_half_x_half_y = 0.25d0 * (rho(i,j) + rho(i+1,j) + rho(i+1,j+1) + rho(i,j+1))

vx(i,j) = vx(i,j) + force_x * DELTAT / rho(i,j)
vy(i,j) = vy(i,j) + force_y * DELTAT / rho_half_x_half_y

! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
vx(1,:) = ZERO
vx(NX,:) = ZERO

vx(:,1) = ZERO
vx(:,NY) = ZERO

vy(1,:) = ZERO
vy(NX,:) = ZERO

vy(:,1) = ZERO
vy(:,NY) = ZERO

! compute total energy in the medium (without the PML layers)

! compute kinetic energy first, defined as 1/2 rho ||v||^2
! in principle we should use rho_half_x_half_y instead of rho for vy
! in order to interpolate density at the right location in the staggered grid cell
! but in a homogeneous medium we can safely ignore it
total_energy_kinetic(it) = 0.5d0 * sum( &
    rho(NPOINTS_PML+1:NX-NPOINTS_PML,NPOINTS_PML+1:NY-NPOINTS_PML)*( &
     vx(NPOINTS_PML+1:NX-NPOINTS_PML,NPOINTS_PML+1:NY-NPOINTS_PML)**2 +  &
     vy(NPOINTS_PML+1:NX-NPOINTS_PML,NPOINTS_PML+1:NY-NPOINTS_PML)**2))

! add potential energy, defined as 1/2 epsilon_ij sigma_ij
! in principle we should interpolate the medium parameters at the right location
! in the staggered grid cell but in a homogeneous medium we can safely ignore it
total_energy_potential(it) = ZERO
do j = NPOINTS_PML+1, NY-NPOINTS_PML
  do i = NPOINTS_PML+1, NX-NPOINTS_PML
    epsilon_xx = ((lambda(i,j) + 2.d0*mu(i,j)) * sigma_xx(i,j) - lambda(i,j) * &
      sigma_yy(i,j)) / (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)))
    epsilon_yy = ((lambda(i,j) + 2.d0*mu(i,j)) * sigma_yy(i,j) - lambda(i,j) * &
      sigma_xx(i,j)) / (4.d0 * mu(i,j) * (lambda(i,j) + mu(i,j)))
    epsilon_xy = sigma_xy(i,j) / (2.d0 * mu(i,j))
    total_energy_potential(it) = total_energy_potential(it) + &
      0.5d0 * (epsilon_xx * sigma_xx(i,j) + epsilon_yy * sigma_yy(i,j) + 2.d0 * epsilon_xy * sigma_xy(i,j))
  enddo
enddo

enddo   ! end of time loop

end program seismic_CPML_2D_iso_second