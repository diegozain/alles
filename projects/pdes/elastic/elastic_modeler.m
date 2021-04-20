clear
close all
% ------------------------------------------------------------------------------
%  2d elastic SV wave
% 
% rho*dt(v) = grad * T + F
% s*dt(u) = grad.' * v - tau*T
%
% where,
% 
% v = [vx ; vz]       particle velocity
% T = [sxx; sxz; szz] stress
% 
% grad = [dx 0  dz]
%        [0  dz dx]
% 
% inv(s) = [lam+2*mu lam        0]
%           [0         0        mu]
%           [lam       lam+2*mu  0]
% 
% ------------------------------------------------------------------------------
% vp = sqrt((lam + 2*mu)./rho);
% vs = sqrt(mu./rho);
% 
% lam = (vp.^2-2*vs.^2).*rho;
% mu  = rho.*vs^2;
% ------------------------------------------------------------------------------
% 
% the wave solver used here is identical to:
% 
% https://github.com/geodynamics/seismic_cpml/blob/master/seismic_CPML_2D_isotropic_second_order.f90,
% 
% but in Matlab & with free surface boundary conditions at the top.
% 
% ------------------------------------------------------------------------------
% --- pde parameters
load('outside-code/lam.mat');
load('outside-code/mu.mat');
load('outside-code/rho.mat');

vp = sqrt((lam + 2*mu)./rho);
vs = sqrt(mu./rho);
% ------------------------------------------------------------------------------
% --- spatial constraints
% -- slow - fast, fast - slow, homogeneous 1
X= 15;
Z= 6;
T= 0.42;
% -- source parameters
fo = 40;
to = 1/fo;
to = 5 * (to/pi);
wo = 2*pi*fo;

fo = 40;
tau = 1/fo;
tau= tau/pi;
to = 5 * tau;
wo = 2*pi*fo;
% ------------------------------------------------------------------------------
% --- stability
vp_min = min(vp(:));
vs_min = min(vs(:));

vp_max = max(vp(:));
vs_max = max(vs(:));

% taking min/max of vs AND vp
if vs_min==0
  vel_min = vp_min;
else
  vel_min = min([vs_min,vp_min]);
end
vel_max = max([vs_max,vp_max]);

% % Virieux, 1987 says this is not necessary,
% % and only vp is relevant. (he's wrong in general)
% vel_min = vp_min;
% vel_max = vp_max;

fmax = 2.2*fo;
% - space
l_min = vel_min / fmax;
no_p_wa = 10;
dx = l_min/no_p_wa;
dx=min([0.05; dx]);
dz = dx;
% - time
courant_factor = 0.3;
dt = 1/(vel_max * sqrt((1/dx^2)+(1/dz^2)));
dt = courant_factor * dt;
% ------------------------------------------------------------------------------
fprintf('\n  space discretization = %2.2d\n',dx)
fprintf('  time  discretization = %2.2d\n',dt)
fprintf('\n   2nd order in time, 4th order in space is \n ~~~ less stable than 2nd and 2nd order ~~~\n\n')
fprintf('  ...therefore, courant = %2.2d\n\n',courant_factor)
% ------------------------------------------------------------------------------
% --- discretization
x=(0:dx:X).';
z=(0:dz:Z).';
t=(0:dt:T).';
nx=numel(x);
nz=numel(z);
nt=numel(t);
% ------------------------------------------------------------------------------
fprintf('\n  end of length = %2.2d\n',x(nx))
fprintf('  end of depth  = %1.2d\n',z(nz))
fprintf('  end of time   = %2.2d\n\n',t(nt))

fprintf('  nx = %i\n',nx)
fprintf('  nz = %i\n',nz)
fprintf('  nt = %i\n\n',nt)

fprintf('  vp min = %2.2d\n',vp_min)
fprintf('  vp max = %2.2d\n',vp_max)

fprintf('  vs min = %2.2d\n',vs_min)
fprintf('  vs max = %2.2d\n\n',vs_max)
% ------------------------------------------------------------------------------
fprintf('wavecube of size (doubles): %2.2d Gb\n\n',nx*nz*nt*8*1e-9)
% ------------------------------------------------------------------------------
figure;
subplot(231)
fancy_imagesc(lam,x,z)
colormap(rainbow2_cb(1))
% colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Lame #1')
simple_figure()

subplot(232)
fancy_imagesc(mu,x,z)
colormap(rainbow2_cb(1))
% colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Lame #2')
simple_figure()

subplot(233)
fancy_imagesc(rho,x,z)
colormap(rainbow2_cb(1))
% colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Density')
simple_figure()

subplot(234)
fancy_imagesc(vp,x,z)
colormap(rainbow2_cb(1))
% colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('P velocity')
simple_figure()

subplot(236)
fancy_imagesc(vs,x,z)
colormap(rainbow2_cb(1))
% colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('S velocity')
simple_figure()
% ------------------------------------------------------------------------------
figure;
subplot(121)
fancy_imagesc(vp/fo,x,z)
xlabel('Length (m)')
ylabel('Depth (m)')
title('P wavelength (m)')
simple_figure()

subplot(122)
fancy_imagesc(vs/fo,x,z)
xlabel('Length (m)')
ylabel('Depth (m)')
title('S wavelength (m)')
simple_figure()
% ------------------------------------------------------------------------------
% -- receivers (real coordinates)
receivers = [x zeros(nx,1)];
% - receiver locations (index coordinates)
irecs_x = binning(x,receivers(:,1));
irecs_z = binning(x,receivers(:,2));
% ------------------------------------------------------------------------------
% -- source function (real coordinates)
src_xz = [x(fix(nx*0.15)) , z(1)]; 

% - source location (index coordinates)
isrc_x = binning(x,src_xz(1));
isrc_z = binning(z,src_xz(2));

% - sources in x and z
fx=zeros(nt,1);
fz=zeros(nt,1);

prompt = '\n    do you want ricker, hammer, or hammer delayed? (r, h, hd):  ';
src_type = input(prompt,'s');
if strcmp(src_type,'r')
  % ricker
  fz=( 1-0.5*(wo^2)*(t-to).^2 ) .* exp( -0.25*(wo^2)*(t-to).^2 );
end
if strcmp(src_type,'h')
  % 'hammer' according to the germans
  Fo = 1; % kg*m/s^2
  fz = Fo * (1/(dx*dz)) * sin((pi*t)./to).^3;
  fz(t>=to) = 0;
end
if strcmp(src_type,'hd')
  % 'hammer' according to the germans - but delayed at 0.5*to
  Fo = 1; % kg*m/s^2
  a=1;
  fz = Fo * (1/(dx*dz)) * sin((pi*(t-a*to))./to).^3;
  fz(t >= (1+a)*to) = 0;
  fz(t < a*to) = 0;
end

% NOTE: the staggered fd scheme actually outputs integral(f,dt)!!!
%       so, if you want an output source f, you need to input dt_(f,dt)
fx_=fx;
fz_=fz;
% NOTE: the hammer is special: Lisa Groos does not differentiate source (why??)
% NOTE: taking the numerical derivative of the hammer introduces lots of 
% higher freqs AND a bunch of notches in the freq domain, so better do the 
% analytical derivative for this one.
% Worth saying, even the analytical derivative has plenty of notches.
% derivative for the hammer
if ~strcmp(src_type(1),'h')
  fx = differentiate_o6(fx,dt);
  fz = differentiate_o6(fz,dt);
else
  fz = (1/(dx*dz)) * Fo*(3*pi/to)*(sin(pi*t/to).^2).*cos(pi*t/to);
  fz(t>=to) = 0;
end
% ------------------------------------------------------------------------------
% power spectra of source
figure;
subplot(121)
hold on
plot(t(1:fix(nt*0.5)),normali(fz_(1:fix(nt*0.5))),'r','linewidth',4)
plot(t(1:fix(nt*0.5)),normali(fz(1:fix(nt*0.5))),'k','linewidth',4)
hold off
axis tight
axis square
xlabel('Time (s)')
ylabel('Norm. amplitude')
set(gca,'ytick',[])
simple_figure()

[fz_fou,freq,~] = fourier_rt(fz,dt);
[fz_fou_,freq,~]= fourier_rt(fz_,dt);
fz_fou = normali(abs(fz_fou));
fz_fou_= normali(abs(fz_fou_));

subplot(122)
hold on
plot(freq(1:fix(numel(freq)*0.25)),fz_fou_(1:fix(numel(freq)*0.25)),'r','linewidth',4)
plot(freq(1:fix(numel(freq)*0.25)),fz_fou(1:fix(numel(freq)*0.25)),'k','linewidth',4)
hold off
axis tight
axis square
xlabel('Frequency (Hz)')
ylabel('Norm. power')
set(gca,'ytick',[])
legend({'source','dt(source)'})
simple_figure()
% ------------------------------------------------------------------------------
% -- init pml
n_points_pml= 20;
n_power_pml = 2; 
k_max_pml = 2; 
alpha_max_pml = 2*pi*(fo/2); 

thickness_PML_x = n_points_pml * dx;
thickness_PML_z = n_points_pml * dz;
Rcoef = 1e-3;
% ------------------------------------------------------------------------------
% -- free surface ghost nodes
n_ghost = 2;
% ------------------------------------------------------------------------------
% -- expand everything to pml & free surface
% lam, mu, rho, vp, vs, x, z, nx, nz

% x axis
lam= [repmat(lam(:,1),1,n_points_pml), lam , repmat(lam(:,nx),1,n_points_pml)];
% z axis
lam= [repmat(lam(1,:),n_ghost,1); lam ; repmat(lam(nz,:),n_points_pml,1)];
% x axis
mu = [repmat(mu(:,1),1,n_points_pml), mu , repmat(mu(:,nx),1,n_points_pml)];
% z axis
mu = [repmat(mu(1,:),n_ghost,1); mu ; repmat(mu(nz,:),n_points_pml,1)];
% x axis
rho= [repmat(rho(:,1),1,n_points_pml), rho , repmat(rho(:,nx),1,n_points_pml)];
% z axis
rho= [repmat(rho(1,:),n_ghost,1); rho ; repmat(rho(nz,:),n_points_pml,1)];

% % ?
% lam(1:n_ghost,:) = 0;
% mu(1:n_ghost,:)  = 0;
% rho(1:n_ghost,:) = 0;

x = [x; (dx*(1:2*n_points_pml)+x(nx)).' ];
z = [z; (dz*(1:(n_points_pml+n_ghost))+z(nz)).' ];

nx = nx+2*n_points_pml;
nz = nz + n_points_pml + n_ghost;

vp = sqrt((lam + 2*mu)./rho);
vs = sqrt(mu./rho);
% ------------------------------------------------------------------------------
% -- receivers new index because of pml & free surface
% - receoivers (index coordinates)
irecs_x= irecs_x+n_points_pml;
irecs_z= irecs_z+n_ghost;
irecs  = sub2ind([nz,nx],irecs_z,irecs_x);
nr = numel(irecs);
% ------------------------------------------------------------------------------
% -- source function new index because of pml & free surface
% - source location (index coordinates)
isrc_x = isrc_x+n_points_pml;
isrc_z = isrc_z+n_ghost;
% ------------------------------------------------------------------------------
% -- init fields

% - particle velocity
vx = zeros(nz,nx);
vz = zeros(nz,nx);

% - stress
sxx=zeros(nz,nx);
szz=zeros(nz,nx);
sxz=zeros(nz,nx);

% - data
d = zeros(nt,nr);
% ------------------------------------------------------------------------------
% - wavefield recorder
prompt = '\n    do you want to save the wave-cube? (y or n):  ';
wave_cube = input(prompt,'s');
if strcmp(wave_cube,'y')
 vz_=zeros(nz,nx,nt);
end
% ------------------------------------------------------------------------------
% -- init PML
a_x = zeros(nz,nx);
a_x_half = zeros(nz,nx);
b_x = zeros(nz,nx);
b_x_half = zeros(nz,nx);
K_x = zeros(nz,nx);
K_x_half = zeros(nz,nx);

a_z = zeros(nz,nx);
a_z_half = zeros(nz,nx);
b_z = zeros(nz,nx);
b_z_half = zeros(nz,nx);
K_z = zeros(nz,nx);
K_z_half = zeros(nz,nx);
% ------------------------------------------------------------------------------
%
%
%             pml init but requires computing (no defs anymore)
%
%
% ------------------------------------------------------------------------------
% % -- init temporal variables used in the static pml computation
% xval = 0;
% zval = 0;
% abscissa_in_PML = 0;
% abscissa_normalized = 0;

% d0_x = - (n_power_pml + 1) * vp * log(Rcoef) / (2 * thickness_PML_x); 
% d0_z = - (n_power_pml + 1) * vp * log(Rcoef) / (2 * thickness_PML_z);

% -------------
% x direction
% -------------
xoriginleft = thickness_PML_x;
xoriginright= (nx-1)*dx - thickness_PML_x;
for iz=1:nz
for ix=1:nx
  
  d0_x = - (n_power_pml + 1) * vp(iz,ix) * log(Rcoef) / (2 * thickness_PML_x); 
  
  xval = dx * (ix-1);
  
  % -- left edge
  % damping profile at the grid points
  abscissa_in_PML = xoriginleft - xval;
  if (abscissa_in_PML >= 0)
    abscissa_normalized = abscissa_in_PML / thickness_PML_x;
    d_x = d0_x * abscissa_normalized^n_power_pml;
    K_x(iz,ix) = 1 + (k_max_pml - 1) * abscissa_normalized^n_power_pml;
    alpha_x = alpha_max_pml * (1 - abscissa_normalized);
  end
  % damping profile at half the grid points
  abscissa_in_PML = xoriginleft - (xval + dx/2);
  if (abscissa_in_PML >= 0)
    abscissa_normalized = abscissa_in_PML / thickness_PML_x;
    d_x_half = d0_x * abscissa_normalized^n_power_pml;
    K_x_half(iz,ix) = 1 + (k_max_pml - 1) * abscissa_normalized^n_power_pml;
    alpha_x_half = alpha_max_pml * (1 - abscissa_normalized);
  end
  
  % - right edge
  % damping profile at the grid points
  abscissa_in_PML = xval - xoriginright;
  if (abscissa_in_PML >= 0)
    abscissa_normalized = abscissa_in_PML / thickness_PML_x;
    d_x = d0_x * abscissa_normalized^n_power_pml;
    K_x(iz,ix) = 1 + (k_max_pml - 1) * abscissa_normalized^n_power_pml;
    alpha_x = alpha_max_pml * (1 - abscissa_normalized);
  end
  % damping profile at half the grid points
  abscissa_in_PML = xval + dx/2 - xoriginright;
  if (abscissa_in_PML >= 0)
    abscissa_normalized = abscissa_in_PML / thickness_PML_x;
    d_x_half = d0_x * abscissa_normalized^n_power_pml;
    K_x_half(iz,ix) = 1 + (k_max_pml - 1) * abscissa_normalized^n_power_pml;
    alpha_x_half = alpha_max_pml * (1 - abscissa_normalized);
  end

  % 
  b_x(iz,ix) = exp(- ((d_x/(K_x(iz,ix))) + alpha_x) * dt);
  b_x_half(iz,ix) = exp(- ((d_x_half/(K_x_half(iz,ix))) + alpha_x_half) * dt);

  % this to avoid division by zero outside the PML
  if (abs(d_x) > 1e-6) 
    a_x(iz,ix) = d_x*(b_x(iz,ix) - 1)/(K_x(iz,ix)*(d_x + K_x(iz,ix) * alpha_x));
  end
  if (abs(d_x_half) > 1e-6) 
    a_x_half(iz,ix) = d_x_half * (b_x_half(iz,ix) - 1) / (K_x_half(iz,ix) * (d_x_half + K_x_half(iz,ix) * alpha_x_half));
  end
end
end

% -------------
% z direction
% -------------
zoriginbottom = thickness_PML_z;
zorigintop = (nz-1)*dz - thickness_PML_z;
for ix=1:nx
for iz=2:nz
  
  d0_z = - (n_power_pml + 1) * vp(iz,nx) * log(Rcoef) / (2 * thickness_PML_z);
  zval = dz * (iz-1);

  % -- bottom edge
  % damping profile at the grid points
  abscissa_in_PML = zoriginbottom - zval;
  if (abscissa_in_PML >= 0)
    abscissa_normalized = abscissa_in_PML / thickness_PML_z;
    d_z = d0_z * abscissa_normalized.^n_power_pml;
    K_z(iz,ix) = 1 + (k_max_pml - 1) * abscissa_normalized.^n_power_pml;
    alpha_z = alpha_max_pml * (1 - abscissa_normalized);
  end

  % damping profile at half the grid points
  abscissa_in_PML = zoriginbottom - (zval + dz/2);
  if (abscissa_in_PML >= 0)
    abscissa_normalized = abscissa_in_PML / thickness_PML_z;
    d_z_half = d0_z * abscissa_normalized.^n_power_pml;
    K_z_half(iz,ix) = 1 + (k_max_pml - 1) * abscissa_normalized.^n_power_pml;
    alpha_z_half = alpha_max_pml * (1 - abscissa_normalized);
  end

  % -- top edge
  % damping profile at the grid points
  abscissa_in_PML = zval - zorigintop;
  if (abscissa_in_PML >= 0)
    abscissa_normalized = abscissa_in_PML / thickness_PML_z;
    d_z = d0_z * abscissa_normalized.^n_power_pml;
    K_z(iz,ix) = 1 + (k_max_pml - 1) * abscissa_normalized.^n_power_pml;
    alpha_z = alpha_max_pml * (1 - abscissa_normalized);
  end
  
  % damping profile at half the grid points
  abscissa_in_PML = zval + dz/2 - zorigintop;
  if (abscissa_in_PML >= 0)
    abscissa_normalized = abscissa_in_PML / thickness_PML_z;
    d_z_half = d0_z * abscissa_normalized.^n_power_pml;
    K_z_half(iz,ix) = 1 + (k_max_pml - 1) * abscissa_normalized.^n_power_pml;
    alpha_z_half = alpha_max_pml * (1 - abscissa_normalized);
  end

  % just in case, for -5 at the end
  b_z(iz,ix) = exp(- ((d_z / (K_z(iz,ix))) + alpha_z) * dt);
  b_z_half(iz,ix) = exp(- ((d_z_half/(K_z_half(iz,ix))) + alpha_z_half) * dt);

  % this to avoid division by zero outside the PML
  if (abs(d_z) > 1e-6) 
    a_z(iz,ix) = d_z*(b_z(iz,ix) - 1)/(K_z(iz,ix)*(d_z + K_z(iz,ix)*alpha_z));
  end
  if (abs(d_z_half) > 1e-6) 
    a_z_half(iz,ix) = d_z_half * (b_z_half(iz,ix) - 1) / (K_z_half(iz,ix) * (d_z_half + K_z_half(iz,ix) * alpha_z_half));
  end
end
end
% these variables have some funky shit going on inside the domain,
% so this next part takes care of that.
% TODO: fix this entire section so the pml coefficients are not so funky.

b_x(:,(n_points_pml+2):(nx-n_points_pml-1)) = 0;
a_x(:,(n_points_pml+1):(nx-n_points_pml-1)) = 0;
K_x(:,(n_points_pml+2):(nx-n_points_pml-1)) = 1;
b_x_half(:,(n_points_pml+1):(nx-n_points_pml-1)) = 0;
a_x_half(:,(n_points_pml+1):(nx-n_points_pml-1)) = 0;
K_x_half(:,(n_points_pml+1):(nx-n_points_pml-1)) = 1;

b_z(1:(nz-n_points_pml-1),:) = 0;
a_z(1:(nz-n_points_pml-1),:) = 0;
K_z(1:(nz-n_points_pml-1),:) = 1;
b_z_half(1:(nz-n_points_pml-1),:) = 0;
a_z_half(1:(nz-n_points_pml-1),:) = 0;
K_z_half(1:(nz-n_points_pml-1),:) = 1;
% ------------------------------------------------------------------------------
% 
% 
%                   wave solver for loop in time
% 
% 
% ------------------------------------------------------------------------------
% % -- init temporal variables used in the time for loop
% lam_half_x = zeros(nz-1,nx-1);
% mu_half_x = zeros(nz-1,nx-1);
% lam2mu_half_x = zeros(nz-1,nx-1);
% mu_half_z = zeros(nz-1,nx-1);
% rho_half_x_half_z = 0;
% dvx_dx = zeros(nz-1,nx-1);
% dvx_dz = zeros(nz-1,nx-1);
% dvz_dx = zeros(nz-1,nx-1);
% dvz_dz = zeros(nz-1,nx-1);
% dsxx_dx = zeros(nz-1,nx-1);
% dszz_dz = zeros(nz-1,nx-1);
% dsxz_dx = zeros(nz-1,nx-1);
% dsxz_dz = zeros(nz-1,nx-1);

% could be saved only on PML region
dvx_dx_memory = zeros(nz,nx);
dvx_dz_memory = zeros(nz,nx);
dvz_dx_memory = zeros(nz,nx);
dvz_dz_memory = zeros(nz,nx);
dsxx_dx_memory= zeros(nz,nx);
dszz_dz_memory= zeros(nz,nx);
dsxz_dx_memory= zeros(nz,nx);
dsxz_dz_memory= zeros(nz,nx);
% ------------------------------------------------------------------------------
%
%
%                                wave solver
%
%
% ------------------------------------------------------------------------------
fprintf('\n\n ¡¡¡ solving the wave !!!\n\n');
tic;
for it=1:nt
  % ----------------------------------------------------------------------------
  %  print what time step we are in so we dont get anxious
  % ----------------------------------------------------------------------------
  if mod(it,10000)==0
    fprintf('we are at iteration %i of a total of %i \n',it,nt)
  end
  % ----------------------------------------------------------------------------
  %  compute stress sigma and update memory variables for C-PML
  % ----------------------------------------------------------------------------
  % compute dvx and dvz,
  % for sxx and szz
  iz=3:(nz-1);
  ix=2:(nx-2);
  % interpolate at the right location in the staggered grid cell
  lam_half_x= 0.5 * (lam(iz,ix+1) + lam(iz,ix));
  mu_half_x = 0.5 * (mu(iz,ix+1) + mu(iz,ix));
  lam2mu_half_x = lam_half_x + 2*mu_half_x;
  
  dvx_dx = (27*vx(iz,ix+1)-27*vx(iz,ix)-vx(iz,ix+2)+vx(iz,ix-1)) / (24*dx);
  dvz_dz = (27*vz(iz,ix)-27*vz(iz-1,ix)-vz(iz+1,ix)+vz(iz-2,ix)) / (24*dz);
  
  dvx_dx_memory(iz,ix) = b_x_half(iz,ix).*dvx_dx_memory(iz,ix) + a_x_half(iz,ix).*dvx_dx;
  dvz_dz_memory(iz,ix) = b_z(iz,ix).*dvz_dz_memory(iz,ix) + a_z(iz,ix).*dvz_dz;
  
  dvx_dx = (dvx_dx./K_x_half(iz,ix)) + dvx_dx_memory(iz,ix);
  dvz_dz = (dvz_dz./K_z(iz,ix)) + dvz_dz_memory(iz,ix);
  % ----------------------------------------------------------------------------
  % -- update sxx and szz
  sxx(iz,ix) = sxx(iz,ix) + (lam2mu_half_x.*dvx_dx + lam_half_x.*dvz_dz)*dt;
  szz(iz,ix) = szz(iz,ix) + (lam_half_x.*dvx_dx + lam2mu_half_x.*dvz_dz)*dt;
  % ----------------------------------------------------------------------------
  % -- free surface explicit conditions
  iz=n_ghost+1;
  % ix=2:(nx-2);
  i_ghost = 1:n_ghost;
  % -- boundary condition on sxx and szz
  % interpolate at the right location in the staggered grid cell
  lam_half_x= 0.5 * (lam(iz,ix+1) + lam(iz,ix));
  mu_half_x = 0.5 * (mu(iz,ix+1) + mu(iz,ix));
  lam2mu_half_x = lam_half_x + 2*mu_half_x;
  
  dvx_dx = (27*vx(iz,ix+1)-27*vx(iz,ix)-vx(iz,ix+2)+vx(iz,ix-1)) / (24*dx);
  dvx_dx_memory(iz,ix) = b_x_half(iz,ix).*dvx_dx_memory(iz,ix) + a_x_half(iz,ix).*dvx_dx;
  dvx_dx = (dvx_dx./K_x_half(iz,ix)) + dvx_dx_memory(iz,ix);

  % sxx(iz,ix) = sxx(iz,ix) + 4*(( (lam_half_x.*mu_half_x + mu_half_x.^2) ./ lam2mu_half_x ) .* dvx_dx)*dt;
  sxx(iz,ix) = sxx(iz,ix) + (lam2mu_half_x.*dvx_dx - lam_half_x.*dvx_dx)*dt;
  szz(iz,ix) = 0;
  
  szz(iz-i_ghost,ix)= -szz(iz+i_ghost,ix);
  % ----------------------------------------------------------------------------
  % compute dvx and dvz,
  % for sxz
  iz=2:(nz-2);
  ix=3:(nx-1);
  % interpolate at the right location in the staggered grid cell
  mu_half_z = 0.5 * (mu(iz+1,ix) + mu(iz,ix));

  dvz_dx = (27*vz(iz,ix)-27*vz(iz,ix-1)-vz(iz,ix+1)+vz(iz,ix-2)) / (24*dx);
  dvx_dz = (27*vx(iz+1,ix)-27*vx(iz,ix)-vx(iz+2,ix)+vx(iz-1,ix)) / (24*dz);
  
  dvz_dx_memory(iz,ix) = b_x(iz,ix).*dvz_dx_memory(iz,ix) + a_x(iz,ix).*dvz_dx;
  dvx_dz_memory(iz,ix) = b_z_half(iz,ix).* dvx_dz_memory(iz,ix) + a_z_half(iz,ix).* dvx_dz;
  
  dvz_dx = (dvz_dx./K_x(iz,ix)) + dvz_dx_memory(iz,ix);
  dvx_dz = (dvx_dz./K_z_half(iz,ix)) + dvx_dz_memory(iz,ix);
  % ----------------------------------------------------------------------------
  % -- update stress xz
  sxz(iz,ix) = sxz(iz,ix) + mu_half_z.*(dvz_dx + dvx_dz)*dt;
  % ----------------------------------------------------------------------------
  % -- free surface explicit conditions
  iz=n_ghost+1;
  % ix=1:nx;
  i_ghost = 1:n_ghost;
  % -- boundary condition on stress xz
  sxz(iz-i_ghost,ix) = -sxz(iz+i_ghost-1,ix);
  % ----------------------------------------------------------------------------
  %  compute velocity and update memory variables for C-PML
  % ----------------------------------------------------------------------------
  % compute dsxx and dsxz,
  % for vx
  iz=3:(nz-1);
  ix=3:(nx-1);
  
  dsxx_dx = (27*sxx(iz,ix)-27*sxx(iz,ix-1)-sxx(iz,ix+1)+sxx(iz,ix-2)) / (24*dx);
  dsxz_dz = (27*sxz(iz,ix)-27*sxz(iz-1,ix)-sxz(iz+1,ix)+sxz(iz-2,ix)) / (24*dz);
  
  dsxx_dx_memory(iz,ix) = b_x(iz,ix).*dsxx_dx_memory(iz,ix) + a_x(iz,ix).*dsxx_dx;
  dsxz_dz_memory(iz,ix) = b_z(iz,ix).*dsxz_dz_memory(iz,ix) + a_z(iz,ix).*dsxz_dz;
  
  dsxx_dx = (dsxx_dx./K_x(iz,ix)) + dsxx_dx_memory(iz,ix);
  dsxz_dz = (dsxz_dz./K_z(iz,ix)) + dsxz_dz_memory(iz,ix);
  % ----------------------------------------------------------------------------
  % -- update velocity vx
  vx(iz,ix) = vx(iz,ix) + (dsxx_dx+dsxz_dz)*dt./rho(iz,ix);
  % ----------------------------------------------------------------------------
  % -- free surface explicit conditions
  iz=n_ghost+1;
  % ix=1:nx;
  i_ghost = 1:n_ghost;
  % -- boundary condition on velocity vx
  vx(iz-i_ghost,ix) = 0;
  % ----------------------------------------------------------------------------
  % compute dsxz and dszz,
  % for vz
  iz=2:(nz-2);
  ix=2:(nx-2);
  % interpolate at the right location in the staggered grid cell
  rho_half_x_half_z = 0.25 * (rho(iz,ix) + rho(iz,ix+1) + rho(iz+1,ix+1) + rho(iz+1,ix));

  dsxz_dx = (27*sxz(iz,ix+1)-27*sxz(iz,ix)-sxz(iz,ix+2)+sxz(iz,ix-1)) / (24*dx);
  dszz_dz = (27*szz(iz+1,ix)-27*szz(iz,ix)-szz(iz+2,ix)+szz(iz-1,ix)) / (24*dz);
  
  dsxz_dx_memory(iz,ix) = (b_x_half(iz,ix).*dsxz_dx_memory(iz,ix)) + (a_x_half(iz,ix).*dsxz_dx);
  dszz_dz_memory(iz,ix) = (b_z_half(iz,ix).*dszz_dz_memory(iz,ix)) + (a_z_half(iz,ix).*dszz_dz);
  
  dsxz_dx = (dsxz_dx./K_x_half(iz,ix)) + dsxz_dx_memory(iz,ix);
  dszz_dz = (dszz_dz./K_z_half(iz,ix)) + dszz_dz_memory(iz,ix);
  % ----------------------------------------------------------------------------
  % -- update velocity vz
  vz(iz,ix) = vz(iz,ix) + (dsxz_dx + dszz_dz)*dt ./ rho_half_x_half_z;
  % ----------------------------------------------------------------------------
  % -- free surface explicit conditions
  iz=n_ghost+1;
  % ix=1:nx;
  i_ghost = 1:n_ghost;
  % -- boundary condition on velocity vz
  vz(iz-i_ghost,ix) = 0;
  % ----------------------------------------------------------------------------
  %  source update
  % ----------------------------------------------------------------------------
  % interpolate density at the right location in the staggered grid cell
  rho_half_x_half_z = 0.25 * (rho(isrc_z,isrc_x) + rho(isrc_z,isrc_x+1) + rho(isrc_z+1,isrc_x+1) + rho(isrc_z+1,isrc_x));
  
  vx(isrc_z,isrc_x) = vx(isrc_z,isrc_x) + fx(it)*dt/rho(isrc_z,isrc_x);
  vz(isrc_z,isrc_x) = vz(isrc_z,isrc_x) + fz(it)*dt/rho_half_x_half_z;
  % ----------------------------------------------------------------------------
  % % source as a boundary condition
  % vx(isrc_z,isrc_x) = fx(it);
  % vz(isrc_z,isrc_x) = fz(it);
  % ----------------------------------------------------------------------------
  %  store data
  % ----------------------------------------------------------------------------
  % d(it,:) = vz(n_ghost+1,:); % receivers
  d(it,:) = vz(irecs); % receivers
  % ----------------------------------------------------------------------------
  %  store wavefield
  % ----------------------------------------------------------------------------
  if strcmp(wave_cube,'y')
    vz_(:,:,it) = vz;
  end
end
toc;
% ------------------------------------------------------------------------------
if strcmp(wave_cube,'y')
  vz_min = min(vz_(:));
  vz_max = max(vz_(:));
  vz_min = max([abs(vz_min) vz_max]);
  vz_min = vz_min*0.03;

  t1=to*2.25;
  t2=to*3;
  t3=to*5;

  figure;

  subplot(2,3,1)
  fancy_imagesc(vz_(:,:,binning(t,t1)),x,z)
  colormap(rainbow2_cb(1))
  caxis([-vz_min vz_min])
  hold on
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'k.','markersize',60)
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'w.','markersize',40)
  hold off
  colorbar off
  % xlabel('Length')
  % ylabel('Depth')
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  % title('First')
  simple_figure()

  subplot(2,3,2)
  fancy_imagesc(vz_(:,:,binning(t,t2)),x,z)
  colormap(rainbow2_cb(1))
  caxis([-vz_min vz_min])
  hold on
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'k.','markersize',60)
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'w.','markersize',40)
  hold off
  colorbar off
  % xlabel('Length')
  % ylabel('Depth')
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  title('Wavefield snapshots')
  % title('Middle')
  simple_figure()

  subplot(2,3,3)
  fancy_imagesc(vz_(:,:,binning(t,t3)),x,z)
  colormap(rainbow2_cb(1))
  caxis([-vz_min vz_min]);
  hold on
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'k.','markersize',60)
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'w.','markersize',40)
  hold off
  colorbar off
  % xlabel('Length')
  % ylabel('Depth')
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  % title('Final')
  simple_figure()

  subplot(2,3,[4 5 6])
  hold on
  plot(t,fz_,'k','linewidth',2)
  plot(t1*ones(3,1),linspace(min(fz_),max(fz_),3),'linewidth',3)
  plot(t2*ones(3,1),linspace(min(fz_),max(fz_),3),'linewidth',3)
  plot(t3*ones(3,1),linspace(min(fz_),max(fz_),3),'linewidth',3)
  hold off
  axis tight
  set(gca,'ytick',[])
  xlabel('Time')
  title('Source')
  simple_figure()
end
% ------------------------------------------------------------------------------
d_min = min(d(:));
d_max = max(d(:));
d_min = max([abs(d_min) d_max]);
d_min = d_min*0.8;

figure;
fancy_imagesc(d,x,t);
axis normal;
colorbar off
caxis(1e-1*[-d_min d_min])
set(gca,'xtick',[])
set(gca,'ytick',[])
% ylabel('Time')
% xlabel('Length')
title('Data')
simple_figure()
% ------------------------------------------------------------------------------
figure;

subplot(231)
fancy_imagesc(sxx,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Stress xx')
simple_figure()

subplot(232)
fancy_imagesc(szz,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Stress zz')
simple_figure()

subplot(233)
fancy_imagesc(sxz,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Stress xz')
simple_figure()

subplot(234)
fancy_imagesc(vx,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Velocity x')
simple_figure()

subplot(236)
fancy_imagesc(vz,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_ghost,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Velocity z')
simple_figure()
% ------------------------------------------------------------------------------
figure;
fancy_imagesc(d,x,t)
axis normal
colorbar off
caxis(1e-1*[-d_min d_min])
xlabel('Length (m)')
ylabel('Time (s)')
title('Surface receivers')
simple_figure()
% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
fancy_imagesc(vp,x,z)
colormap(rainbow2_cb(1))
xlabel('Length (m)')
ylabel('Depth (m)')
title('P velocity (m/s)')
simple_figure()

subplot(1,2,2)
fancy_imagesc(vs,x,z)
colormap(rainbow2_cb(1))
xlabel('Length (m)')
ylabel('Depth (m)')
title('S velocity (m/s)')
simple_figure()
% ------------------------------------------------------------------------------
%{
% linear semblance
d_onesided = d(:,(isrc_x-n_points_pml):(nx-2*n_points_pml));
rx= x((isrc_x-n_points_pml):(nx-2*n_points_pml)) - x(isrc_x);

velos = linspace(vel_min*0.8,vel_max*1.2,1e+2);

v_analy = v_linear(d_onesided,t,rx,velos);
semblance_v = sum(v_analy,1);

[~,iv] = max(semblance_v);
fprintf('\n  max linear velocity achieved at %2.2d m/s',velos(iv))
fprintf('\n  this velocity is %2.2d percent of Vs\n\n',velos(iv)*100/vs(n_ghost+1))
% ------------------------------------------------------------------------------
figure;
subplot(2,2,1)
fancy_imagesc(v_analy,velos,t);
colorbar off
axis normal
title('Linear semblance')
xlabel('Velocity (m/s)')
ylabel('Time (s)')
simple_figure()

subplot(2,2,3)
plot(velos,semblance_v)
set(gca,'ytick',[])
axis tight
xlabel('Velocity (m/s)')
ylabel('Lin. semblance')
simple_figure()

subplot(2,2,[2,4])
fancy_imagesc(d_onesided,rx,t);
colorbar off
caxis(1e-1*[-d_min d_min])
axis normal
title('One-sided data') 
xlabel('Receivers (m)')
ylabel('Time (s)')
simple_figure()
% ------------------------------------------------------------------------------
% % dispersion analysis
% [d_onesided_,f,df] = fourier_rt(d_onesided,dt);
% 
% f_disp = 0:df:30;
% dsx= (1/vel_min) - (1/vel_max);
% dsx= dsx/1e4;
% sx = (1/vel_max):dsx:(1/vel_min);
% 
% [disper_vxf,disper_sxf] = masw(d_onesided_,rx,sx,f,f_disp);
% 
% figure;
% fancy_imagesc(flip(disper_vxf,1),f,flip(1./sx));
% colormap(rainbow2_cb(1))
% axis normal
% xlabel('Frequency (Hz)')
% ylabel('Phase velocity (m/s)')
% title('MASW')
% simple_figure()
% ------------------------------------------------------------------------------
%}
