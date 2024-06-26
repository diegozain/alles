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
% but in Matlab.
% 
% ------------------------------------------------------------------------------
% - air at 20C and 1 atm (101.325 kPa): 
% lam= 1.4e+5 Kg/m/s^2
% mu = 0
% rho= 1.2041 Kg/m^3
%
% yet another velocity for air:
% vp = 300 m/s;
% vs = 0.0 m/s;
% rho = 1.25 kg/m^3
% lam = 112500 kg/m^3
% mu = 0 kg/m^3
% 
% - some rock
% lam= 1.0163e+10
% mu = 1.0165e+10
% rho= 1.7e+3
% ------------------------------------------------------------------------------
% --- pde parameters
% lam_= [2; 5];
% mu_ = [3; 9];
% rho_= [1; 6];

lam_= [2; 2];
mu_ = [3; 3];
rho_= [6; 6];

% lam_= [2; 2];
% mu_ = [0; 3];
% rho_= [6; 6];

% lam_= [1e+5; 1e+9];
% mu_ = [0; 1e+10];
% rho_= [1; 2e+3];
% --- spatial constraints
X=3;
Z=2.5;
T=5;
% --- source parameters
to = 0.5;
fo = 3;
wo = 2*pi*fo;
% ------------------------------------------------------------------------------
% --- stability
vp = sqrt((lam_ + 2*mu_)./rho_);
vs = sqrt(mu_./rho_);

vs_min = min(vs);
vp_min = min(vp);
vs_max = max(vs);
vp_max = max(vp);

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
no_p_wa = 10; % 10 50
dx = l_min/no_p_wa;
dx=min([0.1; dx]);
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
fprintf('wavecube of size: %2.2d Gb\n\n',nx*nz*nt*8*1e-9)
% ------------------------------------------------------------------------------
% --- pde parameters
lam=lam_(1)*ones(nz,nx);
mu =mu_(1)*ones(nz,nx);
rho=rho_(1)*ones(nz,nx);

% % -- box in the middle
% lam(fix(nz*(1/3)):fix(nz*(2/3)),fix(nx*(1/3)):fix(nx*(2/3)))= lam_(2);
% mu(fix(nz*(1/3)):fix(nz*(2/3)),fix(nx*(1/3)):fix(nx*(2/3))) = mu_(2);
% rho(fix(nz*(1/3)):fix(nz*(2/3)),fix(nx*(1/3)):fix(nx*(2/3)))= rho_(2);

% -- two layers
lam(fix(nz*(1/4)):nz,:)= lam_(2);
mu(fix(nz*(1/4)):nz,:) = mu_(2);
rho(fix(nz*(1/4)):nz,:)= rho_(2);

vp = sqrt((lam + 2*mu)./rho);
vs = sqrt(mu./rho);
% ------------------------------------------------------------------------------
figure;
subplot(331)
fancy_imagesc(lam,x,z)
colormap(rainbow2_cb(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Lame #1')
simple_figure()

subplot(333)
fancy_imagesc(mu,x,z)
colormap(rainbow2_cb(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Lame #2')
simple_figure()

subplot(335)
fancy_imagesc(rho,x,z)
colormap(rainbow2_cb(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Density')
simple_figure()

subplot(337)
fancy_imagesc(vp,x,z)
colormap(rainbow2_cb(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('P velocity')
simple_figure()

subplot(339)
fancy_imagesc(vs,x,z)
colormap(rainbow2_cb(1))
colorbar off
% xlabel('Length')
% ylabel('Depth')
set(gca,'xtick',[])
set(gca,'ytick',[])
title('S velocity')
simple_figure()
% ------------------------------------------------------------------------------
% -- source function (real coordinates)

% - source location 
% src_xz = [(x(end)+x(1))/2 , (z(end)+z(1))*(1/2)];
% src_xz = [(x(end)+x(1))/2 , z(60)];
% src_xz = [x(fix(nx*0.5)) , z(fix(nz*0.5))];
src_xz = [x(fix(nx*0.5)) , z(fix(nz*(1/4)) - 1)];

% - source location (index coordinates)
isrc_x = binning(x,src_xz(1));
isrc_z = binning(z,src_xz(2));

% - sources in x and z
fx=zeros(nt,1);
fz=zeros(nt,1);
% ricker
fz=( 1-0.5*(wo^2)*(t-to).^2 ) .* exp( -0.25*(wo^2)*(t-to).^2 );
% % 'hammer' according to the germans
% Fo = 1; % kg*m/s^2
% fz = Fo * (1/(dx*dz)) * sin((pi*t)./(to)).^3;
% fz(t>=to) = 0;
% % 'hammer' according to the germans - but delayed at 0.5*to
% fz = Fo * (1/(dx*dz)) * sin((pi*(t-0.5*to))./(to)).^3;
% fz(t >= 1.5*to) = 0;
% fz(t < 0.5*to) = 0;
% NOTE: the staggered fd scheme actually outputs integral(f,dt)!!!
%       so, if you want an output source f, you need to input dt_(f,dt)
fx_=fx;
fz_=fz;
fx = differentiate_line(fx,dt);
fz = differentiate_line(fz,dt);
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
n_points_pml= 10;
n_power_pml = 2; % 2
k_max_pml = 2; % 1 80
alpha_max_pml = 2*pi*(fo/2); % *1

thickness_PML_x = n_points_pml * dx;
thickness_PML_z = n_points_pml * dz;
Rcoef = 1e-3; % 1e-3
% ------------------------------------------------------------------------------
% -- expand everything to pml
% lam, mu, rho, vp, vs, x, z, nx, nz
lam= [repmat(lam(:,1),1,n_points_pml), lam , repmat(lam(:,nx),1,n_points_pml)];
lam= [repmat(lam(1,:),n_points_pml,1); lam ; repmat(lam(nz,:),n_points_pml,1)];
mu = [repmat(mu(:,1),1,n_points_pml), mu , repmat(mu(:,nx),1,n_points_pml)];
mu = [repmat(mu(1,:),n_points_pml,1); mu ; repmat(mu(nz,:),n_points_pml,1)];
rho= [repmat(rho(:,1),1,n_points_pml), rho , repmat(rho(:,nx),1,n_points_pml)];
rho= [repmat(rho(1,:),n_points_pml,1); rho ; repmat(rho(nz,:),n_points_pml,1)];
vp = sqrt((lam + 2*mu)./rho);
vs = sqrt(mu./rho);

x = [x; (dx*(1:2*n_points_pml)+x(nx)).' ];
z = [z; (dz*(1:2*n_points_pml)+z(nz)).' ];

nx = nx+2*n_points_pml;
nz = nz+2*n_points_pml;
% ------------------------------------------------------------------------------
% -- source function new index because of pml
% - source location (index coordinates)
isrc_x = isrc_x+n_points_pml;
isrc_z = isrc_z+n_points_pml;
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
d = zeros(nt,nx);

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
for iz=1:nz
  
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

b_z((n_points_pml+2):(nz-n_points_pml-1),:) = 0;
a_z((n_points_pml+1):(nz-n_points_pml-1),:) = 0;
K_z((n_points_pml+2):(nz-n_points_pml-1),:) = 1;
b_z_half((n_points_pml+1):(nz-n_points_pml-1),:) = 0;
a_z_half((n_points_pml+1):(nz-n_points_pml-1),:) = 0;
K_z_half((n_points_pml+1):(nz-n_points_pml-1),:) = 1;
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
% value_dvx_dx = zeros(nz-1,nx-1);
% value_dvx_dz = zeros(nz-1,nx-1);
% value_dvz_dx = zeros(nz-1,nx-1);
% value_dvz_dz = zeros(nz-1,nx-1);
% value_dsxx_dx = zeros(nz-1,nx-1);
% value_dszz_dz = zeros(nz-1,nx-1);
% value_dsxz_dx = zeros(nz-1,nx-1);
% value_dsxz_dz = zeros(nz-1,nx-1);

% could be saved only on PML region
memory_dvx_dx = zeros(nz,nx);
memory_dvx_dz = zeros(nz,nx);
memory_dvz_dx = zeros(nz,nx);
memory_dvz_dz = zeros(nz,nx);
memory_dsxx_dx = zeros(nz,nx);
memory_dszz_dz = zeros(nz,nx);
memory_dsxz_dx = zeros(nz,nx);
memory_dsxz_dz = zeros(nz,nx);
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
  
  value_dvx_dx = (27*vx(iz,ix+1)-27*vx(iz,ix)-vx(iz,ix+2)+vx(iz,ix-1)) / (24*dx);
  value_dvz_dz = (27*vz(iz,ix)-27*vz(iz-1,ix)-vz(iz+1,ix)+vz(iz-2,ix)) / (24*dz);
  
  memory_dvx_dx(iz,ix) = b_x_half(iz,ix).*memory_dvx_dx(iz,ix) + a_x_half(iz,ix).*value_dvx_dx;
  memory_dvz_dz(iz,ix) = b_z(iz,ix).*memory_dvz_dz(iz,ix) + a_z(iz,ix).*value_dvz_dz;
  
  value_dvx_dx = (value_dvx_dx./K_x_half(iz,ix)) + memory_dvx_dx(iz,ix);
  value_dvz_dz = (value_dvz_dz./K_z(iz,ix)) + memory_dvz_dz(iz,ix);
  % ----------------------------------------------------------------------------
  % -- update sxx and szz
  sxx(iz,ix) = sxx(iz,ix) + (lam2mu_half_x.*value_dvx_dx + lam_half_x.*value_dvz_dz)*dt;
  
  szz(iz,ix) = szz(iz,ix) + (lam_half_x.*value_dvx_dx + lam2mu_half_x.*value_dvz_dz)*dt;
  % ----------------------------------------------------------------------------
  % compute dvx and dvz,
  % for sxz
  iz=2:(nz-2);
  ix=3:(nx-1);
  % interpolate at the right location in the staggered grid cell
  mu_half_z = 0.5 * (mu(iz+1,ix) + mu(iz,ix));

  value_dvz_dx = (27*vz(iz,ix)-27*vz(iz,ix-1)-vz(iz,ix+1)+vz(iz,ix-2)) / (24*dx);
  value_dvx_dz = (27*vx(iz+1,ix)-27*vx(iz,ix)-vx(iz+2,ix)+vx(iz-1,ix)) / (24*dz);
  
  memory_dvz_dx(iz,ix) = b_x(iz,ix).*memory_dvz_dx(iz,ix) + a_x(iz,ix).*value_dvz_dx;
  memory_dvx_dz(iz,ix) = b_z_half(iz,ix).* memory_dvx_dz(iz,ix) + a_z_half(iz,ix).* value_dvx_dz;
  
  value_dvz_dx = (value_dvz_dx./K_x(iz,ix)) + memory_dvz_dx(iz,ix);
  value_dvx_dz = (value_dvx_dz./K_z_half(iz,ix)) + memory_dvx_dz(iz,ix);
  % ----------------------------------------------------------------------------
  % -- update stress xz
  sxz(iz,ix) = sxz(iz,ix) + mu_half_z.*(value_dvz_dx + value_dvx_dz)*dt;
  % ----------------------------------------------------------------------------
  %  compute velocity and update memory variables for C-PML
  % ----------------------------------------------------------------------------
  % compute dsxx and dsxz,
  % for vx
  iz=3:(nz-1);
  ix=3:(nx-1);
  
  value_dsxx_dx = (27*sxx(iz,ix)-27*sxx(iz,ix-1)-sxx(iz,ix+1)+sxx(iz,ix-2)) / (24*dx);
  value_dsxz_dz = (27*sxz(iz,ix)-27*sxz(iz-1,ix)-sxz(iz+1,ix)+sxz(iz-2,ix)) / (24*dz);
  
  memory_dsxx_dx(iz,ix) = b_x(iz,ix).*memory_dsxx_dx(iz,ix) + a_x(iz,ix).*value_dsxx_dx;
  memory_dsxz_dz(iz,ix) = b_z(iz,ix).*memory_dsxz_dz(iz,ix) + a_z(iz,ix).*value_dsxz_dz;
  
  value_dsxx_dx = (value_dsxx_dx./K_x(iz,ix)) + memory_dsxx_dx(iz,ix);
  value_dsxz_dz = (value_dsxz_dz./K_z(iz,ix)) + memory_dsxz_dz(iz,ix);
  % ----------------------------------------------------------------------------
  % -- update velocity vx
  vx(iz,ix) = vx(iz,ix) + (value_dsxx_dx+value_dsxz_dz)*dt./rho(iz,ix);
  % ----------------------------------------------------------------------------
  % compute dsxz and dszz,
  % for vz
  iz = 2:(nz-2);
  ix = 2:(nx-2);
  % interpolate at the right location in the staggered grid cell
  rho_half_x_half_z = 0.25 * (rho(iz,ix) + rho(iz,ix+1) + rho(iz+1,ix+1) + rho(iz+1,ix));

  value_dsxz_dx = (27*sxz(iz,ix+1)-27*sxz(iz,ix)-sxz(iz,ix+2)+sxz(iz,ix-1)) / (24*dx);
  value_dszz_dz = (27*szz(iz+1,ix)-27*szz(iz,ix)-szz(iz+2,ix)+szz(iz-1,ix)) / (24*dz);
  
  memory_dsxz_dx(iz,ix) = (b_x_half(iz,ix).*memory_dsxz_dx(iz,ix)) + (a_x_half(iz,ix).*value_dsxz_dx);
  memory_dszz_dz(iz,ix) = (b_z_half(iz,ix).*memory_dszz_dz(iz,ix)) + (a_z_half(iz,ix).*value_dszz_dz);
  
  value_dsxz_dx = (value_dsxz_dx./K_x_half(iz,ix)) + memory_dsxz_dx(iz,ix);
  value_dszz_dz = (value_dszz_dz./K_z_half(iz,ix)) + memory_dszz_dz(iz,ix);
  % ----------------------------------------------------------------------------
  % -- update velocity vz
  vz(iz,ix) = vz(iz,ix) + (value_dsxz_dx + value_dszz_dz)*dt ./ rho_half_x_half_z;
  % ----------------------------------------------------------------------------
  %  source update
  % ----------------------------------------------------------------------------
  % interpolate density at the right location in the staggered grid cell
  rho_half_x_half_z = 0.25 * (rho(isrc_z,isrc_x) + rho(isrc_z,isrc_x+1) + rho(isrc_z+1,isrc_x+1) + rho(isrc_z+1,isrc_x));
  
  vx(isrc_z,isrc_x) = vx(isrc_z,isrc_x) + fx(it)*dt/rho(isrc_z,isrc_x);
  vz(isrc_z,isrc_x) = vz(isrc_z,isrc_x) + fz(it)*dt/rho_half_x_half_z;
  % ----------------------------------------------------------------------------
  %  wrap up dirichlet bc on edges
  % ----------------------------------------------------------------------------
  vx(:,1) = 0;
  vx(:,nx)= 0;

  vx(1,:) = 0;
  vx(nz,:)= 0;

  vz(:,1) = 0;
  vz(:,nx)= 0;

  vz(1,:) = 0;
  vz(nz,:)= 0;
  % ----------------------------------------------------------------------------
  %  store data
  % ----------------------------------------------------------------------------
  d(it,:) = vz(n_points_pml+2+isrc_z,:);
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

  t1=to*1.25;
  t2=to*1.65;
  t3=to*2.25;

  figure;

  subplot(231)
  fancy_imagesc(vz_(:,:,binning(t,t1)),x,z)
  colormap(rainbow2_cb(1))
  caxis([-vz_min vz_min])
  hold on
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'k.','markersize',60)
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'w.','markersize',40)
  hold off
  colorbar off
  % xlabel('Length')
  % ylabel('Depth')
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  % title('First')
  simple_figure()

  subplot(232)
  fancy_imagesc(vz_(:,:,binning(t,t2)),x,z)
  colormap(rainbow2_cb(1))
  caxis([-vz_min vz_min])
  hold on
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'k.','markersize',60)
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'w.','markersize',40)
  hold off
  colorbar off
  % xlabel('Length')
  % ylabel('Depth')
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  title('Wavefield snapshots')
  % title('Middle')
  simple_figure()

  subplot(233)
  fancy_imagesc(vz_(:,:,binning(t,t3)),x,z)
  colormap(rainbow2_cb(1))
  caxis([-vz_min vz_min]);
  hold on
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'k.','markersize',60)
  plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'w.','markersize',40)
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

subplot(331)
fancy_imagesc(sxx,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Stress xx')
simple_figure()

subplot(333)
fancy_imagesc(szz,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Stress zz')
simple_figure()

subplot(335)
fancy_imagesc(sxz,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Stress xz')
simple_figure()

subplot(337)
fancy_imagesc(vx,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Velocity x')
simple_figure()

subplot(339)
fancy_imagesc(vz,x,z)
colormap(rainbow2_cb(1))
hold on
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'k.','markersize',60)
plot(src_xz(1)+dx*n_points_pml,src_xz(2)+dz*n_points_pml,'w.','markersize',40)
hold off
colorbar off
% xlabel('Length')
% ylabel('Depth')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
title('Velocity z')
simple_figure()
% ------------------------------------------------------------------------------
d_onesided = d(:,isrc_x:(nx-n_points_pml));                                               
rx= x(isrc_x:(nx-n_points_pml)) - x(isrc_x);                                              
figure;hold on;plot(t,d_onesided(:,10));plot(t,d_onesided(:,50));plot(t,d_onesided(:,100))
%}