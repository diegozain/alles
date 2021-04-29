clear
close all
% ------------------------------------------------------------------------------
nx=200;
ny=200;
% ------------------------------------------------------------------------------
y=linspace(-1, 1,ny);
x=linspace(-1, 1,nx);

dx=x(2)-x(1);
dy=y(2)-y(1);
% ------------------------------------------------------------------------------
[Dy,Dx] = Dx_Dz(ny,nx);
Dx=(1/dx)*Dx;
Dy=(1/dy)*Dy;
% ------------------------------------------------------------------------------
[xx,yy] = meshgrid(x,y);

ux =   sin(pi*xx).*cos(pi*yy) + sin(pi*yy).*cos(pi*xx) + 0.5*ones(ny,nx);
% uy = - sin(pi*yy).*cos(pi*xx) + sin(pi*xx).*cos(pi*yy) - ones(ny,nx);
uy =   -sin(pi*xx).*cos(pi*yy) - sin(pi*yy).*cos(pi*xx) + 0.5*ones(ny,nx);

% ------------------------------------------------------------------------------
figure;

subplot(1,2,1)
fancy_imagesc(ux,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('True ux')
simple_figure()

subplot(1,2,2)
fancy_imagesc(uy,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('True uy')
simple_figure()
% ------------------------------------------------------------------------------
ux = ux(:);
uy = uy(:);
% ------------------------------------------------------------------------------
uxx = Dx*ux;
uyy = Dy*uy;

uxy = Dy*ux;
uyx = Dx*uy;

divu = uxx + uyy;
rotu = uyx - uxy;
% ------------------------------------------------------------------------------
% solve
% ------------------------------------------------------------------------------
% 
% unfourtantely, this one does not work :(
% 
% L = Dx*Dx + Dy*Dy;
% 
% phi_ = L \ divu;
% psi_ = L \ rotu;
% ------------------------------------------------------------------------------
phi_=zeros(ny,nx);
psi_=zeros(ny,nx);

for ix_=1:nx
 for iy_=1:ny
  xo=x(ix_);
  yo=y(iy_);
  g = - (1/2/pi)*log( sqrt( ( xx-xo ).^2 + ( yy-yo ).^2 ) );

  g(iy_,ix_) = 0;

  a = g(:).*divu;
  phi_(iy_,ix_) = sum(a)*dx*dy;

  b = g(:).*rotu;
  psi_(iy_,ix_) = sum(b)*dx*dy;
 end
end

phi_=phi_(:);
psi_=psi_(:);
% ------------------------------------------------------------------------------
% rebuild
phi_x = Dx*phi_;
phi_y = Dy*phi_;

psi_x = Dx*psi_;
psi_y = Dy*psi_;

hx = ux - phi_x - psi_y;
hy = uy - phi_y + psi_x;

ux = phi_x + psi_y + hx;
uy = phi_y - psi_x + hy;
% ------------------------------------------------------------------------------
hxx = Dx*hx;
hyy = Dy*hy;

hxy = Dy*hx;
hyx = Dx*hy;

divh = hxx + hyy;
roth = hyx - hxy;
% ------------------------------------------------------------------------------
ux=reshape(ux,ny,nx);
uy=reshape(uy,ny,nx);

phi_=reshape(phi_,ny,nx);
psi_=reshape(psi_,ny,nx);

phi_x=reshape(phi_x,ny,nx);
phi_y=reshape(phi_y,ny,nx);

psi_x=reshape(psi_x,ny,nx);
psi_y=reshape(psi_y,ny,nx);

hx=reshape(hx,ny,nx);
hy=reshape(hy,ny,nx);
% ------------------------------------------------------------------------------
figure;

subplot(4,2,1)
fancy_imagesc(ux,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('ux')
simple_figure()

subplot(4,2,2)
fancy_imagesc(uy,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('uy')
simple_figure()

subplot(4,2,3)
fancy_imagesc(phi_x,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Curl-free x')
simple_figure()

subplot(4,2,4)
fancy_imagesc(phi_y,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Curl-free y')
simple_figure()

subplot(4,2,5)
fancy_imagesc(psi_y,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Div-free x')
simple_figure()

subplot(4,2,6)
fancy_imagesc(-psi_x,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Div-free y')
simple_figure()

subplot(4,2,7)
fancy_imagesc(hx,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Harmonic x')
simple_figure()

subplot(4,2,8)
fancy_imagesc(hy,x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Harmonic y')
simple_figure()
% ------------------------------------------------------------------------------
figure;

subplot(2,4,1)
fancy_imagesc(ux(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,2)
fancy_imagesc(phi_x(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,3)
fancy_imagesc(psi_y(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,4)
fancy_imagesc(hx(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,5)
fancy_imagesc(uy(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,6)
fancy_imagesc(phi_y(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,7)
fancy_imagesc(-psi_x(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()

subplot(2,4,8)
fancy_imagesc(hy(3:(ny-3),3:(nx-3)),x,y);
colorbar('off')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
simple_figure()
% ------------------------------------------------------------------------------
% % }