close all
clear
clc
% ------------------------------------------------------------------------------
x=linspace(0,20,100);
z=linspace(0,10,50);
% ------------------------------------------------------------------------------
ar=0.5; % 1.5
ac=[binning(z,z(end)*0.45),binning(x,x(end)*0.5)];
a = gaussi(x,z,ac,ar,2);
% ------------------------------------------------------------------------------
box_ix_ = 8;
box_ix__= 12;
box_iz_ = 3;
box_iz__= 7;
% ------------------------------------------------------------------------------
xbox = [box_ix_; box_ix__; box_ix__; box_ix_; box_ix_];
zbox = [box_iz_; box_iz_; box_iz__; box_iz__; box_iz_];
box_ix_  = binning(x,box_ix_);
box_ix__ = binning(x,box_ix__);
box_iz_  = binning(z,box_iz_);
box_iz__ = binning(x,box_iz__);
b = zeros(numel(z),numel(x));
b(box_iz_:box_iz__,box_ix_:box_ix__) = 1;
% ------------------------------------------------------------------------------
thet = linspace(0,2*pi,100);
circ_x = cos(thet);
circ_z = sin(thet);
circ_ax = 1.5*ar*circ_x+x(ac(2));
circ_az = 1.5*ar*circ_z+z(ac(1));
% ------------------------------------------------------------------------------
% compute xgradient
[nz,nx] = size(a);
[Dz,Dx] = Dx_Dz(nz,nx);
xab = cross2d(a,b);
% ------------------------------------------------------------------------------
figure;
hold on
fancy_imagesc(xab,x,z)
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
title('\tau');
simple_figure();
% ------------------------------------------------------------------------------
Ja = crossJa(b);
Jb = crossJb(a);
J=Ja+Jb;
% ------------------------------------------------------------------------------
ga=Ja*xab(:);
ga = normali(ga);
ga=reshape(ga,size(a));
% ------------------------------------------------------------------------------
gb=Jb*xab(:);
gb = normali(gb);
gb=reshape(gb,size(a));
% ------------------------------------------------------------------------------
figure;
subplot(2,1,1)
hold on
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
fancy_imagesc(ga,x,z)
hold on
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
set(gca,'YTick',[])
set(gca,'XTick',[])
title('gradient \tau_a')
simple_figure()

subplot(2,1,2)
hold on
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
fancy_imagesc(gb,x,z)
hold on
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
title('gradient \tau_b')
simple_figure()
% ------------------------------------------------------------------------------
ka_ = -1*1e+3;
ka__=  1*1e+3;
ka = linspace(ka_,ka__,100).';
% ka = [-1; 0; 1];
kb=ka;
Ea = zeros(numel(ka),1); Eb=Ea;
for i_=1:numel(ka)
  % fwd & obj fnc
  a_ = a - ka(i_)*ga;
  Ea_ = cross2d(a_,b);
  Ea_ = sum(Ea_(:).^2) / numel(Ea_);
  Ea(i_) = Ea_;
  % ---
  % fwd & obj fnc
  b_ = b - kb(i_)*gb;
  Eb_ = cross2d(a,b_);
  Eb_ = sum(Eb_(:).^2) / numel(Eb_);
  Eb(i_) = Eb_;
end
% parabola approx
p = polyfit(ka,Ea,2);
% find zero of parabola (update = -step*gradient)
step_a = -p(2)/(2*p(1));
% parabola approx
p = polyfit(kb,Eb,2);
% find zero of parabola (update = -step*gradient)
step_b = -p(2)/(2*p(1));
fprintf('  step size for a = %2.2d\n',step_a);
fprintf('  step size for b = %2.2d\n',step_b);
% ------------------------------------------------------------------------------
figure;
hold on
plot(ka,Ea,'k.-','markersize',30)
plot(kb,Eb,'r.-','markersize',30)
hold off
legend({'a','b'})
xlabel('Gradient perturbation value')
ylabel('Norm of xgradient')
title('Line search')
simple_figure()
%%{
% ------------------------------------------------------------------------------
%
%         xgradients
%
% ------------------------------------------------------------------------------
a_=a; b_=b;
ka = linspace(-1,2,10).';
kb = linspace(-1,2,10).';
ni = 500;
tol_=1e-20;
da_ = zeros(size(a));
db_ = zeros(size(b));
% ------------------------------------------------------------------------------
prompt = '\n\n    do you want a and b to be like each other (ab)\n    or a to be like b (a)\n    or b to be like a (b)\n\n    choice: \n';
a_b   = input(prompt,'s');
% ------------------------------------------------------------------------------
tic;
if strcmp(a_b,'ab')
  [a_,b_,da_,db_] = cross_ab(a,b,ni,tol_,ka,kb,Dx,Dz);
elseif strcmp(a_b,'a')
  [a_,da_] = cross_a(a,b,ni,tol_,ka,Dx,Dz);
elseif strcmp(a_b,'b')
  [b_,db_] = cross_b(a,b,ni,tol_,kb,Dx,Dz);
end
toc;
% ------------------------------------------------------------------------------
% new values 

figure;
subplot(2,1,1)
hold on
fancy_imagesc(a_,x,z)
c=min([min(a(:)) min(a_(:))]);
c_=max([max(a(:)) max(a_(:))]);
caxis([c c_]);
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
title('new a');
simple_figure();

subplot(2,1,2)
hold on
fancy_imagesc(b_,x,z)
c=min([min(b(:)) min(b_(:))]);
c_=max([max(b(:)) max(b_(:))]);
caxis([c c_]);
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
simple_figure();
title('new b');
% ------------------------------------------------------------------------------
% true values

figure;
subplot(2,1,1)
hold on
fancy_imagesc(a,x,z)
c=min([min(a(:)) min(a_(:))]);
c_=max([max(a(:)) max(a_(:))]);
caxis([c c_]);
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
set(gca,'YTick',[])
set(gca,'XTick',[])
title('a')
simple_figure()

subplot(2,1,2)
hold on
fancy_imagesc(b,x,z)
c=min([min(b(:)) min(b_(:))]);
c_=max([max(b(:)) max(b_(:))]);
caxis([c c_]);
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
title('b')
simple_figure()
% ------------------------------------------------------------------------------
% updates

figure;
subplot(2,1,1)
hold on
fancy_imagesc(da_,x,z)
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
title('stack of a-updates');
simple_figure();

subplot(2,1,2)
hold on
fancy_imagesc(db_,x,z)
plot(circ_ax,circ_az,'c--','linewidth',2)
plot(xbox,zbox,'c--','linewidth',2)
hold off
axis ij
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
title('stack of b-updates');
simple_figure();
% ------------------------------------------------------------------------------
% % borehole
% figure;
% hold on;
% plot(a_(:,binning(x,10)),z,'k');
% plot(b_(:,binning(x,10)),z,'b');
% plot(a(:,binning(x,10)),z,'k--');
% plot(b(:,binning(x,10)),z,'b--');
% xlim([-0.1,1.1])
% axis ij
% xlabel('Value')
% ylabel('Depth')
% title('New black')
% simple_figure();
% ------------------------------------------------------------------------------
% fig = gcf;fig.PaperPositionMode = 'auto'; print(gcf,'xgrad-vu','-dpng','-r600')
% ------------------------------------------------------------------------------
%%}