function sensitiv_plot(electrodes,abmn,pseud,nrows,ncols,prct)
% diego domenzain
% julio 2021
% Aarhus Uni
% ------------------------------------------------------------------------------
% 
% sensitiv_plot.m
% 
% plots random examples of some of the abmn sensitivities
% good for getting an idea of what the survey looks like.
%
% electrodes
% abmn sequence
% pseudo locations
% # of rows and columns of the final plot
% percentage of colormap to plot
% ------------------------------------------------------------------------------
nx = 6e2;
nz = 6e2;

nabmn = size(abmn,1);

xmin = min(electrodes(:,1));
xmax = max(electrodes(:,1));

zmin = min(electrodes(:,2));
zmax = max(electrodes(:,2));

x=linspace(xmin-1,xmax+1,nx);
z=linspace(zmin-1,zmax+1,nz);

[X,Z] = meshgrid(x,z);

dx=x(2)-x(1);
dz=z(2)-z(1);

iabmn_=randi(nabmn,nrows*ncols,1);
% ------------------------------------------------------------------------------
figure;
for iexample=1:(nrows*ncols)
  % choose abmn configuration
  iabmn=iabmn_(iexample);

  source_p = electrodes(abmn(iabmn,1),:);
  source_n = electrodes(abmn(iabmn,2),:);

  rec_p = electrodes(abmn(iabmn,3),:);
  rec_n = electrodes(abmn(iabmn,4),:);
  % ----------------------------------------------------------------------------
  psi_  = sensitivity_3dDC(X,Z,dx,dz,source_p,source_n,rec_p,rec_n);
  % ----------------------------------------------------------------------------
  % % to smooth within a wavelength lo,
  % % ax=1/lo; az=1/lo;
  % % ax=ax*dx;
  % % az=az*dz;
  % ax=1;
  % az=1;
  % ax=ax*dx;
  % az=az*dz;
  % psi_ = smooth2d(psi_,ax,az);
  % ----------------------------------------------------------------------------
  subplot(nrows,ncols,iexample)
  fancy_imagesc(psi_,x,z)
  colormap(bone)
  caxis(prct*[-max(psi_(:)) max(psi_(:))])
  colorbar('off')
  %
  hold on;
  plot(pseud(iabmn,1),pseud(iabmn,2),'.','markersize',40,'color',[0.1529,0.7686,0.0980])
  plot(pseud(iabmn,1),pseud(iabmn,2),'w.','markersize',20)
  hold off;
  %
  set(gca,'xtick',[])
  set(gca,'ytick',[])
end
simple_figure()
end