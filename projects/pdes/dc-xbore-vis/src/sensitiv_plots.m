function sensitiv_plots(electrodes,s_i_r_d_std,pseud,nrows,ncols,prct)
% diego domenzain
% julio 2021
% Aarhus Uni
% ------------------------------------------------------------------------------
%
% sensitiv_plots.m
%
% plots random examples of some of the abmn sensitivities...
%                                              !! for all mn with a common ab !!
% good for getting an idea of what the survey looks like.
%
% electrodes
% abmn sequence
% pseudo locations
% # of rows and columns of the final plot
% percentage of colormap to plot
% ------------------------------------------------------------------------------
nx = 1e3;
nz = 1e3;

nsources = size(s_i_r_d_std,2);

xmin = min(electrodes(:,1));
xmax = max(electrodes(:,1));

zmin = min(electrodes(:,2));
zmax = max(electrodes(:,2));

x=linspace(xmin-1,xmax+1,nx);
z=linspace(zmin-1,zmax+1,nz);

[X,Z] = meshgrid(x,z);

dx=x(2)-x(1);
dz=z(2)-z(1);

isources=randi(nsources,nrows*ncols,1);
% ------------------------------------------------------------------------------
figure;
for iexample=1:(nrows*ncols)
  % choose abmn configuration
  isource=isources(iexample);
  src = s_i_r_d_std{ isource }{ 1 }(1:2);
  recs= s_i_r_d_std{ isource }{ 2 }(:,1:2);
  psi_  = sensitivity_3dDC(X,Z,dx,dz,electrodes(src(1),:),electrodes(src(2),:),electrodes(recs(:,1),:),electrodes(recs(:,2),:));
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
  % get geometric median of all of them together
  srcs=repmat(src,size(recs,1),1);
  pseud_ = geom_median([electrodes(srcs(:,1),:); electrodes(srcs(:,2),:); electrodes(recs(:,1),:); electrodes(recs(:,2),:)]);
  % ----------------------------------------------------------------------------
  subplot(nrows,ncols,iexample)
  fancy_imagesc(psi_,x,z)
  % rgb=berlin();
  % rgb=crazy_mellow();
  % rgb=crazymellow();
  % rgb=cuatrocolo();
  % rgb=qualitcolor();
  % rgb=rojonegro();
  rgb=cytwombly();

  colormap(rgb);
  caxis(prct*[-max(psi_(:)) max(psi_(:))])
  colorbar('off')
  % hold on;
  % plot(pseud_(1),pseud_(2),'.','markersize',40,'color',[0.1529,0.7686,0.0980])
  % plot(pseud_(1),pseud_(2),'w.','markersize',20)
  % hold off;
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'visible','off')
  simple_figure()
end
simple_figure()
end
