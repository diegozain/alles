function pseudats_plot3d(electrodes_xyz,pseud,data2plot,ptsz,units_name)
% diego domenzain
% abril 2022
% ------------------------------------------------------------------------------
%
% pseulocs_plot3d.m
%
% plot electrodes_xyz & all pseudo locs (no fancy scheme for repetitions)
% needs:
%    electrode positions in xyz
%    pseud
% ------------------------------------------------------------------------------
nelectrodes = size(electrodes_xyz,1);
nabmn = size(pseud,1);
cmiax = [min(data2plot), max(data2plot)];
% ------------------------------------------------------------------------------
xmin = min(electrodes_xyz(:,1));
xmax = max(electrodes_xyz(:,1));

ymin = min(electrodes_xyz(:,2));
ymax = max(electrodes_xyz(:,2));

zmin = min(electrodes_xyz(:,3));
zmax = max(electrodes_xyz(:,3));

xmin_= min(pseud(:,1));
xmax_= max(pseud(:,1));

ymin_= min(pseud(:,2));
ymax_= max(pseud(:,2));

zmin_= min(pseud(:,3));
zmax_= max(pseud(:,3));

xmax = max([xmax , xmax_]);
ymax = max([ymax , ymax_]);
zmax = max([zmax , zmax_]);
xmin = min([xmin , xmin_]);
zmin = min([zmin , zmin_]);
ymin = min([ymin , ymin_]);
% ------------------------------------------------------------------------------
colo = [0.2392,0.5647,0.9294];
colo_= [0.1412,0.1373,0.7569];

colo = [0.0824,0.8392,0.0549];
colo_= [0.1098,0.8392,0.0549];

colo = [0.8392,0.2000,0.0549];
colo_= [0.8392,0.2000,0.0549];
% ------------------------------------------------------------------------------
rgb = qualitcolor();
% ------------------------------------------------------------------------------
figure;
subplot(1,3,1);
hold on;
% -
plot(electrodes_xyz(:,1),electrodes_xyz(:,3),'.','color',[0.3,0.3,0.3],'markersize',30)
% plot all pseudo positions (including repetitions)
scatter(pseud(:,1),pseud(:,3),ptsz,data2plot,'filled')
% -
hold off;
caxis([cmiax(1),cmiax(2)]);
% colormap(rainbow3_cb(1));
colormap(rgb);
hcb = colorbar('southoutside');
hcb.TickLength = 0;
ylabel(hcb,units_name);
axis ij;
axis image;
xlim([xmin-1 xmax+1])
ylim([zmin-1 zmax+1])
xlabel('Length (m)')
ylabel('Depth (m)')
simple_figure()

subplot(1,3,2);
hold on;
% -
plot(electrodes_xyz(:,2),electrodes_xyz(:,3),'.','color',[0.3,0.3,0.3],'markersize',30)
% plot all pseudo positions (including repetitions)
scatter(pseud(:,2),pseud(:,3),ptsz,data2plot,'filled')
% -
hold off;
caxis([cmiax(1),cmiax(2)]);
% colormap(rainbow3_cb(1));
colormap(rgb);
hcb = colorbar('southoutside');
hcb.TickLength = 0;
ylabel(hcb,units_name);
axis ij;
axis image;
xlim([ymin-1 ymax+1])
ylim([zmin-1 zmax+1])
xlabel('Width (m)')
ylabel('Depth (m)')
simple_figure()

subplot(1,3,3);
hold on;
% -
plot(pseud(:,1),pseud(:,2),'.','color',colo_,'markersize',3)
plot(electrodes_xyz(:,1),electrodes_xyz(:,2),'.','color',[0.3,0.3,0.3],'markersize',40)
% -
hold off;
axis ij;
axis image;
xlim([xmin-1 xmax+1])
ylim([ymin-1 ymax+1])
xlabel('Length (m)')
ylabel('Width (m)')
simple_figure()
% ------------------------------------------------------------------------------
% set(gcf,'units','normalized','outerposition',[0 0 0.8 1]);
set(gcf,'units','normalized','outerposition',[0 0 0.8 0.8]);
end
