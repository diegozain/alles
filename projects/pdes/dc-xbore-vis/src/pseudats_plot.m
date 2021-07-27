function pseudats_plot(electrodes,abmn,pseud,klusters_,data2plot,units_name,cci_path,thics,cmiax)
% ------------------------------------------------------------------------------
% 
% pseudats_plot.m
% 
% plot data in pseudo sections.
% needs:
%    electrode positions in xz
%    pseud
%    klusters_
%    data to plot
%    label for colormap
%    min & max of caxis (optional)
%    path to the nice dense circle packings
%    thickness of circles that plot:
%                                 electrodes                   40
%                                 pseudo locs - not clustered  90
%                                 pseudo locs - clustered      70
%                                 radius of clusters           0.2
% ------------------------------------------------------------------------------
if nargin < 8
  cmiax = [min(data2plot), max(data2plot)];
end

xmin = min(electrodes(:,1));
xmax = max(electrodes(:,1));

zmin = min(electrodes(:,2));
zmax = max(electrodes(:,2));

nklu = size(klusters_,1);
% ------------------------------------------------------------------------------
figure;
hold on;
% plot electrodes
plot(electrodes(:,1),electrodes(:,2),'.','color',[0.5,0.5,0.5],'markersize',thics(1))
% plot all pseudo positions (including repetitions)
scatter(pseud(:,1),pseud(:,2),thics(2),data2plot,'filled')
% plot only repeated positions (clusters)
for iklu=1:nklu 
  nklu_=size(klusters_{iklu},1);
  % first, let's paint white all the repeated locs we already plotted
  scatter(pseud(klusters_{iklu},1),pseud(klusters_{iklu},2),0.4*thics(2),'w','filled')
  % ----------------------------------------------------------------------------
  path_circs   = strcat(cci_path,num2str(nklu_),'.txt');
  dense_circle = load(path_circs);
  % the guy who wrote these put an index up front, so we have to remove that
  dense_circle = dense_circle(:,2:3);
  
  pseud_    = pseud(klusters_{iklu},:) + thics(4)*dense_circle;
  data2plot_= data2plot(klusters_{iklu});
  
  scatter(pseud_(:,1),pseud_(:,2),thics(3),data2plot_,'filled')
  % ----------------------------------------------------------------------------
end
hold off;
colormap(rainbow2_cb(1));
caxis([cmiax(1),cmiax(2)]);
hcb = colorbar;
ylabel(hcb,units_name);
axis ij;
axis image;
xlim([xmin-1 xmax+1]);
ylim([zmin-1 zmax+1]);
xlabel('Length (m)');
ylabel('Depth (m)');
title('Pseudo section');
simple_figure()
end