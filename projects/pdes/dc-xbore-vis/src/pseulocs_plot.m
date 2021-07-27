function pseulocs_plot(electrodes,pseud)
% diego domenzain
% julio 2021
% Aarhus Uni
% ------------------------------------------------------------------------------
% 
% pseulocs_plot.m
%
% plot electrodes & all pseudo locs (no fancy scheme for repetitions)
% needs:
%    electrode positions in xz
%    pseud
% ------------------------------------------------------------------------------
nelectrodes = size(electrodes,1);
nabmn = size(pseud,1);
% ------------------------------------------------------------------------------
xmin = min(electrodes(:,1));
xmax = max(electrodes(:,1));

zmin = min(electrodes(:,2));
zmax = max(electrodes(:,2));
% ------------------------------------------------------------------------------
figure;
hold on;
% -
plot(electrodes(:,1),electrodes(:,2),'.','color',[0.5,0.5,0.5],'markersize',40)
plot(pseud(:,1),pseud(:,2),'.','color',[0.9412,0.2706,0.1686],'markersize',10)
% -
hold off;
axis ij;
axis image;
xlim([xmin-1 xmax+1])
ylim([zmin-1 zmax+1])
xlabel('Length (m)')
ylabel('Depth (m)')
title('ab.mn pseudo locations')
simple_figure()
end