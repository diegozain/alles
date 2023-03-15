function plot_abmn(abmn)
% diego domenzain
% fall 2021
% ------------------------------------------------------------------------------
% % try sorting before plotting. looks sweet 🍫
% % 🔢 sort again by common ab
% [ab,mn,ab_mn,iabmn_sort] = abmn_bundle(abmn);
% abmn = abmn(iabmn_sort,:);
% % uo_dc= uo_dc(iabmn_sort);
% % uo_appres = uo_appres(iabmn_sort);
% ------------------------------------------------------------------------------
verde   = [0.4660 0.6740 0.1880];
azul    = [0.3010 0.7450 0.9330];
purpura = [0.4940 0.1840 0.5560];
naranja = [0.8500 0.3250 0.0980];
rojo    = [0.6350 0.0780 0.1840];
% ------------------------------------------------------------------------------
% figure;
hold on;
plot(abmn(:,1),'.','color',naranja,'markersize',20);
plot(abmn(:,2),'.','color',verde,'markersize',15);
plot(abmn(:,3),'.','color',azul,'markersize',10);
plot(abmn(:,4),'.','color',purpura,'markersize',5);
hold off;
% legend({'a','b','m','n'})
axis tight
xlabel('abmn #')
ylabel('Electrode #')
simple_figure()
end
