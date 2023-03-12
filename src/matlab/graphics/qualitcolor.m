function rgb = qualitcolor(ncolo)
% diego domenzain
% •◧⬣·
% ------------------------------------------------------------------------------
% dmin=0.6; dmid=1; dmax=1.9;
% rgb = qualitcolor();
% rgb = normalizergb(rgb,dmin,dmid,dmax);
% ------------------------------------------------------------------------------
verde   = [0.4660 0.6740 0.1880];
azul    = [0.3010 0.7450 0.9330];
azul_   = [0 0.4470 0.7410];
guinda  = [0.6350 0.0780 0.1840];
purpura = [0.4940 0.1840 0.5560];
yellow  = [0.9290 0.6940 0.1250];
naranja = [0.8500 0.3250 0.0980];
amarillo= [0.9290 0.6940 0.1250];
% ------------------------------------------------------------------------------
coloritos  = [azul_; naranja; amarillo; purpura; verde; azul; guinda];
% ------------------------------------------------------------------------------
x = linspace(0,1,size(coloritos,1));
if nargin<1
  ncolo = 10000;
end
rgb = interp1(x,coloritos,linspace(0,1,ncolo));
end
