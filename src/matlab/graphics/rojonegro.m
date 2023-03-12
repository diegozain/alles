function rgb = rojonegro(ncolo)
% diego domenzain
% •◧⬣·
% ------------------------------------------------------------------------------
% dmin=0.6; dmid=1; dmax=1.9;
% rgb = rojonegro();
% rgb = normalizergb(rgb,dmin,dmid,dmax);
% ------------------------------------------------------------------------------
% coloritos  = [0 0 0; 0.875 0.183 0.121; 0.8828 0.7773 0.0859];
coloritos  = [0 0 0; 0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250];
% ------------------------------------------------------------------------------
x = linspace(0,1,size(coloritos,1));
if nargin<1
  ncolo = 10000;
end
rgb = interp1(x,coloritos,linspace(0,1,ncolo));
end
