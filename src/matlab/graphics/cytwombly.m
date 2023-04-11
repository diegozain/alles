function rgb = cytwombly(ncolo)
% diego domenzain
% •◧⬣·
% ------------------------------------------------------------------------------
% dmin=0.6; dmid=1; dmax=1.9;
% rgb = cytwombly_();
% rgb = normalizergb(rgb,dmin,dmid,dmax);
% ------------------------------------------------------------------------------
coloritos  = [0.0157,0.1294,0.6392; 0.0275,0.7725,0.9569; 1,1,1; 0.9490,0.8941,0.1412; 0.9569,0.4,0.0275];
% ------------------------------------------------------------------------------
x = linspace(0,1,size(coloritos,1));
if nargin<1
  ncolo = 10000;
end
rgb = interp1(x,coloritos,linspace(0,1,ncolo));
end
