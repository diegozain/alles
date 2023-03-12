function rgb=crazy_mellow(ncolo)
% diego domenzain
% •◧⬣·
% ------------------------------------------------------------------------------
% dmin=0.6; dmid=1; dmax=1.9;
% rgb = crazy_mellow();
% rgb = normalizergb(rgb,dmin,dmid,dmax);
% ------------------------------------------------------------------------------
coloritos=[0.0353,0.8824,0.9569; 0.4824,0.1961,0.6588; 0.865,0.865,0.865; 0.9882,0.6,0.0196; 0.0353,0.9569,0.2353];
% ------------------------------------------------------------------------------
x = linspace(0,1,size(coloritos,1));
if nargin<1
  ncolo = 10000;
end
rgb = interp1(x,coloritos,linspace(0,1,ncolo));
end
