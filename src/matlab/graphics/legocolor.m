function rgb=legocolor(ncolo)
% diego domenzain
% •◧⬣·
% ------------------------------------------------------------------------------
% dmin=0.6; dmid=1; dmax=1.9;
% rgb = legocolor();
% rgb = normalizergb(rgb,dmin,dmid,dmax);
% ------------------------------------------------------------------------------
coloritos=[0.9412,0.8941,0.2588; 0.0471,0.4824,0.8627; 0.865,0.865,0.865; 0.865,0.865,0.865; 0.8627,0.1961,0.1255; 0,0,0];
% ------------------------------------------------------------------------------
x = linspace(0,1,size(coloritos,1));
if nargin<1
  ncolo = 10000;
end
rgb = interp1(x,coloritos,linspace(0,1,ncolo));
end