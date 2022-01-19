function [kpdf,x] = kernel_density(datao,bwth,x);
% diego domenzain
% ene 2022

if (nargin<2)
  bwth = 0.5*std(datao);
end
if (nargin<3)
  datao  = sort(datao);
  xextra = max(abs(diff(datao)));
  xmin = min(datao) - xextra;
  xmax = max(datao) + xextra;
  nx= 1e3;
  x = linspace(xmin,xmax,nx);
  x = x.';
end

nx = numel(x);
ndatao = numel(datao);

kpdf = zeros(nx,1);
for idatao=1:ndatao
  kpdf = (1/(sqrt(2*pi*bwth^2))) * exp( -( (x-datao(idatao)).^2/(2*bwth^2) ) ) + kpdf;
end
kpdf = (1/(ndatao*bwth)) * kpdf;

end
