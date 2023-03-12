function [rgb,gfofx] = normalizergb(rgb,dmin,dmid,dmax)
% diego domenzain
% •◧⬣·
% ------------------------------------------------------------------------------
%         f                g
% [0 , 1] ⟶ [dmin , dmax] ⟶ [0 , 1]
%
%                 h = g ∘ f
%
% ------------------------------------------------------------------------------
%                                    🔵🟢
%
% dmin=0.6; dmid=1; dmax=1.9;
% rgb = cytwombly_();
% rgb = normalizergb(rgb,dmin,dmid,dmax);
% data=rand(10)*(dmax-dmin) + dmin;
% figure;imagesc(data);colormap(rgb);colorbar
%
% test with testcolor.m
% ------------------------------------------------------------------------------
nrgb = size(rgb,1);
% ------------------------------------------------------------------------------
%
%
%
% ----------------------------------------------------------------------------
h = [dmin ; dmid ; dmax];
x = [0 ; 0.5 ; 1];
y = linspace(0,1,100);
f = interp1(x,h,y);
% ----------------------------------------------------------------------------
g = linspace(dmin,dmax,nrgb).';
gfofx = interp1(f,y,g);
% ----------------------------------------------------------------------------
z = linspace(0,1,nrgb).';
rgb = interp1(z,rgb,gfofx);
% ------------------------------------------------------------------------------
%
% 🐛
%
% ------------------------------------------------------------------------------
% figure;
% hold on;
% plot(g,gfofx,'linewidth',5);
% for ii=1:numel(f)
% plot(f(ii),y(ii),'.','markersize',40)
% end
% hold off;
% axis tight;
% axis square;
% set(gcf,'units','normalized','outerposition',[0 0 0.3 0.5]);
end
