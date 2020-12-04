x=linspace(0,20,100);
z=linspace(0,10,50);
% ----------------------------------------------------
a = gaussi(x,z,[fix(numel(z)*0.5),fix(numel(x)*0.5)],1,2);
b = 1-gaussi(x,z,[fix(numel(z)*0.5),fix(numel(x)*0.5)],0.5,2);
% ----------------------------------------------------
figure;
subplot(2,1,1)
fancy_imagesc(a,x,z)
set(gca,'YTick',[])
set(gca,'XTick',[])
title('a')
simple_figure()
subplot(2,1,2)
fancy_imagesc(b,x,z)
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
title('b')
simple_figure()
% ----------------------------------------------------
ha = curva_gauss(a);
hb = curva_gauss(b);
he = ha-hb;
he = normali(he);
% ----------------------------------------------------
figure;
subplot(2,1,1)
fancy_imagesc(he,x,z)
set(gca,'YTick',[])
set(gca,'XTick',[])
colormap(rainbow([0 0.85 1]))
title('eh')
simple_figure()
% ----------------------------------------------------
% e_=ones(size(a));
e_=b;
for i_=1:10
  he_ = curva_gauss(e_);
  eh = he_ - he;
  J = curva_gauss_J(e_);
  g = J*eh(:);
  g = normali(g);
  g = reshape(g,size(a));
  e_ = e_ - 0.5*g;
end
% ----------------------------------------------------
% figure;
% fancy_imagesc(g,x,z)
% colormap(rainbow([0 0.85 1]))
% title('g')
subplot(2,1,2)
fancy_imagesc(e_,x,z)
set(gca,'YTick',[])
set(gca,'XTick',[])
colormap(rainbow())
title('e?')
simple_figure()
figure;
fancy_imagesc(a-e_,x,z)
colormap(rainbow())
set(gca,'YTick',[])
set(gca,'XTick',[])
% colormap(rainbow([0 0.6 1]))
title('b?')
simple_figure()
% ----------------------------------------------------