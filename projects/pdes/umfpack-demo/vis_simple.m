clear
close all
clc
% ------------------------------------------------------------------------------
% 
% solve Ax=b when A is sparse
% 
% ------------------------------------------------------------------------------
nz_=12;
n=5;
I=zeros(nz_,1);
J=zeros(nz_,1);
V=zeros(nz_,1);

b=zeros(n,1);
% ------------------------------------------------------------------------------
% fill A & b
I(1) = 0; J(1) = 0; V(1) = 2;
I(2) = 1; J(2) = 0; V(2) = 3;
I(3) = 0; J(3) = 1; V(3) = 3;
I(4) = 2; J(4) = 1; V(4) = -1;
I(5) = 4; J(5) = 1; V(5) = 4;
I(6) = 1; J(6) = 2; V(6) = 4;
I(7) = 2; J(7) = 2; V(7) = -3;
I(8) = 3; J(8) = 2; V(8) = 1;
I(9) = 4; J(9) = 2; V(9) = 2;
I(10) = 2; J(10) = 3; V(10) = 2;
I(11)= 1; J(11)= 4; V(11)= 6;
I(12)= 4; J(12)= 4; V(12)= 1;

b = [8; 45;-3; 3; 19];
% ------------------------------------------------------------------------------
% build P
P = zeros(n+1,1);
P(1)=0;
ij = 1;
for in=1:n
 p_=1;
 while ( ij<nz_ && J(ij)==J(ij+1) )
  p_ = p_+1;
  ij = ij+1;
 end
 ij=ij+1;
 P(in+1) = P(in) + p_;
end
nz = P(n+1);
% ------------------------------------------------------------------------------
% build and solve A
A = sparse(I+1,J+1,V);
x = A\b;
% ------------------------------------------------------------------------------
% visualize 
figure;
fancy_imagesc(A)
colormap(rainbow2_cb(2))
colorbar('off')
title('A')
simple_figure()

figure;
fancy_imagesc(b)
colormap(rainbow2_cb(2))
colorbar('off')
title('b')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()

figure;
fancy_imagesc(x)
colormap(rainbow2_cb(2))
colorbar('off')
title('x')
set(gca,'xtick',[])
set(gca,'ytick',[])
simple_figure()
% ------------------------------------------------------------------------------
