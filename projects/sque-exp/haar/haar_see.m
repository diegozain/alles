clear
close all
% ------------------------------------------------------------------------------
% 
% visualize the results from the c code haar-test.
% 
% the only fancy thing this script does is build the matrices for the 
% Haar basis. Builds both the blocky and the weird one.
% 
% please only input basis of size equal to a power of 2.
% 
% ------------------------------------------------------------------------------
dt= 1e-2;
t = (0:dt:1).';
nt= numel(t);
% ------------------------------------------------------------------------------
jj = 2;
ii = 0;
[phi_ij,phi_ij_] = haar_wvlets(t,ii,jj);

figure;
hold on
plot(t,phi_ij,'linewidth',5)
plot(t,phi_ij_,'linewidth',3)
hold off
xlabel('Time')
ylabel('Amplitude')
title('Haar wavelet pairs example')
simple_figure()
% ------------------------------------------------------------------------------
% coefficients of "block" basis
a = [1.0000 ; 0.5000 ; 0.3333 ; 0.2500 ; 0.2000 ; 0.1667 ; 0.1429 ; 0.1250];
b = [0.3397 ; 0.1250 ; 0.1620 ; 0.0208 ; 0.1811 ; 0.0083 ; 0.0175 ; 0.0045];
% ------------------------------------------------------------------------------
% haar basis parameters
% n_basis_= 4; % (index j)
n_basis_= log2(numel(a)); % (index j)

n_basis = 2^n_basis_;
% ------------------------------------------------------------------------------
% "block" basis
haar_basis = zeros(nt,n_basis);

jj = n_basis_;
for ib=0:n_basis-1
  ii=ib;
  [phi_ij,~] = haar_wvlets(t,ii,jj);
  haar_basis(:,ib+1) = phi_ij;
end
% ------------------------------------------------------------------------------
% "weird" basis
haar_basis_ = zeros(nt,n_basis);

[phi_ij,~] = haar_wvlets(t,0,0);
haar_basis_(:,1) = phi_ij;

for jj=0:(n_basis_-1)
  for ii=0:( 2^jj - 1 )
    ib = 1 + ii + 2^(n_basis_ - 1 - jj) + ii*(2^(n_basis_ - jj) - 1);
    % fprintf('jj %i , ii %i , ib %i\n',jj,ii,ib)
    
    [~,phi_ij_] = haar_wvlets(t,ii,jj);
    haar_basis_(:,ib) = phi_ij_;
  end
end

haar_basis_ = haar_basis_*sqrt(n_basis);
% ------------------------------------------------------------------------------
figure;

subplot(221)
fancy_imagesc(haar_basis,(1:n_basis),t)
colorbar('off')
axis normal
set(gca,'ytick',[])
set(gca,'xtick',[])
ylabel('Time')
xlabel('Base #')
title('Haar block base')
simple_figure()

subplot(222)
fancy_imagesc(haar_basis_,(1:n_basis),t)
colorbar('off')
axis normal
set(gca,'ytick',[])
set(gca,'xtick',[])
ylabel('Time')
xlabel('Base #')
title('Haar weird base')
simple_figure()
% ------------------------------------------------------------------------------
% visualize the functions in both basis. they should be identical.
a = haar_basis * a;
b = haar_basis_* b;

subplot(2,2,[3,4])
hold on
plot(t,a,'linewidth',8,'color','k')
plot(t,b,'linewidth',3,'color','r')
hold off
xlabel('Time')
legend({'block','weird'})
set(gca,'ytick',[])
set(gca,'xtick',[])
simple_figure()
% ------------------------------------------------------------------------------

