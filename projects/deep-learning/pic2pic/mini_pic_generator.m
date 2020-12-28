clear
close all
clc
% ------------------------------------------------------------------------------
% 
% 
% this script produces the training set for the mini project.
%
%
% ------------------------------------------------------------------------------
% generate 4! matrices
a=(1:4).';
a_perms= perms(a);
n_perm = factorial(4);
% ------------------------------------------------------------------------------
% % generate 4^4 matrices
% a_perms=zeros(4^4,4);
% n_perm = 4^4;
% i_=1;
% for i1=1:4
%  for i2=1:4
%   for i3=1:4
%    for i4=1:4
%     a_perms(i_,:) = [i1,i2,i3,i4];
%     i_=i_+1;
%    end
%   end
%  end
% end 
% ------------------------------------------------------------------------------
% set smoothing
% larger number = less smoothing
kx=0.1;
ky=kx;
nx_pad = 36;
ny_pad = 36;
% ------------------------------------------------------------------------------
b =zeros(24,24,n_perm);
b_=zeros(24,24,n_perm);
% ------------------------------------------------------------------------------
% generate all pics and put them in cubes
% ------------------------------------------------------------------------------
% choose permutation
for i_perm = 1:n_perm;
 
 % 1st column
 b(1:12,1:12,i_perm) = a_perms(i_perm,1);
 b(13:24,1:12,i_perm) = a_perms(i_perm,2);
 % 2nd column
 b(1:12,13:24,i_perm) = a_perms(i_perm,3);
 b(13:24,13:24,i_perm) = a_perms(i_perm,4);

 b(:,:,i_perm) = b(:,:,i_perm)/4;

 b_(:,:,i_perm) = image_gaussian_pad(b(:,:,i_perm),kx,ky,'LOW_PASS',nx_pad,ny_pad);
end
% ------------------------------------------------------------------------------
b_mini = b;
b_mini_= b_;
% ------------------------------------------------------------------------------
% save
save('../../../data/pic2pic/b_mini','b_mini')
save('../../../data/pic2pic/b_mini_','b_mini_')
% ------------------------------------------------------------------------------
