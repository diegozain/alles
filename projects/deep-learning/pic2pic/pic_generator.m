clear
close all
clc
% ------------------------------------------------------------------------------
%
%
%         this script is just an example of what should happen.
%
%
% ------------------------------------------------------------------------------
a=(1:9).';
a_perms = perms(a);
% ------------------------------------------------------------------------------
% set smoothing
% larger number = less smoothing
kx=0.1;
ky=kx;
nx_pad = 36;
ny_pad = 36;
% ------------------------------------------------------------------------------
% total permutations
n_perm = factorial(9);
b =zeros(24,24,n_perm);
b_=zeros(24,24,n_perm);
% ------------------------------------------------------------------------------
% generate all pics and put them in cubes
% ------------------------------------------------------------------------------
% choose permutation
for i_perm = 1:n_perm;
 
 % 1st column
 b(1:8,1:8,i_perm) = a_perms(i_perm,1);
 b(8:16,1:8,i_perm) = a_perms(i_perm,2);
 b(16:24,1:8,i_perm) = a_perms(i_perm,3);
 % 2nd column
 b(1:8,8:16,i_perm) = a_perms(i_perm,4);
 b(8:16,8:16,i_perm) = a_perms(i_perm,5);
 b(16:24,8:16,i_perm) = a_perms(i_perm,6);
 % 3rd column
 b(1:8,16:24,i_perm) = a_perms(i_perm,7);
 b(8:16,16:24,i_perm) = a_perms(i_perm,8);
 b(16:24,16:24,i_perm) = a_perms(i_perm,9);

 b(:,:,i_perm) = b(:,:,i_perm)/9;

 b_(:,:,i_perm) = image_gaussian_pad(b(:,:,i_perm),kx,ky,'LOW_PASS',nx_pad,ny_pad);
end
% ------------------------------------------------------------------------------
% save

% ------------------------------------------------------------------------------
