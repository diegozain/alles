function A = array_transfer(r,theta,vo,wo)
% 
% output
% A = array transfer matrix (receiver x direction angle).
% 
% I like to think of this matrix as a sort of instrument response but for arrays , 
% and in the space−frequency domain.
% Because it is in the space−frequency domain , its nature is static , 
% and operates with the time domain.
% 
% input
% r = receiver locations. (n x 2) matrix.
% vfix = fixed velocity. (1 x 1) matrix.
% theta = discretized space angle. (1 x ntheta) matrix.
% vo = velocity of chosen direction. (1 x 1) matrix.
% thetao = space angle of chosen direction. (1 x 1) matrix.
% wo = fixed frequency angle. (1 x 1) matrix.

nr = numel(r(:,1));

% kstar vector 
%
kx = (wo/vo) * cos(theta); 
ky = (wo/vo) * sin(theta); 
kstar = [kx' ky'];
A = (1/nr) * exp(1i * ( r * kstar' )); 
end