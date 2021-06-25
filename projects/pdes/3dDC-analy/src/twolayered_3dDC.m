function u_2l = twolayered_3dDC(sigm,h_,recs_x,recs_z,source_p,source_n,current_)
% ------------------------------------------------------------------------------
% 
% two layer solution within the first layer when source is within first layer.
% 
% 
% Electrical Methods & Geophysical Prosource_pecting, Keller. 1927
% pages 107-111.
% https://archive.org/details/electricalmethod00kell/page/106/mode/1up
% ------------------------------------------------------------------------------
% the idea is that:
% 
%    layer 0 (air usually)
% 
%    this layer continues upward to infinity
% -------------------------
%    layer 1 (where sources and recs are):
% 
%    s at depth
%    layer that is h_ meters deep
% -------------------------
%    layer 2
% 
%    this layer continues downward to infinity
% -------------------------
% 
% ------------------------------------------------------------------------------
% * for each source (or sink), one has to reflect it along the
%   horizontal boundaries at each interface:
% 
%    layer 0 (air usually)
% 
%    source (sink) reflected
% -------------------------
%    layer 1 (where sources and recs are):
% 
%    source (sink)
%    
% -------------------------
%    layer 2
% 
%    source (sink) reflected
% -------------------------
% 
% * and then each reflected image gets reflected again, and again, and so on.
%   for example, the reflected image on layer 2, gets mirrored to layer 0,
%   and the reflected image on layer 0, gets mirrored to layer 2. and so on.
% 
% * at each mirroring, the reflection coefficient needs to be corrected:
%   each time it crosses a boundary, the coefficient is multiplied by the 
%   coefficient of the boundary that it crossed.
% ------------------------------------------------------------------------------
% the code is written so that it computes positions of all the mirrored 
% reflections of the first layer 2, and first layer 0 mirrored reflect.
% 
% these positions are called im_u(pper) and im_l(ower).
% 
% ------------------------------------------------------------------------------
sigm0 = sigm(1);
sigm1 = sigm(2);
sigm2 = sigm(3);

[nrows,ncols] = size(recs_x);

n_accuracy = 20;
% ------------------------------------------------------------------------------
k12 = (sigm1 - sigm2)/(sigm1 + sigm2);
k10 = (sigm1 - sigm0)/(sigm1 + sigm0);

% --- source ( + )
u_source_p = zeros(nrows,ncols);
% image point
im_u= 0;
im_l= 0;
% "reflection/transmission index"
k_u = 1;
k_l = 1;

h = source_p(2);
for in_=1:n_accuracy
  
  if mod(in_,2) == 1
    % -- odd in_
    % image point
    im_u = - im_u - 2*h;
    im_l = - im_l + 2*(h_-h);
    % "reflection/transmission index"
    k_u  = k_u * k10;
    k_l  = k_l * k12;
  else
    % -- even in_
    % image point
    im_u = - im_u + 2*(h_-h);
    im_l = - im_l - 2*h;
    % "reflection/transmission index"
    k_u  = k_u * k12;
    k_l  = k_l * k10;
  end
  
  % -- source & sink image potentials
  
  % source ( + )
  dist_pt_u = sqrt(( recs_x - source_p(1) ).^2 + ( recs_z - (source_p(2)+im_u) ).^2);
  u_im_source_p_u = k_u ./ dist_pt_u;
  
  dist_pt_l = sqrt(( recs_x - source_p(1) ).^2 + ( recs_z - (source_p(2)+im_l) ).^2);
  u_im_source_p_l = k_l ./ dist_pt_l;
  
  u_source_p = u_source_p + (u_im_source_p_u + u_im_source_p_l);
end
% source ( + ) full potential
dist_pt_ = sqrt(( recs_x - source_p(1) ).^2 + ( recs_z - source_p(2) ).^2);
u_source_p_ = 1 ./ dist_pt_;
u_source_p  = (u_source_p_ + u_source_p);

% --- sink   ( - )
u_source_n = zeros(nrows,ncols);
% image point
im_u= 0;
im_l= 0;
% "reflection/transmission index"
k_u = 1;
k_l = 1;

h = source_n(2);
for in_=1:n_accuracy
  
  if mod(in_,2) == 1
    % -- odd in_
    % image point
    im_u = - im_u - 2*h;
    im_l = - im_l + 2*(h_-h);
    % "reflection/transmission index"
    k_u  = k_u * k10;
    k_l  = k_l * k12;
  else
    % -- even in_
    % image point
    im_u = - im_u + 2*(h_-h);
    im_l = - im_l - 2*h;
    % "reflection/transmission index"
    k_u  = k_u * k12;
    k_l  = k_l * k10;
  end
  
  % -- sink image potential
  % sink   ( - )
  dist_pt_u = sqrt(( recs_x - source_n(1) ).^2 + ( recs_z - (source_n(2)+im_u) ).^2);
  u_im_source_n_u = k_u ./ dist_pt_u;
  
  dist_pt_l = sqrt(( recs_x - source_n(1) ).^2 + ( recs_z - (source_n(2)+im_l) ).^2);
  u_im_source_n_l = k_l ./ dist_pt_l;
  
  u_source_n = u_source_n + (u_im_source_n_u + u_im_source_n_l);
end
% sink   ( - ) full potential
dist_pt_ = sqrt(( recs_x - source_n(1) ).^2 + ( recs_z - source_n(2) ).^2);
u_source_n_ = 1 ./ dist_pt_;
u_source_n  = (u_source_n_ + u_source_n);

% source & sink full potential
u_2l = (current_/(4*pi*sigm1)) * (u_source_p - u_source_n);
end