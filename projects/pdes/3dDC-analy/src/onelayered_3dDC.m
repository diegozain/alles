function u_1l = onelayered_3dDC(sigm,recs_x,recs_z,source_,current_)
% diego domenzain @ aarhus uni, summer 2021
% 
% analytic solution for a halfspace.
% 
% -∇⋅∇φ = s
% ------------------------------------------------------------------------------

dist_  = sqrt( (recs_x - source_(1)).^2 + ( recs_z + source_(2) ).^2 );
dist_im= sqrt( (recs_x - source_(1)).^2 + ( recs_z - source_(2) ).^2 );

u_1l = (current_/(4*pi*sigm)) * (1./dist_ + 1./dist_im);
% ------------------------------------------------------------------------------
end