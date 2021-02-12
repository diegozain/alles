function M = g_M(x,z,r)
% diego domenzain
% fall 2017
% Boise State University
% -- in
% x,z = spacial discretizations.
% r = receivers.
% -- out
% M = measuring operator, matrix of size (nd by nu)
% where nd is number of data points, and nu is total number of
% nodes in our discretization
% ---
% easy version, r=x.
% nx = numel(x);
% nz = numel(z);
% nd = numel(r);
% nu = nx*nz;
% M = [eye(nd) zeros(nd,nu-nd)];
% ir = binning(x,r);
% M = [sparse(ir,ir,ones(nd,1)) , sparse(nd,nu-nd)];
% --
% harder version but better. (this code does this version)
% r has to be an (nd by 2) matrix:
% first column are positive receivers, second are negative receivers.
% entries have to be in index coordinate, example:
% positive receiver is (x,z)=(1,0)
% negative receiver is (x,z)=(1,0.5)
% positive receiver in index of x is ix_pos, so x(ix_pos)=1.
% analogous for iz_pos,ix_neg,iz_neg.
% then ir_pos = sub2ind([nx,nz],ix_pos,iz_pos)
% and r=[ir_pos ir_neg].
% ------------------------------------------------------------------------------
nx = numel(x);
nz = numel(z);
nd = size(r,1);
r = r.';
r = r(:);
i_ = (1:nd);
i_ = repmat(i_,[2,1]);
i_ = i_(:);
v = [ones(1,nd); -ones(1,nd)];
v = v(:);
M = sparse(i_,r,v,double(nd),double(nx*nz));
end