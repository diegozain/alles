function geom_median(r)
# diego domenzain
# spring 2018 @ TUDelft
# ------------------------------------------------------------------------------
# calculate geometric median from set of points r.
# r has as many points as rows.
#
# uses Weiszfeld's algorithm,
# which is just a potential energy minimizer.
# ------------------------------------------------------------------------------
# using LinearAlgebra (for the function "norm()")
# ------------------------------------------------------------------------------
nr = size(r,1);
ndim = size(r,2);
normi=0.0;
numerator = zeros(ndim,);
# y ⟵ initial guess
y = (sum(r)/nr)*ones(1,ndim);
if nr==1
  y = r;
  return;
end
iter = 0;
c_ = 0;
c__ = Inf;
while (c__ > 1e-12 && iter < 1e+4)
  c = 0;
  numerator[:].=0.0;
  denominator = 0;
  for ii = 1:nr
    normi = norm(r[ii,:].-y);
    numerator = numerator .+ (r[ii,:]./(normi+1e-9));
    denominator = denominator + (1/(normi+1e-9));
    c = c + normi;
  end
  y = numerator / denominator;
  c__ = abs(c-c_);
  c_ = c;
  iter = iter + 1;
end
# ------------------------------------------------------------------------------
return y
end
