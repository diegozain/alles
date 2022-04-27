function binning(x,r)
# diego domenzain
# fall 2020 @ CSM
# ------------------------------------------------------------------------------
# bin numbers x onto receivers r:
# 
# i_xr = binning(x,r);
#
# does this by returning indicies on x 
# that map x into r:
#
# binned x on r = x( ix_r ), i.e.
#
# 				x(ix_r) ~ r
#
# IMPOrTANT: assumes x is sorted. r can have repeated elements.
# ------------------------------------------------------------------------------
# make r and x doubles so abs and Inf play nice,
x = convert(Array{Float64},x);
r = convert(Array{Float64},r);
# sort and get indicies for reversing back to original r,
r,i_rss = unsort(r);
# ------------------------------------------------------------------------------
nx = length(x);
nr = length(r);
# ------------------------------------------------------------------------------
ix_r= ones(nr,1);
d_  = maximum(x)+1; 
d__ = maximum(x)+2;
j   = 1;
j_  = j;
j__ = 0;
# ------------------------------------------------------------------------------
for i_=1:nr
	# if i_>nx
		j=1;
	# end
	bol = 1;
	while bol == 1
		d = abs( r[i_]-x[j] );
		# if previous x(j_) was closer, chose that one
		if d_<=d
			d__ = d_;
			j__ = j_;
			bol = 0;					# go to next i -> exit while
			ix_r[i_] = j_; 		# record index for x
			d_ = Inf;					# reset distance
			j_ = j;						# record current position
		# if previous x(j_) is further, move to next one
		elseif d_>d
			d_ = d;
			j_ = j;
			# next while looks at next x(j)
			if j<nx
				j = j+1;
			end
			# this somehow makes repeated r's get repeated ix_r's
			if d__<d_
				if j_-j__< nx
					j_=j__;
				end
				d__ = Inf;
			end
		end
	end
end
# unsort back to original r
ix_r = ix_r[i_rss];
# matlab makes doubles out of everything, but ix_r are indicies,
# julia does this too, smh
ix_r = convert(Array{UInt64},ix_r);
# if r is just a number
nix_r = size(ix_r);
if nix_r[1]==1
	ix_r=ix_r[1];
end
return ix_r
end