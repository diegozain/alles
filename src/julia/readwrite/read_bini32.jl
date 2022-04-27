function read_bini32(name_,size_)
# ------------------------------------------------------------------------------
# diego domenzain 2021
#
# 💾:
# save_bin(name_,varname,accu); ⚠️ write this!
#
# 👀:
# var_mat = read_binui32(name_,size_);
# ------------------------------------------------------------------------------
fid=open(string(name_,".bin"));
var_mat=Array{Int32}(undef,size_[1],size_[2])
read!(fid,var_mat);
close(fid);
return var_mat;
end
