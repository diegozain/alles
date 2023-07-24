function save_bin(name_,varname)
# ------------------------------------------------------------------------------
# diego domenzain 2021
#
# 💾:
# save_bin(name_,varname,accu); 
#
# 👀:
# var_mat = read_binf64(name_,size_);
# ------------------------------------------------------------------------------
fid=open(string(name_,".bin"),"w");
write(fid,varname);
close(fid);
end
