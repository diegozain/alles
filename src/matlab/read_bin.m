function var_mat = read_bin(name_,size_,accu)
% ------------------------------------------------------------------------------
% diego domenzain 2021
%
% accu = 'uint32'
% accu = 'double'
%
% 💾:
% save_bin(name_,varname,accu);
%
% 👀:
% var_mat = read_bin(name_,size_,accu);
% ------------------------------------------------------------------------------
fid=fopen(strcat(name_,'.bin'));
var_mat=fread(fid,size_,accu);
fclose(fid);
end
