function save_bin(name_,varname,accu)
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
fid = fopen(strcat(name_,'.bin'),'w');
fwrite(fid,varname,accu);
fclose(fid);
end
