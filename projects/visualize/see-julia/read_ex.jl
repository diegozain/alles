include("..\\..\\..\\src\\julia\\readwrite\\read_binf64.jl")
include("..\\..\\..\\src\\julia\\readwrite\\read_bini32.jl")
# ------------------------------------------------------------------------------
path_read="..\\..\\..\\..\\dcip\\field\\ejls\\jan22\\e3\\bin\\";
# ------------------------------------------------------------------------------
nx= read_bini32(string(path_read,"x_size"),[3,1]);
x = read_binf64(string(path_read,"x"),[nx[1],1]);

nmn= read_bini32(string(path_read,"mn_size"),[3,1]);
mn = read_bini32(string(path_read,"mn"),[nmn[1],nmn[2]]);
# ------------------------------------------------------------------------------
