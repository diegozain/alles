function fourier_rt(d,dt)
# diego domenzain
# spring 2018 @ TUDelft
# ------------------------------------------------------------------------------
# fourier: d(r,t), df -> d(r,f), f
# -------------------------------------------------------
nt,nr = size(d);
# fft can't do this by itself so we help it,
nt_ = 2^ceil(log2(nt));
nt_ = convert(Int64,nt_);
nt_extra = nt_-nt;
d = [d; zeros(nt_extra,nr)];
# t data to f domain
f_d_f = fft(d,1);
df = 1/dt/nt_;
f = (-nt_/2:nt_/2-1)*df;
# # NOTE:
# # ifft would go here, BEFORE fftshift,
# d = ifft(f_d_f);
# # trim padded zeros because fft can't do this by itself
# d = d(1:end-nt_extra, : );
# # d would now be (up to machine precission),
# # identical to original d (same size even).
# -
# make it look like something we can actually read from
d = fftshift(f_d_f);
# get rid of negative part
i_start = ceil(nt_/2)+1;
i_start = convert(Int64,i_start);
d = d[ i_start:nt_-1, : ];
f = f[ i_start:nt_-1 ];
# d is of size ( (nt_/2)-1 by nr )
return d,f,df
end
