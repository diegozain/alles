function [f_d_f,f,df] = fourier_rt_(d,dt)
% diego domenzain
% spring 2018 @ TUDelft
% ------------------------------------------------------------------------------
% fourier: d(r,t), df -> d(r,f), f
% ------------------------------------------------------------------------------
[nt,nr] = size(d);
% fft can't do this by itself so we help it,
nt_=2^nextpow2(nt);
nt_extra = nt_-nt;
d = [d; zeros(nt_extra,nr)];
% nt_=nt;
% t data to f domain
f_d_f = fft(d,[],1);
df = 1/dt/nt_;
f = (-nt_/2:nt_/2-1)*df;
f = fftshift(f);
% % NOTE:
% % ifft would go here, BEFORE fftshift,
% d = ifft(f_d_f);
% % trim padded zeros because fft can't do this by itself
% d = d(1:nt_-nt_extra, : );
% % d would now be (up to machine precission),
% % identical to original d (same size even).
end