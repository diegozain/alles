function [stft_, gau_, fspect, tspect] = stft(x, dt, t, alph, width_)
% ------------------------------------------------------------------------------
% diego domenzain
% circa 2019
% ------------------------------------------------------------------------------
x = x(:);
t = t(:);
nt= numel(x);
% ------------------------------------------------------------------------------
nt_=2^nextpow2(nt)+fix(nt*0.1);
df = 1/dt/nt_;
fspect  = (1:ceil(nt_/2)-1)*df;

nt_extra = nt_-nt;
x = [x; zeros(nt_extra,1)];
tspect = [t ; ((t(nt)+dt) : dt : (dt*nt_extra+t(nt))).'];
% ------------------------------------------------------------------------------
gau_ = zeros(nt,numel(x));
nf   = numel(fspect);
stft_= zeros(nf, nt);
for il = 1:nt

    lo = t(il);

    gau = exp( -alph * ( ( ( tspect - lo ) ./ width_ ).^8 ) );

    x_ = x.*gau;

    x_ = fft(x_);

    x_ = x_(1:ceil(nt_/2)-1);

    stft_(:, il) = x_;
    gau_(il,:) = gau;
end
end
