function [nb,nt_,nt__] = hd_nbnt_(nt,fo,dt,nb)
% ------------------------------------------------------------------------------
% diego domenzain. jul 2021
% ------------------------------------------------------------------------------
% nt : total number of time samples
% fo : central frequncy to look for         ( nt__ = ceil(1/fo/dt) )
% dt : dt in time
% nb : desired # of time blocks
%
% example for nt=4000, fo=50, dt=2.5e-4, nb=4
%        [nb,nt_,nt__] = hd_nbnt_(4000,50,2.5e-4,4)
% nb=4, nt_=180, nt__=80.
% ------------------------------------------------------------------------------
% dt_ : desired time duration of time-blocks (dt_=(nt*dt)/(desired # of blocks))
% nt_ : number of samples in the over-lapping time
% ------------------------------------------------------------------------------
% nb  : # of time blocks
% nt_ : number of samples in the over-lapping time
% nt__: overlapping number of time samples   (nt__ = ceil(1/fo/dt))
% ------------------------------------------------------------------------------
%
% the forward model is,
%
% uh = cos_blocs * α + sin_blocs * β
%
% the cos_blocs matrix (nt x nb*nh) looks a little like this:
%
%             nh
%             |
%          _____________________________
%         |        |                    |
% nt_ ->  |    *   |_________     0     |
%         |________| <-nt__  |          |
%         |        |    *    |          | * α
%         |        |_________|          |
%         |   0                  etc    |
%         |_____________________________|
%
%
% each block is of size nt_ x nh.
% they all overlap on nt__ samples.
% this big matrix is of size nt x nb*nh.
% there are nb = (nt-nt__)/(nt_-nt__) blocks.
% α is of size nb*nh x 1.
% ------------------------------------------------------------------------------
% --
% nt_  :: size(1) of blocks
% --
% nt__ :: size(1) of overlap
% take dt * nt__ = 1/fo
%
%     --> nt__ = 1/fo/dt
% --
% nb   :: # of blocks
% each block is of size nt_ x nh
% nt = (# of blocks)*(size(1) of blocks) - collisions
% nt = nb*nt_ - nt__*(nb - 1)
%
% -->    nb = (nt-nt__)/(nt_-nt__)
%
% nb should always be an integer!
%
% assuming nt__ is fixed (because nt__ depends on the frequency to be found),
% we can change nt_ (always an integer) until we get an integer for nb.
% we have,
%
% nt_ = ((nt-nt__) / nb) + nt__
%
% so nb must divide nt-nt__.
%
% we can factor nt-nt__ and then look for an nb which satisfies our estimate
% length for dt_ = nt_*dt,
%
% dt_ ≈ (((nt-nt__)/(some combination of factor(nt-nt__) )) + nt__) * dt
%
% nb = the best combination of factor(nt-nt__) that is close to dt_.
%
% whatever that combination may be, it may change the initial value for nt_!
%
% nt_ = ((nt-nt__) / nb) + nt__.
%
% ------------------------------------------------------------------------------
% overlapping number of time samples
nt__ = ceil(1/fo/dt);
% get number of factors of ntnt
ntnt= nt-nt__;
nfact = 0;
for ifact = 1:floor(sqrt(ntnt))
  if (mod(ntnt,ifact)==0)
    if (ntnt/ifact == ifact)
      nfact = nfact+1;
    else
      nfact = nfact+2;
    end
  end
end
% ------------------------------------------------------------------------------
% get list of these factors
fact_combi = zeros(nfact,1);
ifact_=1;
for ifact = 1:floor(sqrt(ntnt))
  if (mod(ntnt,ifact)==0)
    if (ntnt/ifact == ifact)
      fact_combi(ifact_) = ifact;
      ifact_=ifact_+1;
    else
      fact_combi(ifact_) = ifact;
      fact_combi(ifact_+1) = ntnt/ifact;
      ifact_=ifact_+2;
    end
  end
end
% ------------------------------------------------------------------------------
% get nb by matching to dt_
% desired time duration of time-blocks (dt_=(nt*dt)/(desired # of blocks))
dt_ = (nt*dt)/nb;
% error bucket
err_= zeros(nfact,1);
for ifact=1:nfact
  % estimate dt_
  % dt_esti_ = (((nt-nt__)/( some combination of factor(nt-nt__) )) + nt__) * dt
  dt_esti_ = (((nt-nt__)/( fact_combi(ifact) )) + nt__) * dt;
  err_(ifact) = abs(dt_esti_ - dt_);
end
[~,ifact] = min(err_);
nb = fact_combi(ifact);
% ------------------------------------------------------------------------------
% number of samples in the over-lapping time
nt_ = ((nt-nt__) / nb) + nt__;
end
