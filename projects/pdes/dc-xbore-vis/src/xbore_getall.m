function abmn = xbore_getall(ne,nTx)
% diego domenzain
% summer 2021 @ AU
% ------------------------------------------------------------------------------
% ne = total number of electrodes
% nTx= number of electrodes that can be source-sink pairs
% ------------------------------------------------------------------------------
nRx = ne - nTx;
% ------------------------------------------------------------------------------
%              get all borehole quadruples of the form AB - MN
% 
% assumes both boreholes have the same number of electrodes,
% and that all electrodes from the "source-sink" borehole (Tx) can be AB pairs.
% 
% for the more complicated case when Tx is like in Adapt, look for 
%                                                                xbore_getall_.m
% ------------------------------------------------------------------------------
nabmn= ( ((nTx)*(nTx-1))/2 ) * ( ((nRx)*(nRx-1))/2 );
abmn = zeros(nabmn,4,'uint32');

% nRx_ is where the receivers begin 
nRx_=nTx+1;

iabmn= 1;
ib_  = 2;
% ------------------------------------------------------------------------------
for ia=1:(nTx-1)
  a=ia;
  for ib=ib_:nTx
    b=ib;
    
    im_  = nRx_;
    in_  = nRx_+1;
    for im=im_:(ne-1)
      m=im;
      for in=in_:ne
        n=in;
        if im<in
          abmn(iabmn,:) = [a,b,m,n];
          iabmn=iabmn+1;
        end
      end
    end
    in_=in_+1;
    im_=im_+1;
  end
  ib_=ib_+1;
end
end