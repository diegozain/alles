function abmn = xbore_getall_(Tx_,Rx_)
% diego domenzain
% summer 2021 @ AU
% ------------------------------------------------------------------------------
% Tx_ = indexed electrodes capable of source-sink
% Rx_ = indexed electrodes capable of being receivers
% ------------------------------------------------------------------------------
nRx = numel(Rx_);
nTx = numel(Tx_);
nelectrodes = nRx + nTx;
% ------------------------------------------------------------------------------
%              get all borehole quadruples of the form AB - MN
% ------------------------------------------------------------------------------
nabmn= ( ((nTx)*(nTx-1))/2 ) * ( ((nRx)*(nRx-1))/2 );
abmn = zeros(nabmn,4,'uint32');

iabmn= 1;
ib_  = 2;
% ------------------------------------------------------------------------------
for ia=1:(nTx-1)
  a=Tx_(ia);
  for ib=ib_:nTx
    b=Tx_(ib);
    
    im_  = 1;
    in_  = 2;
    for im=im_:(nRx-1)
      m=Rx_(im);
      for in=in_:nRx
        n=Rx_(in);
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