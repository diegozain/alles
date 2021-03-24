function [vx_rx,vz_rx] = lamb_c(rx,sz,Fo,vp,vs,rho,t)
% diego domenzain
% @ CSM, spring 2021
% ------------------------------------------------------------------------------
% taken from lamb.c:
% Paul Michaels BSU code
%
% special case of vp/vs=sqrt(3) 
%
% ------------------------------------------------------------------------------
% rho=1;vp=1;vs=1;rx=100;Fo=1;sz=1;
% t=linspace(0,3,1525).';
% [vx_rx,vz_rx] = lamb_c(rx,sz,Fo,vp,vs,rho,t);
% ------------------------------------------------------------------------------
nt = numel(t);
dt = t(2) - t(1);

vz_rx = zeros(nt,1);
vx_rx = zeros(nt,1);
% ------------------------------------------------------------------------------
delta = sqrt(3);
mu = rho * vs^2;

% ...rayleigh pole root
gama= sqrt(3 + sqrt(3) )/ 2;

% ???
amp = (delta^2 * Fo) ./ (pi * pi * mu * rx);
amp = 1;

small_number = 1e-8;
% ------------------------------------------------------------------------------
for it=1:nt;
  
  tau= (vs./rx).*t(it);
  % --- fixup for  tau=1 or tau=gama ---
  if (tau == 1)
   tau = tau + small_number;
  end
  if (tau == gama)
   tau = tau + small_number;
  end
  
  % ----------------------------------------------------------------------------
  % set vz_rx(it)
  if (tau < 1 / delta)
   vz_rx(it) = 0;
  end
  if (tau > 1 / delta && tau < 1)
    vz_rx(it) = - (pi / 96 ) * (6 - sqrt(3 * sqrt(3) + 5) / sqrt(gama^2 - tau^2) + sqrt(3 * sqrt(3) - 5) / sqrt(tau^2 + sqrt(3)/4 - 3/4) - sqrt(3) / (sqrt(tau^2 - 1/4)));
  end
  if (tau > 1 && tau < gama)
    vz_rx(it) = -(pi / 48) * (6 - sqrt(3 * sqrt(3) + 5) / sqrt(gama^2 - tau^2));
  end
  if (tau > gama)
    vz_rx(it) = - pi / 8;
  end

  % ...apply scale factor
  vz_rx(it) = vz_rx(it) * amp;
	% ----------------------------------------------------------------------------
  % ----------------------------------------------------------------------------
	tau^2 = tau * tau;
	k2 = (delta^2 * tau^2 - 1) / (delta^2 - 1);
	kappa2 = 1 / k2;
	xvir6 = 16 * sqrt(6);

	% ...compute elliptic integrals
	nup = 20 + 12 * sqrt(3);
	num = 20 - 12 * sqrt(3);
	nua8 = 8 * k2;
	nuap = nup * k2;
	nuam = num * k2;
	k1 = sqrt(k2);
	kappa1 = sqrt(kappa2);
	fp = 6 + D4 * sqrt(3);
	fm = 6 - D4 * sqrt(3);
	
	if (tau < (1 / delta))
		vx_rx(it) = 0.;
 	end

	if ((tau > (1 / delta)) && (tau < 1))
		% ellipticF(phi,m)
		% ellipticPi(n,phi,m)
		ka  =ellipticF(pio2,k1);
		pia8=ellipticPi(pio2,k1,nua8);
		piap=ellipticPi(pio2,k1,nuap);
		piam=ellipticPi(pio2,k1,nuam);

		vx_rx(it) = - (tau / xvir6) * (6 * ka - 18 * pia8 + fm * piam + fp * piap);
	end
	if ((tau > 1) && (tau < gama))
		kb  =ellipticF(pio2,kappa1);
		pib8=ellipticPi(pio2,kappa1,8);
		pibp=ellipticPi(pio2,kappa1,nup);
		pibm=ellipticPi(pio2,kappa1,num);

		vx_rx(it) = - ((tau * kappa1) / xvir6) * (6 * kb - 18 * pib8 + fm * pibm + fp * pibp);
	end
	if (tau > gama)
		kb  =ellipticF(pio2,kappa1);
		pib8=ellipticPi(pio2,kappa1,8);
		pibp=ellipticPi(pio2,kappa1,nup);
		pibm=ellipticPi(pio2,kappa1,num);

		vx_rx(it) = - (((tau * kappa1) / xvir6) * (6 * kb - 18 * pib8 + fm * pibm + fp * pibp) + (pi * tau) / (D24 * sqrt(tau^2 - gama2)));
	end

	% ...apply scale factor
	vx_rx(it) = vx_rx(it) * amp;
end
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

end