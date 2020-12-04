function r = covid_neigh_(r,d)

nr= size(r,1);
updates = zeros(nr,2);

iter_=0;
niter= 1e+1;

% R = zeros(nr,2,niter);

while iter_ < niter
  r_=r;
  for ir=1:nr
    g = grad_Ei(r,d,ir);
    g = g.';
    H = Hess_Ei(r,d,ir);
    % H = H + 1e-3*eye(2);
    update_ = -H\g;
    update_ = update_.';

    Ei_ = Ei(r,d,ir);
    k=[1e-4; 1e-3; 1e-2; 1e-1];
    step_ = step_Ei(r,-update_,d,Ei_,k,ir);
    if step_<0
      step_=0;
    end
    r(ir,:) = r(ir,:) + step_*update_;
    % updates(ir,:) = step_*update_;
  end
  % r = r + updates;
  iter_=iter_+1;

  if or(sum(isnan(r(:)))>0 , sum(isinf(r(:)))>0)
    r=r_;
    return;
  end

  if mod(iter_,10*(log10(niter)+1))==0
    fprintf('\n just finished iteration # %i',iter_);
  end
  % R(:,:,iter_+1) = r;
end
fprintf('\n\nbye bye\n\n')
end