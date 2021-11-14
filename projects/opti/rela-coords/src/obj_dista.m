function [objfnc_,err_] = obj_dista(dat_, r, OBJ)
% diego domenzain
% 2021

err_   = dat_ - r;
objfnc_= sum(err_.^2);

if strcmp(OBJ,'rms_log')
  objfnc_= log(sum(err_.^2));
  err_   = (1/objfnc_) * err_;
end

end
