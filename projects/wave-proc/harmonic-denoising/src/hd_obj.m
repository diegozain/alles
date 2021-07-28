function [obj,error_] = hd_obj(data_o,data_,type_)
% diego domenzain. jul 2021 @ AU
% 'sse'   = sum of squared errors
% 'lnsse' = ln( sum of squared errors )
if strcmp(type_,'sse')
  % sum of squared errors
  error_= data_-data_o;
  obj   = sum(error_.^2);
elseif strcmp(type_,'lnsse')
  % log( sum of squared errors )
  error_= data_-data_o;
  obj   = log(sum(error_.^2));
  % obj   = sum(error_.^2);
  error_= (1/obj) * error_;
end
end
