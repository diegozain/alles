function [ia,ib,ic,id] = get_iabcd(iabcd,na,nb,nc,nd)
% diego domenzain
% august 2021
% ------------------------------------------------------------------------------
% mat4d = zeros(na,nb,nc,nd);
% mat4d_= zeros(na*nb*nc*nd,5);
%
% mat4d(ia,ib,ic,id) = iabcd;
% mat4d_(iabcd,:)  = [ia,ib,ic,id,iabcd];
%
%         _____________
%        /             /|
%       /             / |
%      /_____________/  |
%     |  1           |  |
%     |  2           |  |
% nd  |  3           |  /
%     |  etc         | / nb
%     |______________|/
%           nc
%              .
%              .                  na
%              .
%              .
%         _____________
%        /             /|
%       /             / |
%      /_____________/  |
%     |  1           |  |
%     |  2           |  |
% nd  |  3           |  /
%     |  etc         | / nb
%     |______________|/
%           nc
%
% get nd × nc matrix at very first level,
% squeeze(mat4d(1,1,:,:)).'
%
% for cartesian:
% x y z t
% c b d a
% ------------------------------------------------------------------------------
% get d coordinate
id = mod(iabcd,nd);
if (id==0)
  id=nd;
end
% iabcd = (ia-1)*nb*nc*nd + (ib-1)*nc*nd + (ic-1)*nd + id  ... (*)
icd = mod(iabcd,nc*nd);
if (icd==0)
  icd=nc*nd;
end
% icd = (ic-1)*nd + id from (*)
% get c coordinate
ic = ((icd-id)/nd) + 1;
% ibcd = (ib-1)*nc*nd + icd
ibcd = mod(iabcd,nb*nc*nd);
if (ibcd==0)
  ibcd=nb*nc*nd;
end
% get b
ib = ((ibcd - icd) / (nc*nd)) + 1;
% get a
ia = ((iabcd-ibcd)/(nb*nc*nd)) + 1;
end
