function pseud = xbore_pseudo(electrodes,abmn)
nabmn = size(abmn,1);
pseud = zeros(nabmn,size(electrodes,2));
% ------------------------------------------------------------------------------
for iabmn=1:nabmn
  pseud(iabmn,:) = geom_median(electrodes(abmn(iabmn,:),:));
end
end
