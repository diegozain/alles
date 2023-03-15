function [ab,mn,ab_mn,iabmn_sort] = abmn_bundle(abmn)
% ------------------------------------------------------------------------------
% diego domenzain. jul 2021, Aarhus University
% ------------------------------------------------------------------------------
% given a list of abmn quadrupoles, we want to rearrange them by common ab.
%
% instead of doing a weird cell struct like dc_bundle.m, we will build a new
% vector ab_mn that encodes which receivers belong to which source pair.
%
% while we do that, we will also get a list of all unique ab pairs,
% and also a list of all mn pairs ordered in a way that ab_mn tells us which
% of these mn pairs belong to which ab source pair in ab.
%
%         abmn              ab      mn    ab_mn
%        _______           ____    ____    __
%       |      |          |   |   |   |   | |
%       |      |          |   |   |   |   | |  nab
% nabmn |      |   →      |   |   |   |   | |
%       |      |          -----   |   |   ---
%       |      |            2     |   |    1
%       |      |                  |   |
%        ------                   -----
%          4                        2
%
% ------------------------------------------------------------------------------
% so for example, the pair ab(iab,:) has these many ab_mn(iab) receivers in mn.
%
% if we want to know which receivers,
%
% imn__= sum( ab_mn(1:iab) );
% imn_ = imn__ - ab_mn(iab) + 1;
%
% mn receivers for ab_mn(iab)   =   mn(imn_:imn__,:)
%
% ------------------------------------------------------------------------------
% example:
%
% abmn=[1 2 3 4; 2 3 5 6; 1 2 4 5; 2 3 4 5; 1 2 7 8; 5 4 2 1; 5 4 3 1; 4 5 1 2; 4 5 3 4; 4 5 6 7]
%
% [ab,mn,ab_mn,iabmn_sort] = abmn_bundle(abmn);
%
% iab=1;
% % source pair 1 2
% imn__= sum( ab_mn(1:iab) );
% imn_ = imn__ - ab_mn(iab) + 1;
%
% mn(imn_:imn__,:)
%                 3     4
%                 4     5
%                 7     8
% ------------------------------------------------------------------------------
% 👌 see also ab2mn.m 👍
% [ab_,mn_,iabmn_] = ab2mn(iab,ab,mn,ab_mn);
%
% NOTE: the iabmn_ work with pseud too after sorting pseud with iabmn_sort!
% ------------------------------------------------------------------------------
% in order for this entire scheme to work, we need the abmn list
% to be sorted by rows.
% ------------------------------------------------------------------------------
nabmn = size(abmn,1);
[abmn,iabmn_sort] = sortrows(abmn);

ab = abmn(:,1:2);
% ------------------------------------------------------------------------------
% first, we need the size of ab_mn.
% the idea is to go through ab and count each time ab changes.
% each time ab changes we get the # of mn per ab pair.
% ------------------------------------------------------------------------------
nmn_ab = 0;
nmn_ab_= 0;
nab    = 1;
iabmn  = 1;

% we need to shrink the list ab each time we pass an ab pair.
while (iabmn < (nabmn + 1))
  iabmn = nmn_ab_ + 1;
  nmn_ab = 0;
  for iabmn_ = iabmn:nabmn
    % is this ab pair still in our current list ab(iabmn_:nabmn,:) ?
    [iexist,~] = ismember(ab(iabmn,:),ab(iabmn_:nabmn,:),'rows');
    if (iexist==1)
      % no change means the ab pair still has receivers,
      % so we count them (# of mn per ab pair)
      nmn_ab = nmn_ab + 1;
    else
      % a change means ab switched to the next ab pair,
      % so we count it.
      nab = nab + 1;
      % a change also means we're done with this list,
      % so let's go to the new one.
      break;
    end
  end
  nmn_ab_ = nmn_ab_ + nmn_ab;
end
% ------------------------------------------------------------------------------
% now we need to fill in ab_mn and get our new ab and mn.
%
% ------------------------------------------------------------------------------
ab_mn  = zeros(nab,1,'uint32');
ab_    = zeros(nab,2,'uint32');

nmn_ab = 0;
nmn_ab_= 0;
nab    = 1;
iabmn  = 1;
iab    = 1;

while (iabmn < (nabmn + 1))
  iabmn = nmn_ab_ + 1;
  nmn_ab = 0;
  for iabmn_ = iabmn:nabmn
    [iexist,~] = ismember(ab(iabmn,:),ab(iabmn_:nabmn,:),'rows');
    if (iexist==1)
      nmn_ab = nmn_ab + 1;
      ab_mn(iab) = nmn_ab;
      ab_(iab,:) = ab(iabmn,:);
    else
      nab = nab + 1;
      break;
    end
  end
  iab = iab + 1;
  nmn_ab_ = nmn_ab_ + nmn_ab;
end
% lets throw away the temp ab_
ab = ab_;
clear ab_;
% get the mn pairs
mn = abmn(:,3:4);
% ------------------------------------------------------------------------------
end
