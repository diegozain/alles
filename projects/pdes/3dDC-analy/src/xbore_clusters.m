function [klusters_,repeated] = xbore_clusters(abmn,pseud)
% diego domenzain
% summer 2021 @ AU
% ------------------------------------------------------------------------------
nabmn = size(abmn,1);
% ------------------------------------------------------------------------------
repeated=[];
nrepeat=0;

nklu =0;
nklu_=2;
klu  =0;
klu_tmp=0;
klusters=[];

iabmn_=2;
for iabmn=1:(nabmn-1)
 for iiabmn=iabmn_:nabmn
    % this ensures we dont count repetitions twice
    klu = find(repeated==iiabmn);
    klu = numel(klu);
    if (norm(pseud(iabmn,:)-pseud(iiabmn,:)) < 1e-3 && klu < 1)
      % ------------------------------------------------------------------------
      % repeated elements (step 1)
      % ------------------------------------------------------------------------
      % the idea is to fill in, and then remove if they are already inside 
      % the list.
      nrepeat=nrepeat+2;
      
      repeated = [repeated; iabmn; iiabmn];
      nrep_=find(repeated==iabmn);
      nrep_=numel(nrep_);
      
      nrep__=find(repeated==iiabmn);
      nrep__=numel(nrep__);
      
      % prune list
      if (nrep_ > 1)
        nrepeat = nrepeat-1;
        repeated(end-1)=[];
      end
      if (nrep__ > 1)
        nrepeat = nrepeat-1;
        repeated(end)=[];
      end
      % ------------------------------------------------------------------------
      % cluster build up (steps 2 & 3)
      % ------------------------------------------------------------------------
      % use switch of previous iabmn to the next iabmn to count clusters.
      % do a small example and print next statement to see why this works:
      % 
      % fprintf('iabmn & iiabmn %i %i\n',iabmn,iiabmn)
      % 
      % iabmn & iiabmn 4 8
      % iabmn & iiabmn 4 19
      % iabmn & iiabmn 6 11  <-- switch
      % iabmn & iiabmn 6 15
      % iabmn & iiabmn 6 22
      % iabmn & iiabmn 6 26
      % iabmn & iiabmn 6 31
      % iabmn & iiabmn 24 29 <-- switch
      % iabmn & iiabmn 24 34
      
      if (klu_tmp~=iabmn)
        nklu    = nklu+1;
        klu_tmp = iabmn;
        
        klusters=[klusters; nklu_];
        nklu_=2;
      else
        nklu_=nklu_+1;
      end
      % ------------------------------------------------------------------------
    end
  end
  iabmn_=iabmn_+1;
end
% ------------------------------------------------------------------------------
% cluster build up (step 4)
% ------------------------------------------------------------------------------
klusters=[klusters; nklu_];
klusters(1)=[];

% repeated
% klusters

klusters_=cell(nklu,1);
iklu_ =1;
iklu__=1;
for iklu=1:nklu
  klusters_{iklu} = repeated(iklu_:(klusters(iklu__)+iklu_-1));
  iklu_ =iklu_+klusters(iklu__);
  iklu__=iklu__+1;
  % an error here means the double precission when comparing 
  % pseud locations is not enough.
  % in that case, you need to space the electrodes further appart.
end
% to check, this should output all coordinates equal
% pseud(klusters_{iklu},:)
end