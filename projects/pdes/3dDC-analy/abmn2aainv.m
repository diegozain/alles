clc
clear
close all
% ------------------------------------------------------------------------------
addpath('src');
% ------------------------------------------------------------------------------
% total number of electrodes
nelectrodes=32; % 32; 
% declare all electrode positions
electrodes=[1*ones(nelectrodes/2,1) , zeros(nelectrodes/2,1), 9+(linspace(1,nelectrodes*0.6,nelectrodes/2)).'; 4*ones(nelectrodes/2,1) , zeros(nelectrodes/2,1), 9+(linspace(1,nelectrodes*0.6,nelectrodes/2)).'];
% # of electrodes that will be Tx
nTx=nelectrodes/2;
% ------------------------------------------------------------------------------
xmin = min(electrodes(:,1));
zmin = min(electrodes(:,3));

xmax = max(electrodes(:,1));
zmax = max(electrodes(:,3));
% ------------------------------------------------------------------------------
% get all abmn pairs
abmn = xbore_getall(nelectrodes,nTx);
nabmn= size(abmn,1);
% ------------------------------------------------------------------------------
figure;
plot(electrodes(:,1),electrodes(:,3),'k.','markersize',40);
axis ij
xlim([xmin-1, xmax+1])
ylim([0, zmax+1])
xlabel('Length (m)')
ylabel('Depth (m)')
title('Survey electrodes')
simple_figure()
% ------------------------------------------------------------------------------
% the idea is to run the "analytical" solution using the 1d forward model 
% in AarhusInv.
% we compare this "analytical" solution to the 3d forward model living in the 
% Line inversion.
% ------------------------------------------------------------------------------
% print to file the way forward-1ddc AarhusInv likes it:
% xa ya za xm ym zm xn yn zn xb yb zb
abmn_xyz_1d = zeros(nabmn,12);
for iabmn=1:nabmn
 abmn_xyz_1d(iabmn,:) = [electrodes(abmn(iabmn,1),:), electrodes(abmn(iabmn,3),:), electrodes(abmn(iabmn,4),:), electrodes(abmn(iabmn,2),:)];
end
% ------------------------------------------------------------------------------
dcp_meat = zeros(nabmn,14);
dcp_meat(:,1:12) = abmn_xyz_1d;
% data (units?)
dcp_meat(:,13) = 1;
% std (why not = 0 ?)
dcp_meat(:,14) = 2;
% ------------------------------------------------------------------------------
% ok, so the .dcp has a weird header that i am not in the mood of writing.
% so this part is the meat and bones.
% will have to fix this shit later.
fileid = fopen('dcp.txt','w');
for iabmn=1:nabmn
  fprintf(fileid,'%1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f\n',dcp_meat(iabmn,:));
end
fclose(fileid);
% ------------------------------------------------------------------------------
% !cat dcp.txt
% ------------------------------------------------------------------------------
% print to file the way forward-3ddc AarhusInv likes it:
% a b m n
% NOTE: this is already done in the matrix abmn.
% ... however, the .rip file likes the electrodes with the depth coordinate 
% negative. whoever wrote all this code was on some crazy acid trip.
electrodes_rip = zeros(nelectrodes,6);
electrodes_rip(:,1) = (1:nelectrodes).';
electrodes_rip(:,2:4) = electrodes;
electrodes_rip(:,4)   = - electrodes(:,3);
% ------------------------------------------------------------------------------
% headers. idk what these do.
header_1 = [nelectrodes , 1 , 1 , 0 , 1 , 0];
header_2 = [1.00, 3.00, 1.00, -1.00, 1.00, 0.00, 2.00, 4.00];
header_3 = [1.00, 1.00, 0.00, 0.00, 2.00, 4.00, 8.00, 0.00, 0.00, 0.00, 1.00, 2.00, 1.00, 2.00, -1.00, 4.00, -1.00, 4.00, 1.00, 6.00, 1.00, 6.00, 0.00];
header_4 = [1.00, 1.00, 0.00, 0.00, 1871.60, 1911.60, 0.00];
% ------------------------------------------------------------------------------
protocol = zeros(nabmn,16);
% gate id ??
protocol(:,1) = 1.0;
% filter id ??
protocol(:,2) = 0.0;
% # of ab per Tx
protocol(:,3) = 2;
% # of mn per Rx
protocol(:,4) = 2;
% index of a electrode
protocol(:,5) = abmn(:,1);
% input current (A)
protocol(:,6) = 1.0;
% WFID ??
protocol(:,7) = 1.0;
% index of b electrode
protocol(:,8) = abmn(:,2);
% output current (A)
protocol(:,9) = -1.0;
% WFID ??
protocol(:,10) = 1.0;
% index of m electrode
protocol(:,11) = abmn(:,3);
% index of n electrode
protocol(:,12) = abmn(:,4);
% measured voltage (V) AarhusInv wont like it if this is 0.
protocol(:,13) = 1.0;
% measured app. res. (ohm.m)
protocol(:,14) = 1.0;
% std
protocol(:,15) = 0.0;
% flag ??
protocol(:,16) = 0.0;
% ------------------------------------------------------------------------------
fileid = fopen('rip.txt','w');
% first header 
fprintf(fileid,'%%Version number\n1\n%%Domain\nTD\n'); 
fprintf(fileid,'%% #Electrode	 #WaveForms	 #GateSet	 #Filters	 StacSet	 #Lines\n');
fprintf(fileid,' %i	 %i	 %i	 %i	 %i	 %i\n\n',header_1);
% electrode positions
fprintf(fileid,'%% ElecID	 UTMx	 UTMy	 Surface	 Depth	 LineID\n');
for ielectrode=1:nelectrodes
 fprintf(fileid,'%2.4f	 %2.4f	 %2.4f	 %2.4f	 %2.4f	 %2.4f \n',electrodes_rip(ielectrode,:));
end
% some weird shit idk what it does #1
fprintf(fileid,'%%StacSetID	 #Pulses	 Amplitudes	 Times\n');
fprintf(fileid,' %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f\n',header_2);
% some weird shit idk what it does #2
fprintf(fileid,'%% WFID	 StackID	 CurrGateType	 CurrGateS	 CurrGateE	 WaveType	 #Points	 t1	 a1	 t2	 a2	 t3	 a3	 t4	 a4	 t5	 a5	 t6	 a6	 t7	 a7	 t8	 a8\n');
fprintf(fileid,' %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f\n',header_3);
% some weird shit idk what it does #3
fprintf(fileid,'%% GateID	 DCStacID	 TransType	 DCGateType	 DCgateS	 DCgateE	 #IPGates\n');
fprintf(fileid,' %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f\n',header_4);
% protocol
fprintf(fileid,'%%GateID	 FilterID	 NTx	 NRx	 Tx1	 Curr	 WFID	 Tx2	 Curr	 WFID	 Rx1	 Rx2	 Volt	 rhoa	 ResSTD	 ResFlag\n');
for iabmn=1:nabmn
 fprintf(fileid,' %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f	 %1.2f\n',protocol(iabmn,:));
end
fclose(fileid);
% ------------------------------------------------------------------------------
