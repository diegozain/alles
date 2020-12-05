function h = beamPlot( ev, px, py, pmax, w )
% matt haney
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% get the slowness-vectors for plotting
pxArray = px(1,:);
pyArray = py(:,1);
freq = w / 2 / pi;

% prepare for polar plotting instead of cartesian
rp = sqrt( px.^2 + py.^2 ); % compute ray parameter at each node in grid
ev( rp >= pmax ) = 0; % set these to zero so they do not show up

% load('yet_white'); % A colorscale from Elmer Ruigrok

% change over to the slowness grid
h = figure('Color','w');
pcolor(pyArray,pxArray,rot90(real(ev)',2)); shading('interp'); hold on;
% colormap(yet_white); % set better colormap
c = colorbar; ylabel(c,'Beam Power'); % colorbar title
axis xy; axis('square'); 
axis(gca,'off'); % turn off x and y labels
axis( gca, pmax*[-1.15 1.15 -1.15 1.15] ); % set axis limits

% % plot Haney's estimate of velocity
% % apparent velocity to plot in figures that have slowness axes
% velocity = 3.4; % [km/s]
% % The circle in (Kx,Ky) space that corresponds to the apparent velocity
% t  = linspace( 0, 2 * pi );
% xv = ( 1 / velocity ) * cos(t);
% yv = ( 1 / velocity ) * sin(t);
% % plot a circle indicating surface wave slowness
% plot( xv, yv, 'w', 'LineWidth', 5 ); % plot velocity circle

% plot radial spokes
th =  (0:5) * 2 * pi / 12; % we have 6 spokes
cst = cos(th);
snt = sin(th);
cs = [-cst; cst];
sn = [-snt; snt];
plot( pmax*cs, pmax*sn, 'k--', 'LineWidth', 2 );

% annotate spokes in degrees
rt = 1.1 * pmax; % rayparameter max limit for annotations
for ii = 1 : length(th)
    text( rt*cos(th(ii)), rt*sin(th(ii)), int2str(th(ii)*180/pi), 'HorizontalAlignment', 'Center' );
    th2 = pi + th(ii); % shift by 180 to plot on other side as well
    text( rt*cos(th2), rt*sin(th2), int2str(th2*180/pi), 'HorizontalAlignment', 'Center' );
end
text( 1.25*pmax, -0.25*pmax, 'Backazimuth [deg]' );

%----------------------------------------------------------------------
% Ray-parameter values
drp = pmax / 3; % we plot three circles at different rayparameters
rpArray = 0 : drp : pmax;

% plot the circle for ray-parameter
dth = 1; % Azimuth values
thArray = ( 0 : dth : 360 ); % [deg]
for ii = 2 : numel(rpArray)-1
    plot( rpArray(ii).*cosd(thArray), rpArray(ii).*sind(thArray), 'k--', 'LineWidth', 2 );
end
% plot the outer ring solid line
ii = ii + 1;
plot( rpArray(ii).*cosd(thArray), rpArray(ii).*sind(thArray), 'k', 'LineWidth', 4 );

cArray = 1./rpArray; % convert ray parameter to phase velocity for ring annotation
% annotate ray-parameters
cshift = cosd( 45 ); % annotate at 45 degrees
sshift = sind( 45 );
for ii = 2 : numel(rpArray)-1
    text( (rpArray(ii)+drp*0.05)*cshift, (rpArray(ii)+drp*0.05)*sshift, num2str(cArray(ii),'%2.2f'), 'VerticalAlignment', 'Bottom' );
end
% plot the units on the outside ring
ii = ii + 1;
text( (rpArray(ii)+drp*0.05)*cshift, (rpArray(ii)+drp*0.05)*sshift, [num2str(cArray(ii),'%2.2f') ' [km/s]'], 'VerticalAlignment', 'Bottom' );

% set view to clockwise from North
view([90,-90]);

title(['BEAM at ' num2str(freq) ' [Hz]'],'Position',[-1.2*pmax -1.2*pmax -1]);

set( findall( h, '-property', 'FontSize' ), 'FontSize', 14 );
set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Normal' );

end