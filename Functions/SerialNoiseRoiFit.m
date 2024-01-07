function [NFIT,NRES,gof] = SerialNoiseRoiFit( N, degree )
% Apply polynomial fit to ROI used to estimate noise
% -----------------------------------------------------------------------
%
%   SYNTAX  [NFIT,NRES,gof] = SerialNoiseRoiFit( N, degree )
%
%   INPUT   BG      region of image to be used for noise estimation
%           degree  degree of polynomial (1-8)
%
%   OUTPUT  BGF     fitted intensity
%           RES     residuals, used for noise estimation
%
% -----------------------------------------------------------------------

%% Prepare

% image size
[ny,nx] = size( N );

% 2D coordinate system
[X,Y] = meshgrid( 1:nx,1:ny );

% prepare variables
[xData, yData, zData] = prepareSurfaceData( X, Y, N );

% Set up fittype and options.
fstring = sprintf( 'poly%1u%1u',degree,degree );
ft = fittype( fstring );



%% Fit

% apply fit
[fitresult, gof] = fit( [xData, yData], zData, ft );

% fit surface
NFIT = feval( fitresult,X,Y );

% residual
NRES = N - NFIT;



end