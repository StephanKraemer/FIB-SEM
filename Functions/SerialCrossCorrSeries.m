function [ccOut,IR,I,CC] = SerialCrossCorrSeries(ccIn,IMS)
%  Calculate sample drift from image cross correlation
% -----------------------------------------------------------------------
%  . Calculate image-to-image drift in serial section run using image cross
%    correlation.
%  . A surface fit can be applied to the cross correlation peak providing
%    sub-pixel resolution. The statistical error in determining the peak
%    position is calculated.
%  . If a second output variable is defined, the drift will be corrected.
%    Otherwise only the drift values will be provided.
% 
%
%  SYNTAX   [ccOut] = SerialCrossCorrSeries(ccIn)           *1
%           [ccOut] = SerialCrossCorrSeries(ccIn,IMS)       *2
%           [ccOut,I,IR,CC] = ...                           *3
%   
%           *1 use when input images are files in a folder
%           *2 use when image series in memory as variable IMS
%           *3 for testing: save last pair of images and cross correlation
%
%  
%  INPUT    ccIn    MATLAB structure with parameters for cross correlation
%           .inputType  'files' or 'variable'
%
%           .Input      structure with filenames and folder if 'files' used
%            .projectFolder
%            .sliceName
%            .sliceFolder
%            .iStart
%            .iEnd
%
%           .YNGauss    'yes'   apply Gauss filter
%           .YNgrad     'yes'   calculate gradient image
%           .YNleastCommonArea  'yes' use only common non-0 area for cc *1
%           .YNfit      'yes'   apply sub-pixel surface fit             *2
%
%               *1 in case images were shifted during the run
%               *2 two polynomial fits are applied
%                  fit 1: general fit, plotted over given sub-pixel grid
%                 fit 2: peak postion is parameterized, results of fit 1 are used as
%                 starting values
%
%           .idxRef     index of reference image
%           .idxPlot    image for which fit of cross correlation peak to be
%                       shown (doesn't work yet)
%
%           .Gauss.h    Gaussian kernel size
%           .Gauss.s    Standard deviation of Gaussian
%
%
%           Input parameters for FastCrossCorr (peak fitting done here)
%           .mode       use of the CrossCorrelation
%                       'template'  Template matching. Taper only image
%                           with width w of ~>=2*diameter of template.
%                       'shift'     Shift between two images. Taper both images.
%           .w          taper witdth (taper if w>0)
%
%           .fit        structure with fit parameters
%            .deg           degree of polynomial (2 or 4]
%            .sx2, .sy2     half width of area cropped around
%                           pixel-resolved cc center
%            .dpf           sub-pixel resolution of first fit [0.1,0.01]
%            .plotfitnum        1 or 2
%
%           Image alignment
%           .alignmode  'fl' use floating alignment, otherwise
%           .usefit     1 or 2 (sub-pixel resolution fit type)
%   
%
%  OUTPUT   ccOut       updated MATLAB structure with shift information
%           .dfl1, .dfl2    relative shifts for floating alignment for
%                           fit 1 (fixed sub-pixel grid) and 
%                           fit 2 (parametrized of maximum)
%           .dabs1, .dabs2  corresponding shift relative to reference image
%
%
%  NOTES
%
%  TO DO
%
% -----------------------------------------------------------------------

%% Dummie input


% % number of slices
% [~,~,nz] = size( IMS ); 
% 
% ccIn.YNGauss  = 'yes';
% ccIn.Gauss.h  = 10;
% ccIn.Gauss.s  = 3;
% ccIn.YNgrad   = 'yes';
% ccIn.YNfillBg = 'yes';
% ccIn.Ibg      = 0;
% ccIn.YNleastCommonArea = 'yes';
% ccIn.idxRef   = round( nz/2 );   
% ccIn.idxPlot  = ccIn.idxRef;
% 
% ccIn.mode = 'shift';
% ccIn.w = 0;
% ccIn.YNfit = 'yes';
% ccIn.fit.deg = 4;
% ccIn.fit.sx2 = 3;
% ccIn.fit.sy2 = 3;
% ccIn.fit.dpf = 0.01;
% ccIn.plotfitnum = 2;
% 
% ccIn.usefit = 1;



%% Prepare data

% % test input: Don't need varagarin, Matlab recognizes how many variables
% % are given and does not complain
% 
% fprintf('test\n');
% ccOut = ccIn;
% if nargin>1
%     stest = size(IMS)
% end
% return


if strcmpi( ccIn.inputType,'files' )

    iStart = ccIn.Input.iStart;
    iEnd   = ccIn.Input.iEnd;
    ni     = iEnd - iStart + 1;

    projectFolder  = ccIn.Input.projectFolder;
    sliceFolder    = ccIn.Input.sliceFolder;
    sliceName      = ccIn.Input.sliceName;

else
    if nargin > 1
        % image series dimensions
        [~,~,ni] = size( IMS );
        iStart   = 1;
        iEnd     = ni;
    else
        error('Error . SerialCrossCorrSeries . Please provide image series variable, otherwise change to input type "files" ')
    end
end


% prepare Gaussian filter
if strcmpi( ccIn.YNGauss,'yes' )
    H = fspecial( 'gaussian',ccIn.Gauss.h, ccIn.Gauss.s );
end


% CROSSBEAM
%
% This function considers cases where the XROI has been shifted during a
% run. During the extraction of the XROI the user can choose to either crop
% out the least common area or use the common envelope of all XROIs. In the
% latter case images will contain stripes of 0 intensity. This function
% uses morphological operations to find the least common area of adjacent
% images.
%
% Close examinaiton of the images shows that there are low-intensity
% spill-over pixels that reach beyond the mean edge of the XROI (???). As a
% result the morphologically determined area is larger than general image
% area leading to sharp edges that can affect the precision of the cross
% correlation.
%
% In order to counteract this artifact, the function will taper the edges
% of the found non-0 intensity areas. In order to activate this process,
% please set 
%
%   ccIn.YNleastCommonArea = 'yes'

% prepare taper
if strcmpi( ccIn.YNleastCommonArea,'yes' )

    w     = 65;
    w2    = floor(w/2)+1;
    [X,Y] = meshgrid(1:w,1:w);
    R = sqrt( (X-w2).^2 + (Y-w2).^2 );

    % taper function
    T = cosd( 360/w*R ) + 1;
    T( R > floor(w/2) ) = 0;
    T = T / sum(T(:));

    % structural element for imerode
    s  = strel( 'square',w2+2 );

end

% initialize output data
ccOut.dfl1  = zeros(2,ni);               % floating shifts 1
ccOut.dfl2  = zeros(2,ni);               % floating shifts 2 (after fit)

ccOut.dfl2conf  = zeros(2,2,ni);         % confidence intervals from fit 2



%% Prepare first reference iamge


if strcmpi( ccIn.inputType,'files' )

    % read first slice image
    fname = fullfile( projectFolder,sliceFolder,sliceName(iStart) );
    IR = double(imread( fname ));


else
    
    % select first image from data cube
    IR = double( IMS(:,:,1) );

    % NOTE  existence of IMS was checked above
end


if strcmpi( ccIn.YNpartial,'yes' )

    % crop image
    ixS = ccIn.Crop.ixS;
    ixE = ccIn.Crop.ixE;
    iyS = ccIn.Crop.iyS;
    iyE = ccIn.Crop.iyE;

    IR = IR(iyS:iyE,ixS:ixE);
end

    
if strcmpi( ccIn.YNleastCommonArea,'yes' )

    % reduce region with I>0 by half the filter kernel
    IRM = imerode( uint8( IR>0 ),s );

    % applying filter tapers all edges
    IRMT = imfilter( double(IRM),T,'replicate' );

    % scale image with tapered mask
    IR = IR .* IRMT;

end


if strcmpi( ccIn.YNGauss,'yes' )

    % apply Gaussian filter
    IR  = imfilter( IR,H,'replicate' );

end


if strcmpi( ccIn.YNgrad,'yes' )

    % calculate gradient
    [DX,DY] = gradient( IR );
    IR      = sqrt( DX.^2 + DY.^2 );

end



%% floating cross correlation

% moving index
k = 2;

for i = iStart+1:iEnd
    PrintIndex( i,iEnd );


    if strcmpi( ccIn.inputType,'files' )

        % read current image
        fname = fullfile( projectFolder,sliceFolder,sliceName(i) );
        I = double(imread( fname ));

    else
        I = double( IMS(:,:,k) );
    end


    if strcmpi( ccIn.YNpartial,'yes' )

        ixS = ccIn.Crop.ixS;
        ixE = ccIn.Crop.ixE;
        iyS = ccIn.Crop.iyS;
        iyE = ccIn.Crop.iyE;

        I = I(iyS:iyE,ixS:ixE);

    end

    
    if strcmpi( ccIn.YNleastCommonArea,'yes' )

        % reduce region with I>0 by half the filter kernel
        IM = imerode( uint8( I>0 ),s );

        % applying filter tapers all edges
        IMT = imfilter( double(IM),T,'replicate' );

        % scale image with tapered mask
        I = I .* IMT;

    end


    if strcmpi( ccIn.YNGauss,'yes' )     % Gaussian filter ?

        I = imfilter( I,H,'replicate' );

    end


    if strcmpi( ccIn.YNgrad,'yes' )     % gradient ?

        [DX,DY] = gradient( I );
        I       = sqrt( DX.^2 + DY.^2 );

    end

    % plot fit ?
    if i==ccIn.idxPlot
        plotfit = 'yes';
    else
        plotfit = 'no';
    end

    
    % cross correlate
    ccFl = FastCrossCorr(I,IR,ccIn,plotfit);


    if nargout==4 && i==iEnd

        % export last cross correlation
        CC = ccFl.CC;

    end


    % save shifts
    ccOut.dfl1(:,k) = ccFl.df1;
    ccOut.dfl2(:,k) = ccFl.df2;

    % confidence intervals
    ccOut.dfl2conf(:,:,k) = ccFl.df2conf;

    % make current image new reference
    IR = I;    

    % increase running index
    k  = k+1;

end


% cummulative drift (floating drift measurement starts with the
% first image as reference)

if ccIn.usefit==1
    shft = cumsum(ccOut.dfl1,2);
else
    shft = cumsum(ccOut.dfl2,2);
end

% shift relative to reference image
iRef = ccIn.idxRef - iStart + 1;
shft(1,:) = shft(1,iRef) - shft(1,:);
shft(2,:) = shft(2,iRef) - shft(2,:);

ccOut.dflcum = shft;
ccOut.iref   = ccIn.idxRef;



end


