function TomoSaveMRC(IMS,fname,mode,varargin)
%  Save image series to MRC file for visualization in imageJ
% -----------------------------------------------------------------------
%
%  Permute first two axes (PC) or rotate image 90 degrees (Mac) in order
%  to show the images in ImageJ with the same orientation as in MATLAB
%
%  SYNTAX   TomoSaveMRC(IMS,fname,mode)
%           TomoSaveMRC(IMS,fname,mode,eq)
%           TomoSaveMRC(IMS,fname,mode,eq,cut)
%
%  INPUT    IMS     image series
%           fname   string with filename (including .mrc)
%           mode        data class
%                       0    8 bit integer signed (-128 -> 127)
%                       1   16 bit integer signed (-32768 -> 32767)
%                       2   32 bit real signed (single)
%                       3   16 bit integer complex *
%                       4   32 bit real complex *
%                       6   16 bit integer unsigned (0 - 65535) *
%           eq      (o) equilibrate intensities to bring intensities of
%                       all image to the one of the 0-degree image. This
%                       facilitates the inspection with ImageJ for samples
%                       that strongly attenuate at higher tilt angles.
%           cut     (o) [cmin,cmax]: fraction of cumsum max to be cut off
%                       ([0.01,0.99]) in case of strong outliers in image
%
% LATER     Incorporate pixel resolution...
%
%       *   uint16 seems to be the best solution for attenuation data.
%           Allows the largest dynamic range for 16-bit data. Need to keep
%           scale factor to transform back to actual attenuation values.
%
% -----------------------------------------------------------------------

[~,~,ni] = size(IMS);



%% Check data type

% TASK  Make sure that Matlab array has data type that can be handled by
%       MRC reader in ImageJ

% TODO !!!



%% Equilibrate and set dyanmic range if desired

% TASK  In case of TEM  tomograms of thick films the average intensities
%       reduce with increasing tilt angle. Here the dynamic intensity range
%       in each image will be equailibrated to the one near 0 degrees (
%       here assumed to be image number = round(ni/2) );

if nargin>3
    
    fprintf(' -> Equilibrate all histogram ranges to 0-degree image ... ')

    
    % find all min/max values
    mimas = zeros(2,ni);
    
    
    if nargin<5
        
        % just use absolute minima and maxiam
        
        for i=1:ni
            PrintIndex(i,ni);
            mimas(:,i) = MiMa(IMS(:,:,i));
        end
        
    else
        
        % cut off a certain fraction of the cummulative sum of the
        % histogram at the low and high end
        
        coff  = varargin{2};
        coffl = coff(1);
        coffu = coff(2);
        
        parfor i=1:ni
            PrintIndex(i,ni);
            [ival, nval] = TomoHist(IMS(:,:,i));
            nvc = cumsum(nval);
            nvc = nvc/max(nvc);
            ivalc = ival(nvc>coffl&nvc<coffu);
            mimas(:,i) = MiMa(ivalc);
        end
    end
    
    
    % adjust intensities relative to near 0-degree image
    
    i0 = round(ni/2);
    mr = mimas(:,i0);
    for i=1:ni
        IMS(:,:,i) = mr(1) + ...
            (mr(2)-mr(1))/(mimas(2,i)-mimas(1,i))*(IMS(:,:,i)-mimas(1,i));
    end
end



%%  Rotate and save images

% TASK  MRC files are saved with the goal to inspect with Imagej.
%       Imagej loads MRC images rotated. Compensate to have same
%       orientation as in Matlab

fprintf('\n -> Rotate images for proper orientation in Imagej ... ')

% !!! using rot90 on full data volume needs a lot of memory. Try it image
% by image instead

IMS = permute(IMS,[2,1,3]);
IMS = flip(IMS,2);

% parfor i=1:ni
%     PrintIndex(i,ni);
%     IMS(:,:,i) = rot90(IMS(:,:,i),-1);
% end



%% Convert dynamic range

% TASK  WriteMRC saves 8-bit images as int8 (dynamic range -128... 127)
%       Generate image series of proper class

if isa(IMS,'uint8')
    
    fprintf('\n -> Change uint8 to int8 for WriteMRC ... ')
    
    IMSM = zeros(size(IMS),'int8');
    
    parfor i=1:ni
        IT = double(IMS(:,:,i));
        IMSM(:,:,i) = IT-128;
    end
    
end



%% save images

% ImageJ loads data with intensity shift. compensate
% IMS = IMS - 32768;

fprintf('\n -> Save image series ...')

if exist('IMSM','var')
    WriteMRCVarMode(IMSM,1,mode,fname);
else
    WriteMRCVarMode(IMS,1,mode,fname);
end


% Done
fprintf('\n')

end

