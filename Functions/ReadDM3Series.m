function [IMS dx dy] = ReadDM3Series(n1,inc,n2,file,part,o,w)
% Read image series saved as DM3 files from a dedicated folder
% __________________________________________________________________
% Summary:
%    To read complete series set n1=1, inc=1, n2=n. To read part
%    series use n1 and n2 for start and end number of image. Set inc
%    to 2,3,.. to skip images.
%
%    For the first image a log file is created that contains all
%    DM3 parameters, inclusive microscope parameters.
%
% Input:
%    n1     number of first image in list  [1]
%    inc    increment (to skip images)  [1]
%    n2     number of last image in list  [n]
%           Use 1 and n for reading all images 
%    file   Array with filenames (n,nc) (string must have same length nc)
%    part   Array of indices of top left and bottom right corner of
%           section to be extracted [tlx,tly,brx,bry]
%    o      'outliers' Remove outliers, else load as is
%    w      width of tile to search for outliers (reduce outliers to
%           5 x standard deviation of image   [32]
%
% Output:
%    IMS   data cube (sx,sy,n) of n images
%    dx    scaling factor pixel -> nm on x-axis
%    dy    scaling factor pixel -> nm on x-axis
% __________________________________________________________________

% Image dimensions
sx = part(3)-part(1)+1;
sy = part(4)-part(2)+1;

n = floor((n2-n1)/inc)+1;
fprintf('%3u images extracted from series\n', n);
IMS = zeros(sx,sy,n);


% ___________________________________________________________________
% Read individual images and remove outliers if desired
% ___________________________________________________________________
for i=n1:inc:n2
    fprintf('\r %u  %9s\r',i,file(i,:));
    if i==n1
        logfile = strrep(file(i,:),'.dm3','.log');
        fprintf('Save DM3 info to file %s', logfile);
        [M sxt syt] = MyReadDM3(file(i,:),logfile);
        dx = sxt;
        dy = syt;
    else
        [M,~,~] = MyReadDM3(file(i,:));
    end
    IM = M(part(1):part(3),part(2):part(4));
    if strcmp(o,'outliers')
       IM = RemoveOutlier(IM,w,5,5,'');
    end
    indx = floor((i-n1)/inc)+1;
    fprintf('Save as image %3u\n',indx);
    IMS(:,:,indx) = IM;   
end
 