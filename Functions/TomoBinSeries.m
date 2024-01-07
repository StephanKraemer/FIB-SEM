function [ISB, po] = TomoBinSeries(IS,bin,verbose)
%  function ISB = TomoBinSeries(IS,bin)
% -----------------------------------------------------------------------
%  bin image series
%
%  Input
%   IS          image series nx x ny x np
%   bin         bin factor
%   verbose     (true/false) show image index
%
%  Output
%   ISB         binned image series
%   po          origin of ISB within IS *
%   
%  (*)  for images with uneven nx and ny an area will be cut out based on
%       smallest number of bins that fit. The origine of that image will be 
%       floor(mod(nxy,bin)/2)+1
%
% -----------------------------------------------------------------------


%% dummie input variables from workspace

% IS = IA;
% bin = 4;



%% prepare data

[nx, ny, np] = size(IS);

nbx = floor(nx/bin);
nby = floor(ny/bin);
po  = [floor(mod(nx,bin)/2)+1, floor(mod(ny,bin)/2)+1];

ISB = zeros([nbx, nby, np], class(IS));
IB  = zeros(nbx, nby);



%% run through images

% fprintf('\n-> Bin image series\n')

for i=1:np

    if verbose == true
        PrintIndex(i,np);
    end

    IT = double(IS(po(1):po(1)-1+nbx*bin, po(2):po(2)-1+nby*bin, i));
    
    for k=1:bin
        for l=1:bin
            IB = imadd(IB,IT(k:bin:end,l:bin:end));            
        end
    end
    
    ISB(:,:,i) = IB./(bin*bin);
    IB(:,:) = 0;
end



