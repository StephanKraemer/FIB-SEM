function SaveImage(IM,file,type,varargin)
%  function SaveImage(IM,file,range,type)
% -----------------------------------------------------------------------
%  Save image to JPG, TIF, using specific contrast range
%
% SYNTAX    SaveImage(IM,file,type)
%           SaveImage(IM,file,type,range)
%
% INPUT     IM     image
%           file   filename
%           type   'JPG', 'TIF'
%           range  (o) [low high] contrast range
%
% -----------------------------------------------------------------------

if nargin>3
    % predefined dyanmic range
    range = varargin{1};
else
    % measure dynamic range
    range = MiMa(IM);
end

% bring intensity to dynamic range [0,1]
mi = range(1);
ma = range(2);

IMT = double(IM);
IM2 = (IMT-mi)/(ma-mi);

% write file
imwrite(IM2,file,type);

end


