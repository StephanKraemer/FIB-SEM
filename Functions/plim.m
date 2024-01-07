function plim(IM,varargin)
%  plim(IM,varargin)
% -----------------------------------------------------------------------
%  plot image with imtool using proper intensity scaling
%
%  (standard)   plim(IM)
%
%  (optional)   plim(IM,title)
%               plim(IM,title,sfact)
%
%  Input
%   IM          image
%   varargin    title   figure title
%               sfact   initial magnification in %
%                       
%
% -----------------------------------------------------------------------

optargin = size(varargin,2);        % parse optional arguments

sfact = 100;
if optargin>1
    sfact = varargin{2};
end


mima = [min(IM(:)),max(IM(:))];
hf   = imtool(IM,mima,'InitialMagnification',sfact);

if optargin,    set(hf,'Name',varargin{1});    end


