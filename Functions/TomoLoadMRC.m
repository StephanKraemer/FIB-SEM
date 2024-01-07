function [IMS,s] = TomoLoadMRC(fmrc,shift)
%  Read MRC file and rotate 90 counter-clockwise
% -----------------------------------------------------------------------
%
%   Load either original MRC image series as recorded via SerialEM or FEI
%   tomo
%
%   SYNTAX  [IMS,s] = TomoLoadMRC(fmrc,shift)
%
%   INPUT   fmrc    text string with file name including '.mrc'
%           shift   32768   (FEI data)
%                   0       (SerialEM)
%           
%           
%   OUTPUT  IMS     Image series
%           s       Matlab structure with header information
%
%   
%   NOTES   
%
% -----------------------------------------------------------------------

fprintf('\n TomoLoadMRC\n')
fprintf(' -> Load data ... ')

[IMS,s] = ReadMRCVarMode(fmrc);

%   Loading original SerialEM data come out rotated in Matlab. 
%   Compensate by rotating counter-clockwise.
%
%   Case    - images recorded with SerialEM on Windows 7, Matlab on Windows
%             10.

fprintf('\n  -> Rotate images ... ')

IMS = rot90(IMS);

IMS = IMS + shift;


% finished
fprintf('\n')

end
