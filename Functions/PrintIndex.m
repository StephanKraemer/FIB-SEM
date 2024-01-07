function PrintIndex(i,imax)
%  function PrintIndex(i,imax)
%
%  print index in a for loop to know where it currently is
%
%  IN   i     current index
%       imax  biggest index
%       mode  number of indices per line
%

nl = 60;                            % characters per line
wn = ceil(log10(imax))+1;           % width of largest number + space
nn = floor(nl/wn);                  % number of numbers fitting a line

coder.extrinsic('sprintf');

ns = sprintf('%u',wn);
form = ['%',ns,'u'];

fprintf(form,i)

if ~mod(i,nn) || i==imax
    fprintf('\n')
end;



end
