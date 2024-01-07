function mima = MiMa(IM)
% Calculate min and max values of image

mi = double(min(IM(:)));
ma = double(max(IM(:)));
mima = [mi,ma];

