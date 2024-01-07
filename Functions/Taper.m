function im2 = Taper(im,w)
% taper(im,w): Taper width w of image im, setting corner pixels to avg

imt = double(im);
[sx, sy] = size(imt);

% Scaling factor for image
% Scaling factor for average is (1-fim)
fim = ones(sx,sy,'double');

for k = 1:w
    fim(k,:) = (0.5 - 0.5*cos(pi/(w-1)*(k-1)))*fim(k,:);
    fim(sx+1-k,:) = (0.5 - 0.5*cos(pi/(w-1)*(k-1)))*fim(sx+1-k,:);
    fim(:,k) = (0.5 - 0.5*cos(pi/(w-1)*(k-1)))*fim(:,k);
    fim(:,sy+1-k) = (0.5 - 0.5*cos(pi/(w-1)*(k-1)))*fim(:,sy+1-k);
end

% Weighted sum of image with average intensity
avg = sum(imt(:))/(sx*sy);
im2 = fim.*imt + (1-fim).*avg;