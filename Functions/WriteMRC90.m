function WriteMRC90(IMS,vol,file)
%
%
%  write MRC file using TOM WriteMRC90 with correct image rotation
%  when viewed with imagej

si = size(IMS);
sx = si(1);
sy = si(2);
sn = si(3);

IMST = zeros(si,class(IMS));

fprintf('-> Rotate images 90 degree\n');
for i=1:sn
    PrintIndex(i,sn);
    IMST(:,:,i) = imrotate(IMS(:,:,i),90);
end

WriteMRC(IMST,1,'ZylinderRotate.mrc');
