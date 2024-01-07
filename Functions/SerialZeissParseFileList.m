function dataOut = SerialZeissParseFileList( dataIn )
% extract information from Zeiss file list
% -----------------------------------------------------------------------
%
%   SYNTAX  dataOut = SerialZeissParseFileList( dataIn )
%
%   dataIn  .projectFolder
%           .sliceFolder
%
%   dataOut .nSlicesRecorded
%           .sliceNumber
%           .slicePosition
%           .sliceThickness
%           .sliceThicknessMean
%           .sliceThicknessStd
%
% -----------------------------------------------------------------------


%% input data

projectFolder = dataIn.projectFolder;
sliceFolder   = dataIn.sliceFolder;


%% Extract information

% files in slice directory
dList = dir(fullfile( projectFolder,sliceFolder ));

if isempty( dList )
    error('Files not found. Please check if path is correct.')
end


% exclude folders '.', '..' and 'chips'
isDir = cat( 1,dList.isdir );
dList( isDir==1 ) = [];

% number of slices
nSlicesRecorded = max(size( dList ));


% file names without '.' and '..'
% search for start and end file (assume uniqueness)
sliceName = strings(nSlicesRecorded,1);
for i=1:nSlicesRecorded
    
    sliceName(i) = convertCharsToStrings( dList(i).name );

end

% parse slice numbers and positions
sliceNumber   = zeros(nSlicesRecorded,1);
slicePosition = zeros(nSlicesRecorded,1);

for i=1:nSlicesRecorded
    
    sni = sliceName(i);
    sliceNumber(i) = str2double(extractBetween( sni,"slice_","_" ));
    slicePosition(i) = str2double(extractBetween( sni,"z=","um" ));
    
end


% slice thicknesses
sliceThickness = zeros(nSlicesRecorded,1);

for i=2:nSlicesRecorded
    
    sliceThickness(i) = slicePosition(i) - slicePosition(i-1);
    
end


% mean slice thickness and standard deviation
sliceThicknessMean = mean( sliceThickness );
sliceThicknessStd  = std( sliceThickness );



%% Set output data

dataOut.nSlicesRecorded = nSlicesRecorded;
dataOut.sliceName       = sliceName;
dataOut.sliceNumber     = sliceNumber;
dataOut.slicePosition   = slicePosition;
dataOut.sliceThickness  = sliceThickness;
dataOut.sliceThicknessMean = sliceThicknessMean;
dataOut.sliceThicknessStd  = sliceThicknessStd;



end
