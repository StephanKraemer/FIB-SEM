%
% NOTE
%
% This template was extracted from the corresponding Matlab script in the
% project demo. At the top it contains all necessary workflow sections. 
% As processes it includes denoising and y-dependent drift correction.
%
% CAREFULLY ADJUST THE UPDATE SECTIONS TO YOUR GIVEN PROJECT
%


%% SESSION 1

% TASK  Set up info structure and start first session

% NEED  - go to 'Align' folder in project folder
%       - make sure that all data from previous sessions are saved 
%       - cleaning the workspace with clearvars command recommended


if exist( 'info','var')

    % warn if info variable exists
    error( 'Workspace contains info structure. Make sure that all relevant data are saved and delete the variable')

else

    % initialize info structure
    info = struct;

end


% UPDATE -- INFO --------------------------------------------------------

% project identifier
info.Name = "<user name> . Run <date> . <run name>";

% current session number
info.Session = 1;

% description
info.Description = "Align images of ESB series";

% Detector signal used in this session
info.Detector = "ESB";


% UPDATE -- WORKFLOW ----------------------------------------------------

% print image indeces during calculations?
verbose = true;

% automatically save info file at end of each process (true/false)
autoSaveInfo = true;


% -----------------------------------------------------------------------


% date created/modified/saved
timeNow             = datetime( 'now' );
info.Dates.Created  = timeNow;
info.Dates.Modified = timeNow;
info.Dates.Loaded   = 0;
info.Dates.Saved    = 0;



%% RUN DATA

% TASK  Collect information about current run

% UPDATE -- RUN ---------------------------------------------------------

info.Run.Instrument         = "Crossbeam";
info.Run.Software           = "Atlas Tomography";
info.Run.Date               = datetime("2022-10-17");

info.Run.Project            = "< project name >";
info.Run.Name               = "< sample name >";

info.Run.Units.Lenght       = "nm";
info.Run.Units.Energy       = "kV";
info.Run.Units.Current      = "pA";

info.Run.pixelSize          = 5;
info.Run.sliceThickness     = 10;

info.Run.Detectors          = "ESB";   % multiple: string array

info.Run.Ebeam.Voltage      = 1.5;
info.Run.Ebeam.Current      = 2000;
info.Run.Ibeam.Voltage      = 30;
info.Run.Ibeam.Current      = 700;

info.Run.Projection         = "Corrected";


% number of digits of slice number in file string
nDigits                     = 5;

% is there anything particular to this run?
info.Run.Notes              = "";

% -----------------------------------------------------------------------

% slice number format
info.Run.sliceNameFormat.nDigits     = nDigits;
info.Run.sliceNameFormat.numberFormatString = ['%0',num2str(nDigits),'d'];



%%  PROJECT FOLDERS

% TASK  . Set project/storage/data folders for loading/storing data

% UPDATE ----------------------------------------------------------------

% project/storage folders ( work folder optional )
info.Project.projectFolder = '< absolute path to project folder >';
info.Project.storageFolder = '< absolute path to storage folder, if used >'; 

% -----------------------------------------------------------------------


% info file name (global path)
infoFile = sprintf( 'info-Session%d.mat',info.Session );
infoFile = fullfile( info.Project.projectFolder,'Align',infoFile );



%% PROJECT NOTES

% TASK  Note particularities of the run and important results throughout
%       the analysis

% PROC  Comment out previous entry to omit multiple entries. Add new note.


% UPDATE ----------------------------------------------------------------

% TASK  ( 'add note','save to Excel', copy from here! )
noteTask = 'add note';

% process number for note ( if applicable, set 0 for general note )
noteProcess = 3;

% first note (comment out when more notes added)
Note = "Start Analysis";


% -----------------------------------------------------------------------

if  ~isfield( info,'Notes' )

    % generate table on first call
    info.Notes = table;

else
    switch noteTask

        case 'add note'

            if ~isempty( Note )

                % count number of notes
                if isempty( info.Notes )

                    newLine = 1;
                else

                    [nl,~] = size( info.Notes );
                    newLine = nl + 1;
                end

                % add new entry
                newNote         = table;
                newNote.Session = info.Session;
                newNote.Step    = noteProcess;
                newNote.Note    = Note;
                newNote.Date    = datetime( 'now' );

                info.Notes(newLine,:) = newNote;

            end


        case 'save to Excel'

            % filename
            fXls = sprintf( 'Notes_Session%u.xlsx', info.Session );
            fXls = fullfile( info.Project.projectFolder,'Align',fXls );

            % save
            writetable( info.Notes,fXls );

        otherwise

            error( 'Please choose the proper task ' )

    end
end



%% PROCESS LIST

% TASK  Initialize process list or save as Excel spread sheet into the 
%       'Align' folder of the project

% UPDATE ----------------------------------------------------------------

% save to Excel spread sheet ( true/false )
saveToExcel = false;

% ----------------------------------------------------------------------

if ~isfield( info,'ProcessList' )

    % generate process list when a session is started
    info.ProcessList = table;


else
    if saveToExcel == true

        % filename
        fXls = sprintf( 'ProcessList_Session%u.xlsx', info.Session );
        fXls = fullfile( info.Project.projectFolder,'Align',fXls );

        % save to Excel file
        writetable( info.ProcessList,fXls );

    end
end


        
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% INFO . Load

% TASK  reload project processing data to continue the current or copy
%       data for a new session

% NEED  go to Align folder

% UPDATE -- WORKFLOW ----------------------------------------------------

% info structure of current session (true/false)
% (otherwise assumed to be a previous session that will be built on)
currentSession = true;

% session number
session = 1;


% PROCESS CONTROL PARAMETERES
% ( suppression of warning ok: process control parameters can be set in
% section 'SESSION X' when starting a session or here when resuming an
% existing one )

% print image indeces during calculations (true/false)
verbose = true; %#ok<NASGU> 

% automatically save info file at end of each process (true/false)
autoSaveInfo = true; %#ok<NASGU>


% -----------------------------------------------------------------------

% info file name
infoFile = sprintf( 'info-Session%d.mat',session );

% load info structure
ldata = load( infoFile );


if currentSession == true

    % extract info structure
    info = ldata.info;

    % save loading time (to compare with saving time)
    info.Dates.Loaded = datetime( 'now' );

    % set global path for info file
    infoFile = fullfile( info.Project.projectFolder,'Align',infoFile );

else

    % extract previous info structure
    infoPrev = ldata.info;

    % save loading time
    infoPrev.Dates.Loaded = datetime( 'now' );

end


% clean up
clearvars ldata



%% INFO . Save

% TASK  Save workflow parameters

% UPDATE ----------------------------------------------------------------
% -----------------------------------------------------------------------

% mark saving time
info.Dates.Saved = datetime( 'now' );

% save
save( infoFile,'info','-v7.3' )



%% INFO . Edit

% TASK  remove last n process steps in order to redo an analysis within
%       current session in case an error occured, a different order of
%       tasks is desired without generating a new session, etc.

% PROC  Example: 8 tasks have been performed and the last three should be
%       removed
%
%       1. Check in the 'Data' column of the process list which info sub
%          structures were added or modified by the last n tasks.
%
%       2. Delete sub structures or relevant structure component. For
%          example if multiple crop tasks were performed, remove the
%          relevant component
%
%          e.g. if only one crop event should be deleted:
%                   info.Crop(2) = [];
%   
%               if entire structure should be deleted:
%                   info = rmfield( info,'Crop');
%
%       3. Delete row in process list
%          
%           e.g. info.ProcessList(6:8,:) = [];



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% IMAGES . Load from folder

% TASK  Load single image, (partial) image series, skipped series
%       Allow transformation to 8 bit

% NEED  sessions    PROJECT DATA
%                   PREPARE . Parse file list

% UPDATE ----------------------------------------------------------------

% input source
projectFolderRead = info.Project.projectFolder;
sliceName         = info.Project.sliceName;

% slice folder ( copy from Windows explorer )
sliceFolder = 'ESB_ExtCroDenYddBin';

% corresponding process number ( see in info.ProcessList )
loadProcess = 7;



% slice increment for saving 
% (1 if single image or non-skipped series loaded )
nInc = 1;

% Flag . reduce dynamic range ( true/false )
dr8bit = false;

% save to info structure ( true/false ) 
saveToInfo = false;

% type of load folder ( 'Project','Storage' )
loadFolderType = 'Project';


% -----------------------------------------------------------------------

% start and end slice
% ( identical if single image loaded )
iStart = info.ProcessList.iStart( loadProcess );
iEnd   = info.ProcessList.iEnd( loadProcess );

% all images to be loaded
iImages = iStart:nInc:iEnd;

% number of images
nImages = max(size( iImages ));


if nImages == 1
    
    % LOAD SINGLE IMAGE

    % read current slice image
    fname = fullfile( projectFolderRead,sliceFolder,sliceName(iStart) );
    
    if dr8bit == true
        
        IT  = double(imread( fname ));
        s8  = double(intmax( 'uint8' ));
        s16 = double(intmax( 'uint16' ));
        IM  = uint8( IT * s8 / s16 );
        
    else
        IM = imread( fname );
    end
    
else
    
    % LOAD IMAGE STACK
    
    % allocate array based on first image
    fname   = fullfile( projectFolderRead,sliceFolder,sliceName(iStart) );
    [ny,nx] = size(imread( fname ));

    if dr8bit == true
        IMS = zeros(ny,nx,nImages,'uint8');
    else
        IMS = zeros(ny,nx,nImages,'uint16');
    end
    
    % running index
    k = 1;
    
    tic
    for i=iImages
        PrintIndex(i,nImages);
        
        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolder,sliceName(i) );
        
        if dr8bit == true
            
            IT  = double(imread( fname ));
            s8  = double(intmax( 'uint8' ));
            s16 = double(intmax( 'uint16' ));
            
            IMS(:,:,k)  = uint8( IT * s8 / s16 );
            
        else
            IMS(:,:,k) = imread( fname );
        end
        
        % update running index
        k = k + 1;
        
    end
    toc
    % 758 sec
    
end


% update current parameters
if saveToInfo == true

    info.CurrentLoad.projectFolderRead     = projectFolderRead;
    info.CurrentLoad.loadFolderType = loadFolderType;
    info.CurrentLoad.sliceFolder    = sliceFolder;

    info.CurrentLoad.iStart         = iStart;
    info.CurrentLoad.iEnd           = iEnd;
    info.CurrentLoad.nInc           = nInc;
    info.CurrentLoad.iRef           = round( 0.5*(iStart+iEnd) );

    info.CurrentLoad.Flag.dr8bit    = dr8bit;
end



%% IMAGES . Save to folder

% TASK  Save image (series) to folder either as 16-bit/8-bit or as tif/jpg
%       Read from folder or use image series stored as variable

% NEED  if image series is in memory and currently named IMSD (which stores
%       the results of the last process) save it first into variable IMS
%       which will be used for saving


% UPDATE -- WORKFLOW ----------------------------------------------------

% use image dimensions from process step
pStep = 1;


% UPDATE -- IN/OUT ------------------------------------------------------


% use input folder ( true/false )
useFolder = false;

if useFolder == true

    % input folders
    projectFolderRead = info.Project.projectFolder;
    sliceFolderLoad   = 'ESB_ExtYddInv';

end


% output
projectFolderSave = info.Project.projectFolder;
sliceName         = info.Project.sliceName;

sliceFolderSave = 'ESB_ExtDenYddInv_8bit';


% UPDATE -- PROCESS -----------------------------------------------------

% slice increment when skipping images
nInc = 1;

% dynamic range of TIFF image ( 'uint8','uint16' )
dynRange = 'uint8';

% image format ( 'TIFF','JPG' )
iFormat = 'TIFF';

% image quality of JPG
jpgQuality = 80;


% -----------------------------------------------------------------------


% start and end slice
iStart = info.ProcessList.iStart( pStep );
iEnd   = info.ProcessList.iEnd( pStep );

% check consistency
if useFolder == true

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end

else

    niFile = iEnd - iStart + 1;
    [~,~,ni] = size( IMS ); 
    if ni~=niFile
        error('Number of images provided does not fit to size of image series variable')
    end

end


% factors for int16 -> int8
s8  = double(intmax( 'uint8' ));
s16 = double(intmax( 'uint16' ));


tic
for i=iStart:iEnd
    PrintIndex(i,iEnd);
    
    if useFolder == true

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderLoad,sliceName(i) );
        IT = imread( fname );

    else

        % current image
        k = i - iStart + 1;
        IT = IMS(:,:,k);

    end
    
    
    if strcmpi(iFormat,'TIFF')
        
        if strcmp(dynRange,'uint8') && ~isa(IT,'uint8')
           
            % make image 8 bit if not
            IT = uint8( double(IT) * s8 / s16 );
            
        end
        
        % save image 
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( IT ,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')
        
    else

        if ~isa(IT,'uint8')

            % make 8 bit
            IT8 = uint8( double(IT) * s8 / s16 );

        else

            % use as is
            IT8 = IT;
        end
        
        jpgName = convertStringsToChars(sliceName(i));
        jpgName(end-2:end) = 'jpg';
        
        % save back
        % ftest = fullfile( projectFolder,sliceFolder,'test.tif' );
        fname = fullfile( projectFolderSave,sliceFolderSave,jpgName );
        imwrite( IT8,fname,iFormat,'Quality',jpgQuality)

    end
end
toc
% __ sec


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%  1 . PREPARE . Crossbeam . Parse Atlas file list

% TASK  . extract filenams and slice thicknesses from filenames

% UPDATE -- WORKFLOW ----------------------------------------------------

processStep = 1;
processName = "PREPARE . Crossbeam . Parse file list";

processNote = "";


% UPDATE -- PROCESS -----------------------------------------------------

% data location
sliceIn.projectFolder = info.Project.storageFolder;
sliceIn.sliceFolder   = convertStringsToChars( info.Detector );


% -----------------------------------------------------------------------

fprintf( 'ProcessStep %u: Crossbeam . Parse file list ...\n',processStep );

% start timer
tic

% extract slice information encoded in file names
sliceOut = SerialZeissParseFileList( sliceIn );

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE  . update info structure

info.Project.Session            = info.Session;
info.Project.processStep        = processStep;
info.Project.processAppliedTo   = 0;
info.Project.Date               = timeFinished;
info.Project.Note               = processNote;

info.Project.nSlicesRecorded    = sliceOut.nSlicesRecorded;
info.Project.sliceName          = sliceOut.sliceName;
info.Project.slicePosition      = sliceOut.slicePosition;
info.Project.sliceThickness     = sliceOut.sliceThickness;
info.Project.sliceThicknessMean = sliceOut.sliceThicknessMean;
info.Project.sliceThicknessStd  = sliceOut.sliceThicknessStd;

% beginning start and end slices
info.Project.iStart = 1;
info.Project.iEnd = info.Project.nSlicesRecorded;
info.Project.iRef = round( 0.5*(info.Project.iStart + info.Project.iEnd) );


% SAVE  . generate process list table

info.ProcessList.Session    = info.Session;
info.ProcessList.Step       = processStep;
info.ProcessList.Process    = processName;
info.ProcessList.Save       = "";

info.ProcessList.iStart     = info.Project.iStart;
info.ProcessList.iEnd       = info.Project.iEnd;
info.ProcessList.nX         = 0;
info.ProcessList.nY         = 0;
info.ProcessList.nZ         = info.Project.nSlicesRecorded;

info.ProcessList.Data       = "Project";
info.ProcessList.Note       = processNote;
info.ProcessList.Date       = timeFinished;
info.ProcessList.Duration   = elapsedTime;


% info structure modified
info.Dates.Modified = timeFinished;

if autoSaveInfo == true

    % update time
    info.Dates.Saved = timeFinished;

    % save to file
    save( infoFile,'info','-v7.3' )

end


% PLOT  . Slice thicknesses

figure('Name','Slice Thickness','NumberTitle','off');
plot( 1000*info.Project.sliceThickness )
xlabel( 'slice' )
ylabel( 'thickness (nm) ' )


% CLEAR large temporary variables
clearvars sliceIn sliceOut

fprintf( '\n ... Finished successfully\n' );



%%  2 . PREPARE . Crossbeam . Extract Manually  >> Ext

% TASK  Extract image boundaries from representative images, define
%       global boundaries and crop manually

% INFO  Useful for runs where the image position is not or infrequently
%       changed or if a sub region is of interest only

% NEED  XROI must NEVER overlap with tracking marker regions if XROI is
%       searched automatically

% PROC  Divided into two steps
%       1. read representative images into variable IRS, measure
%          coordinates of largest area and define bounding box
%       2. extract regions

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number and name
processStep = 2;
processName = "PREPARE . Crossbeam . Extract XROI . Manual";

% process name acronym (3 characters, start with capital letter)
processSave = "Ext";

% process note
processNote = "";

% process step
subProcess = 2;


% UPDATE -- IN/OUT ------------------------------------------------------

% location of original images
projectFolderRead = info.Project.storageFolder;

% use folders (otherwise save images into variable)
useFolders = true;

if useFolders == true

    % output
    projectFolderSave = info.Project.projectFolder;

end


% UPDATE -- PROCESS -----------------------------------------------------

% fiducial markers in image region?
markerInImage = true;


% reference image indeces ( for sub process 1)
refImg.Idx = [20,1200];

% manual bounding box ( for sub process 2)
% origin: pixel top-left above image corner
% x: horizontal, y: vertical
ox = 1862;
oy = 2856;
wx = 4862 - ox;
wy = 5849 - oy;


% -----------------------------------------------------------------------

fprintf( 'Process Step %u: Crossbeam . Extract XROI . Manual ...\n', ...
    processStep );


% slice folder names
[sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
    info.Detector,processStep,info.ProcessList,processSave );

% generate folder
[SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
if SUCCESS == 0
    error( '%s (ID %s)',MESSAGE,MESSAGEID )
end

% slice names
sliceName = info.Project.sliceName;

% number of slices
nSlicesRecorded = info.Project.nSlicesRecorded;

% extract mode
extractMode = 'FixedBoundingBox';

% crop mode
cropMode = 'Manual';


if subProcess == 1

    % number of reference images
    nri = max(size( refImg.Idx ));

    % bounding boxes
    refImg.BoundingBox = zeros(4,nri);

    for i=1:nri
        PrintIndex( i,nri );

        % read file
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName( refImg.Idx(i) ));
        IT = imread( fname );

        if i==1

            % generate reference image series
            [ny,nx] = size( IT );
            IRS = zeros( ny,nx,nri,'like',IT );

        end

        % save image
        IRS(:,:,i) = IT;


        % find bounding box of largest region
        % !!! must not overlap with tracking marker windows

        regI       = regionprops( IT>0,'BoundingBox','Area' );
        [ ~,imax ] = max(cat( 1,regI.Area ));
        bb         = regI( imax ).BoundingBox;

        % generate proper origin
        bb(1:2) = floor( bb(1:2) );

        % save into structure
        refImg.BoundingBox(:,i) = bb';

    end
end


if subProcess == 2

    % start timer
    tic

    if useFolders ~= true

        % allocate memory
        IMS = zeros( wy,wx,nSlicesRecorded,'like',IRS );

    end

    for i=1:nSlicesRecorded

        if verbose == true
            PrintIndex( i,nSlicesRecorded );
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % crop
        ITC = IT( oy+1:oy+wy, ox+1:ox+wx );

        if useFolders == true

            % save back
            fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
            imwrite( ITC,fname,'tiff', ...
                'Compression','packbits','WriteMode','overwrite' )

        else

            IMS(:,:,i) = ITC;

        end
    end

    % end timer
    elapsedTime = toc;
    % 6465 sec

    % current time
    timeFinished = datetime( 'now' );


    % SAVE . update info structure

    info.Extract.Session          = info.Session;
    info.Extract.processStep      = processStep;
    info.Extract.processAppliedTo = processStep - 1;
    info.Extract.Date             = timeFinished;
    info.Extract.Note             = processNote;

    info.Extract.extractMode    = extractMode;
    info.Extract.cropMode       = cropMode;
    info.Extract.markerInImage  = markerInImage;
    info.Extract.RefImages      = refImg;
    info.Extract.BoundingBox    = [ox,oy,wx,wy];

    if exist('IRS','var')
        info.Extract.ImageClass     = class( IRS );
    else
        % last opened image
        info.Extract.ImageClass     = class( IT );
    end


    % SAVE  . update process list

    newProcess            = table;
    newProcess.Session    = info.Session;
    newProcess.Step       = processStep;
    newProcess.Process    = processName;
    newProcess.Save       = processSave;

    newProcess.iStart     = info.Project.iStart;
    newProcess.iEnd       = info.Project.iEnd;
    newProcess.nX         = wx;
    newProcess.nY         = wy;
    newProcess.nZ         = info.Project.nSlicesRecorded;

    newProcess.Data       = "Extract";
    newProcess.Note       = processNote;
    newProcess.Date       = timeFinished;
    newProcess.Duration   = elapsedTime;

    info.ProcessList(processStep,:) = newProcess;


    % info structure modified
    info.Dates.Modified = timeFinished;

    if autoSaveInfo == true

        % update time
        info.Dates.Saved = timeFinished;

        % save to file
        save( infoFile,'info','-v7.3' )

    end

end

fprintf( '\n ... Finished successfully\n' );


%%  3 . HISTOGRAM . Evaluate series

% TASK  Calculate histogram for full series

% NOTE  It is possible to
%       - calculate histogram for any given processing stage for which
%         data were saved. Either load the image series into 'IMS' 
%         variable or provide 'sliceFolderRead' link.
%       - skip saving the results, if only data exploration is the goal


% UPDATE -- WORKFLOW ----------------------------------------------------

% process step / name
processStep = 3;
processName = "HISTOGRAM . Histogram series";

% process note
processNote = "Original";


% save to info structure ? (true/false)
% (set to 'false' if just histogram calculation desired )
saveToInfo = true;

% plot results on screen (true/false)
plotResults = true;

% save results into 'Figures' folder (true/false)
saveResults = true;


% UPDATE -- IN/OUT ------------------------------------------------------

% use data from process step (number, processStep-1 )
iProc = processStep - 1;

% use folders
useFolders = true;

if useFolders == true

    % input
    projectFolderRead = info.Project.projectFolder;

end


% UPDATE -- PROCESS -----------------------------------------------------

% histogram parameters

switch info.Extract.ImageClass
    case 'uint8'
        hpar.nbin = 256;
    case 'uint16'
        hpar.nbin = 2048;
    otherwise
        error( 'Unknown image class.' )
end

hpar.mima   = [0,double(intmax( info.Extract.ImageClass ))];
hpar.YNplot = 'no';


% intensity state of the current series
% ( Original, Denoised, Equilibrated, ... )
intensityState = "Denoised";


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Calculate histogram series\n', processStep );

% start timer
tic

% slice names
sliceName = info.Project.sliceName;

% slice folder name from process list
[sliceFolderRead,~ ] = SerialSliceFolderNames( ...
    info.Detector,iProc+1,info.ProcessList,"" );


% image stack of used process task
iStart = info.ProcessList.iStart(iProc);
iEnd   = info.ProcessList.iEnd(iProc);

% number of slices in series
nSlices = iEnd - iStart + 1;

if useFolders == false

    % check consistency
    [~,~,ni] =size( IMS );
    if ni ~= nSlices
        error( 'Provided number of slices different from size of IMS' )
    end

end

% initialize histogram series
HS = zeros(nSlices,hpar.nbin);


if useFolders == true

    % read first slice image
    fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(iStart) );
    IT = imread( fname );

    % histogram
    [ival, nval] = TomoHist( IT,hpar );

    % save
    HS(1,:) = nval;


    parfor i=2:nSlices

        if verbose == true
            PrintIndex(i,nSlices);
        end

        % read current slice image
        k = iStart + i - 1;
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(k) ); %#ok<PFBNS> 
        IT = imread( fname );

        % histogram
        [~, nval] = TomoHist( IT,hpar );

        % save
        HS(i,:) = nval;

    end

else

    % intensity values from first histogram
    [ival, nval] = TomoHist( IMS(:,:,1),hpar );

    % save
    HS(1,:) = nval;


    parfor i=2:nSlices

        if verbose == true
            PrintIndex(i,nSlices);
        end

        % histogram
        [~, nval] = TomoHist( IMS(:,:,i),hpar );

        % save
        HS(i,:) = nval;

    end

end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% number of histogram calculation events so far
nHist = 0;
if isfield( info,'Histogram' )
    nHist = max(size( info.Histogram ));
end


if saveToInfo == true

    % SAVE . update info structure

    info.Histogram( nHist+1 ).Session          = info.Session;
    info.Histogram( nHist+1 ).processStep      = processStep;
    info.Histogram( nHist+1 ).processAppliedTo = iProc;
    info.Histogram( nHist+1 ).Data             = timeFinished;
    info.Histogram( nHist+1 ).Note             = processNote;

    if useFolders == true
        info.Histogram( nHist+1 ).sliceFolderRead = sliceFolderRead;
    end

    info.Histogram( nHist+1 ).IntensityState   = intensityState;
    info.Histogram( nHist+1 ).hpar             = hpar;
    info.Histogram( nHist+1 ).ival             = ival;
    info.Histogram( nHist+1 ).HS               = HS;


    % SAVE . update process list

    newProcess            = table;
    newProcess.Session    = info.Session;
    newProcess.Step       = processStep;
    newProcess.Process    = processName;
    newProcess.Save       = "";

    newProcess.iStart     = info.ProcessList.iStart(iProc);
    newProcess.iEnd       = info.ProcessList.iEnd(iProc);
    newProcess.nX         = info.ProcessList.nX(iProc);
    newProcess.nY         = info.ProcessList.nY(iProc);
    newProcess.nZ         = info.ProcessList.nZ(iProc);

    newProcess.Data       = "Histogram";
    newProcess.Note       = processNote;     
    newProcess.Date       = timeFinished;
    newProcess.Duration   = elapsedTime;

    info.ProcessList(processStep,:) = newProcess;


    % info structure modified
    info.Dates.Modified = timeFinished;

    if autoSaveInfo == true

        % update time
        info.Dates.Saved = timeFinished;

        % save to file
        save( infoFile,'info','-v7.3' )

    end

end


% PLOT

if plotResults == true

    plim(HS)

end


if saveResults == true

    % movie file
    fname = sprintf('S%u_P%u_Histogram%u.avi', ...
        info.Session,processStep,nHist+1 );
    fname = fullfile( info.Project.projectFolder,'Figures',fname );


    % start video
    v1 = VideoWriter( fname );
    open(v1);
    figure

    for i=1:nSlices
        PrintIndex(i,nSlices);

        % plot data
        plot( ival,HS(i,:) )

        % add frame to video
        frame = getframe;
        writeVideo(v1,frame);

    end

    % stop video
    close(v1);

end

fprintf( '\n ... Finished successfully\n' );



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%  4 . DENOISE - Estimate noise level ...

% TASK  find images with areas containing flat intensity and estimate noise

% PROC  SubProcess 1
%           1. read desired images (SubProcess 1) into cell array 'ims'
%           2. plot each image via plim(ims{i})
%           3. select region with crop tool, right-click on region, select
%              'Copy Region'
%           4. Fill p variable, remove after-comma values
%
%       SubProcess 2
%           1. estimate noise

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number / name
processStep = 4;
processName = "DENOISE . Estimate noise level";

% process note
processNote = "";

% sub process
% 1 load test images and select regions with constant intensity
% 2 extract noise information from those regions
subProcess = 2;


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders
useFolders = true;

if useFolders == true

    % input
    projectFolderRead = info.Project.projectFolder;

end


% UPDATE -- PROCESS -----------------------------------------------------

% test images
% (filename index if folders are used, otherwise image series index)
idx = [680,1307];

% uncomment if above indeces are given in cropped coordinates: 
% transform back to indices of raw dataset
% idx = idx + info.ProcessList.iStart( processStep-1 ) - 1;

% fit degree of polynomial to remove long-range intensity variation
degree = 4;

if subProcess == 2

    % (index 1: test image, index 2: region with near-constant intensity)
    p(1,1,:) = [35 1999 147 250];
    p(2,1,:) = [51 1408 173 203];
%     p(2,2,:) = [1130 1778 71 59];
%     p(2,3,:) = [1227 1975 99 90];

end


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Estimate noise from patches ...\n',processStep );

if useFolders == true

    % slice name
    sliceName = info.Project.sliceName;

    % input folder
    [sliceFolderRead,~ ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

end


% READ IMAGES

if subProcess == 1

    % number of test images
    nIdx = max(size( idx ));

    % cell array for storage (in case images have different size)
    ims = cell( nIdx,1 );
    k   = 0;

    for i=idx
        PrintIndex(i,max(idx));

        % running index
        k=k+1;

        if useFolders == true

            % read test images
            fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
            ims{k} = double( imread( fname ));

        else

            ims{k} = IMS(:,:,i);

        end

    end

end


% EXTRACT AREAS WITH CONSTANT INTENSITY

if subProcess == 2

    % start timer
    tic

    % number of images and background regions
    szp  = size( p );
    nImg = szp(1);
    nBg  = szp(2);

    % extracted ROIs
    rois = cell( nImg,nBg,1 );

    % mean and standard deviations
    roiMean = zeros(nImg,nBg);
    roiStd  = zeros(nImg,nBg);

    for k=1:nImg
        PrintIndex( k,nImg );

        % current test image
        IT = ims{ k };

        for i=1:nBg

            % extract region
            ROI = IT( p(k,i,2)+1:p(k,i,2)+p(k,i,4), ...
                p(k,i,1)+1:p(k,i,1)+p(k,i,3) );

            % remove image gradients by 2D polynomial fit
            [NFIT,NRES] = SerialNoiseRoiFit( ROI, degree );

            % mean
            roiMean(k,i) = mean( ROI(:) );
            roiStd(k,i)  = std( ROI(:)-NFIT(:) );

            % save area
            rois{k,i} = cat( 2,ROI,NFIT,NRES+mean2(NFIT) );

        end
    end

    % stop timer
    elapsedTime = toc;

    % current time
    timeFinished = datetime( 'now' );


    
    % SAVE . update info structure

    iProc                         = processStep-1;

    info.Denoise.Session          = info.Session;
    info.Denoise.processStep      = processStep;
    info.Denoise.processAppliedTo = iProc;
    info.Denoise.Date             = timeFinished;
    info.Denoise.Note             = processNote;

    info.Denoise.Roi.imgIdx       = idx;
    info.Denoise.Roi.position     = p;
    info.Denoise.Roi.fitDegree    = degree;
    info.Denoise.Roi.fits         = rois;
    info.Denoise.ROI.mean         = roiMean;
    info.Denoise.ROI.std          = roiStd;


    % SAVE . update process list

    newProcess            = table;
    newProcess.Session    = info.Session;
    newProcess.Step       = processStep;
    newProcess.Process    = processName;
    newProcess.Save       = "";

    newProcess.iStart     = info.ProcessList.iStart(iProc);
    newProcess.iEnd       = info.ProcessList.iEnd(iProc);
    newProcess.nX         = info.ProcessList.nX(iProc);
    newProcess.nY         = info.ProcessList.nY(iProc);
    newProcess.nZ         = info.ProcessList.nZ(iProc);

    newProcess.Data       = "Denoise";
    newProcess.Note       = processNote;
    newProcess.Date       = timeFinished;
    newProcess.Duration   = elapsedTime;

    info.ProcessList(processStep,:) = newProcess;


    % info structure modified
    info.Dates.Modified = timeFinished;

    if autoSaveInfo == true

        % update time
        info.Dates.Saved = timeFinished;

        % save to file
        save( infoFile,'info','-v7.3' )

    end


    % PLOT

    % all noise levels compiled
    roimean = roiMean(:);
    roistd = roiStd(:);

    % standard deviation over mean (compare with Possion)
    figure
    plot( roiMean',roiStd','+-' )
    xlabel( 'mean intensity' )
    ylabel( 'standard deviation ' )

    % ROI areas with fit and residuals
    for k=1:nImg
        for i=1:nBg
            plim( rois{k,i} )
        end
    end

    fprintf( '\n ... Finished successfully\n' );

end


%% ... DENOISE . Test Parameter Space

% TASK  Check quality of denoising as function of algorithm parameters to
%       make sure that images are not over-processed and residuals do not 
%       contain any systematic intensity features

% NEED  cell array ims from previous section

% PROC  1. Choose one of the images loaded in previous session
%       2. Select small representative region
%       3. Calculate denoised image as function of denoising parameters

% UPDATE -- Process -----------------------------------------------------

% test images in cell array variable 'ims'
iTest = 1;

% bounding box of ROI
p = [1063 1400 744 765];

% non-local means parameter space
% noise level, probe window, averaging window half-width
nlmNoise = [1200,2400,3600];
nlmProbe = [9,7,5];
nlmAvg   = [3,2,1];

% -----------------------------------------------------------------------

% number of parameter combinations
nNoise = max(size(nlmNoise));
nProbe = max(size(nlmProbe));
nAvg   = max(size(nlmAvg));
nPar   = nNoise * nProbe * nAvg; 

% fill parameter space
nlmPar = zeros(nPar,3);

% running index
i = 1;

for k=1:nNoise
    for l=1:nProbe
        for m=1:nAvg

            % parameter set
            nlmPar(i,:) = [ nlmNoise(k),nlmProbe(l),nlmAvg(m) ];

            % update k
            i = i + 1;
        
        end
    end
end


% test region
IT = ims{ iTest };
IT = IT( p(2)+1:p(2)+p(4) , p(1)+1:p(1)+p(3) );

% denoised image series and residuals
[ny,nx] = size( IT );
IDS     = zeros( ny,nx, nPar,'like',IT );
IRS     = zeros( ny,nx, nPar,'like',IT );


parfor i=1:nPar
    PrintIndex(i,nPar);

    % denoise
    IDS(:,:,i) = uint16( FilterNLMeans_mex( double( IT ), ...
            nlmPar(i,2),nlmPar(i,3),nlmPar(i,1) )); %#ok<PFBNS> 

    % residual
    IRS(:,:,i) = IT - IDS(:,:,i);

end


% SAVE . IDS AND IRS mrc files

fImg = sprintf( 'S%u_P%u_NoiseParameterTest_Image.mrc',...
    info.Session,processStep);
fImg = fullfile(info.Project.projectFolder,'Figures',fImg );

fRes = sprintf( 'S%u_P%u_NoiseParameterTest_Residual.mrc',...
    info.Session,processStep);
fRes  = fullfile(info.Project.projectFolder,'Figures',fRes );

TomoSaveMRC( IDS,fImg,6 );
TomoSaveMRC( IRS,fRes,6 );


% PLOT . display parameter sets

nlmIdx  = transpose( 1:nPar );
nlmShow = cat(2,nlmIdx,nlmPar );
disp( nlmShow )


% SAVE to info structure

info.Denoise.Test.nlmNoise = nlmNoise;
info.Denoise.Test.nlmProbe = nlmProbe;
info.Denoise.Test.nlmAvg   = nlmAvg;

% table of image indeces and corresponding 
nlmTable = array2table( nlmShow, ...
    'VariableNames',{'Image','Noise','Probe','Avg'});
info.Denoise.Test.ParameterList = nlmTable;


if autoSaveInfo == true

    % update time
    info.Dates.Saved = datetime( 'now' );

    % save to file
    save( infoFile,'info','-v7.3' )

end



%%  5 . DENOISE . All images  >> Den

% TASK  Remove noise (using non-local means, Buades et. al.)

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number / name
processStep = 5;
processName = "DENOISE . All images";

% process name acronym (3 characters, start with capital letter)
processSave = "Den";

% process note
processNote = "";


% UPDATE -- Input/Output ------------------------------------------------

% use folders (true,false)
useFolders = true;

if useFolders == true

    projectFolderRead      = info.Project.projectFolder;
    projectFolderSave      = projectFolderRead;

end


% UPDATE -- Process -----------------------------------------------------

% noise level
roiStd = 2400;

% denoise type
denoiseType       = 'Non-Local Means';
probeWindow2   = 5;
averageWindow2 = 2;


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Denoise all iamges ...\n', processStep );

% start timer
tic


if useFolders == true

    % slice names
    sliceName  = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end

    iStart = info.ProcessList.iStart( processStep-1);
    iEnd   = info.ProcessList.iEnd( processStep-1);


    parfor i=iStart:iEnd
        % parfor i=1:5

        if verbose == true
            PrintIndex(i,iEnd);
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % remember class
        imgClass = class( IT ); %#ok<PFTUSW> 

        % non-local means
        ID = FilterNLMeans_mex( double( IT ), ...
            probeWindow2,averageWindow2,roiStd );

        % save
        switch imgClass
            case 'uint8'
                fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
                imwrite( uint8(ID),fname,'tiff', ...
                    'Compression','packbits','WriteMode','overwrite')

            case 'uint16'
                fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
                imwrite( uint16(ID),fname,'tiff', ...
                    'Compression','packbits','WriteMode','overwrite')

        end
    end
end


if useFolders == false

    if exist( 'IMSP','var' )

        % reset image series variable from previous process
        IMS = IMSP;
        clearvars IMSP

    end
    
    % reinitialize resulting image series
    IMSP = zeros( size(IMS),'like',IMS );

    % image class
    imgClass = class( IMS );

    % number of images
    [~,~,nSlices] = size( IMS );


    parfor i=1:nSlices
        % parfor i=1:5

        if verbose == true
            PrintIndex(i,nSlices);
        end

        % non-local means
        ID = FilterNLMeans_mex( double( IMS(:,:,i) ), ...
            probeWindow2,averageWindow2,roiStd );
        
        % save
        switch imgClass
            case 'uint8'
                IMSP(:,:,i) = uint8( ID );

            case 'uint16'
                IMSP(:,:,i) = uint16( ID );
        end
    end
end

% stop timer
elapsedTime = toc;
% 1025 sec (40 cores)

% current time
timeFinished = datetime( 'now' );


% SAVE .  update info structure

info.Denoise.Session             = info.Session;
info.Denoise.processStep         = processStep;
info.Denoise.processAppliedTo    = processStep - 1;
info.Denoise.Date                = timeFinished;
info.Denoise.Note                = processNote;

info.Denoise.type                = denoiseType;
info.Denoise.Nlm.probeWindow2    = probeWindow2;
info.Denoise.Nlm.averageWindow2  = averageWindow2;
info.Denoise.roiStd              = roiStd;

info.Denoise.sliceFolderRead     = sliceFolderRead;
info.Denoise.sliceFolderSave     = sliceFolderSave;
info.Denoise.iStart              = iStart;
info.Denoise.iEnd                = iEnd;


% SAVE . update process list

iProc                 = processStep-1;

newProcess            = table;
newProcess.Session    = info.Session;
newProcess.Step       = processStep;
newProcess.Process    = processName;
newProcess.Save       = processSave;

newProcess.iStart     = info.ProcessList.iStart(iProc);
newProcess.iEnd       = info.ProcessList.iEnd(iProc);
newProcess.nX         = info.ProcessList.nX(iProc);
newProcess.nY         = info.ProcessList.nY(iProc);
newProcess.nZ         = info.ProcessList.nZ(iProc);

newProcess.Data       = "Denoise";
newProcess.Note       = processNote;
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

info.ProcessList(processStep,:) = newProcess;


% info structure modified
info.Dates.Modified = timeFinished;

if autoSaveInfo == true

    % update time
    info.Dates.Saved = timeFinished;

    % save to file
    save( infoFile,'info','-v7.3' )

end


fprintf( '\n ... Finished successfully\n' );



%% ... PLOT - Compare

% TASK show degree of denoising in a representative small region

% UPDATE ----------------------------------------------------------------

%in/out
loadFolder        = info.Project.projectFolder;
sliceName         = info.Project.sliceName;

sliceFolderBefore = 'ESB_ExtCro';
sliceFolderAfter  = 'ESB_ExtCroDen';

% test slice
iTest = info.ProcessList.iStart( processStep) + 364;

% iamge region
xS = 1288;
xE = 1898;
yS = 1663;
yE = 2231;


% -----------------------------------------------------------------------

% read images
fname = fullfile( loadFolder,sliceFolderBefore,sliceName(iTest) );
IORIG = double( imread( fname ));
IORIG = IORIG(yS:yE,xS:xE);

fname = fullfile( loadFolder,sliceFolderAfter,sliceName(iTest) );
IDEN  = double( imread( fname ));
IDEN  = IDEN(yS:yE,xS:xE);


% difference
IDIF = IORIG - IDEN;

% scale to image mean
IDIF = IDIF + mean2( IORIG );

% collage
ICOL = cat( 2,IORIG,IDEN,IDIF );


% adjust contrast based on histogram
hpar.YNplot = 'no';
[iOrig,nOrig] = TomoHist( IORIG,hpar );
[iDen,nDen]   = TomoHist( IDEN,hpar );

figure
hold on, plot( iOrig,nOrig );
hold on, plot( iDen,nDen );
hold off
xlabel( 'intensity bin' )
ylabel( 'number of pixels' )
legend('original','denoised')
xlim( [0,40000])

ICOL( ICOL>40000 ) = 40000;

% PLOT
plim( ICOL )


% SAVE . collage to file
fCol = sprintf( 'S%u_P%u_Denoising_ROI_compare.tif',...
    info.Session,processStep);
fCol = fullfile( info.Project.projectFolder,'Figures',fCol );

SaveImage( ICOL,fCol,'TIF' );




%%  6 . HISTOGRAM . Evaluate series

% TASK  Calculate histogram for full series

% NOTE  It is possible to
%       - calculate histogram for any given processing stage for which
%         data were saved. Either load the image series into 'IMS' 
%         variable or provide 'sliceFolderRead' link.
%       - skip saving the results, if only data exploration is the goal


% UPDATE -- WORKFLOW ----------------------------------------------------

% process step / name
processStep = 6;
processName = "HISTOGRAM . Histogram series";

% process note
processNote = "Denoised";


% save to info structure ? (true/false)
% (set to 'false' if just histogram calculation desired )
saveToInfo = true;

% plot results on screen (true/false)
plotResults = true;

% save results into 'Figures' folder (true/false)
saveResults = true;


% UPDATE -- IN/OUT ------------------------------------------------------

% use data from process step (number, processStep-1 )
iProc = processStep - 1;

% use folders
useFolders = true;

if useFolders == true

    % input
    projectFolderRead = info.Project.projectFolder;

end


% UPDATE -- PROCESS -----------------------------------------------------

% histogram parameters

switch info.Extract.ImageClass
    case 'uint8'
        hpar.nbin = 256;
    case 'uint16'
        hpar.nbin = 2048;
    otherwise
        error( 'Unknown image class.' )
end

hpar.mima   = [0,double(intmax( info.Extract.ImageClass ))];
hpar.YNplot = 'no';


% intensity state of the current series
% ( Original, Denoised, Equilibrated, ... )
intensityState = "Denoised";


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Calculate histogram series\n', processStep );

% start timer
tic

% slice names
sliceName = info.Project.sliceName;

% slice folder name from process list
[sliceFolderRead,~ ] = SerialSliceFolderNames( ...
    info.Detector,iProc+1,info.ProcessList,"" );


% image stack of used process task
iStart = info.ProcessList.iStart(iProc);
iEnd   = info.ProcessList.iEnd(iProc);

% number of slices in series
nSlices = iEnd - iStart + 1;

if useFolders == false

    % check consistency
    [~,~,ni] =size( IMS );
    if ni ~= nSlices
        error( 'Provided number of slices different from size of IMS' )
    end

end

% initialize histogram series
HS = zeros(nSlices,hpar.nbin);


if useFolders == true

    % read first slice image
    fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(iStart) );
    IT = imread( fname );

    % histogram
    [ival, nval] = TomoHist( IT,hpar );

    % save
    HS(1,:) = nval;


    parfor i=2:nSlices

        if verbose == true
            PrintIndex(i,nSlices);
        end

        % read current slice image
        k = iStart + i - 1;
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(k) ); %#ok<PFBNS> 
        IT = imread( fname );

        % histogram
        [~, nval] = TomoHist( IT,hpar );

        % save
        HS(i,:) = nval;

    end

else

    % intensity values from first histogram
    [ival, nval] = TomoHist( IMS(:,:,1),hpar );

    % save
    HS(1,:) = nval;


    parfor i=2:nSlices

        if verbose == true
            PrintIndex(i,nSlices);
        end

        % histogram
        [~, nval] = TomoHist( IMS(:,:,i),hpar );

        % save
        HS(i,:) = nval;

    end

end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% number of histogram calculation events so far
nHist = 0;
if isfield( info,'Histogram' )
    nHist = max(size( info.Histogram ));
end


if saveToInfo == true

    % SAVE . update info structure

    info.Histogram( nHist+1 ).Session          = info.Session;
    info.Histogram( nHist+1 ).processStep      = processStep;
    info.Histogram( nHist+1 ).processAppliedTo = iProc;
    info.Histogram( nHist+1 ).Data             = timeFinished;
    info.Histogram( nHist+1 ).Note             = processNote;

    if useFolders == true
        info.Histogram( nHist+1 ).sliceFolderRead = sliceFolderRead;
    end

    info.Histogram( nHist+1 ).IntensityState   = intensityState;
    info.Histogram( nHist+1 ).hpar             = hpar;
    info.Histogram( nHist+1 ).ival             = ival;
    info.Histogram( nHist+1 ).HS               = HS;


    % SAVE . update process list

    newProcess            = table;
    newProcess.Session    = info.Session;
    newProcess.Step       = processStep;
    newProcess.Process    = processName;
    newProcess.Save       = "";

    newProcess.iStart     = info.ProcessList.iStart(iProc);
    newProcess.iEnd       = info.ProcessList.iEnd(iProc);
    newProcess.nX         = info.ProcessList.nX(iProc);
    newProcess.nY         = info.ProcessList.nY(iProc);
    newProcess.nZ         = info.ProcessList.nZ(iProc);

    newProcess.Data       = "Histogram";
    newProcess.Note       = processNote;     
    newProcess.Date       = timeFinished;
    newProcess.Duration   = elapsedTime;

    info.ProcessList(processStep,:) = newProcess;


    % info structure modified
    info.Dates.Modified = timeFinished;

    if autoSaveInfo == true

        % update time
        info.Dates.Saved = timeFinished;

        % save to file
        save( infoFile,'info','-v7.3' )

    end

end


% PLOT

if plotResults == true

    plim(HS)

end


if saveResults == true

    % movie file
    fname = sprintf('S%u_P%u_Histogram%u.avi', ...
        info.Session,processStep,nHist+1 );
    fname = fullfile( info.Project.projectFolder,'Figures',fname );


    % start video
    v1 = VideoWriter( fname );
    open(v1);
    figure

    for i=1:nSlices
        PrintIndex(i,nSlices);

        % plot data
        plot( ival,HS(i,:) )

        % add frame to video
        frame = getframe;
        writeVideo(v1,frame);

    end

    % stop video
    close(v1);

end

fprintf( '\n ... Finished successfully\n' );



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Serial section image series recorded on the Crossbeam 550 FIB show 
% scan distortions which cannot be described with a simple shear transform. 
% Instead of showing a linear dependence, they rather vary as function 
% of the slow scan scan direction (y) showing independent displacements 
% along the x- and y-direction.
% 
% The rationale behind the following group of process steps is that the
% image acquisition can be divided into a fast scan along the horizontal 
% (x) on order of milliseconds and a slow san along the vertical (y) on
% order seconds or minutes. Furthermore it can be be observed that 
%    (a) the absolute values of the scan distortions are samll
%        (a few pixels)
%    (b) vary slowly along the y-direction
%
% The above conditions allow to split the image into a set of thin
% horizontal stripes. The function SerialCrossCorrSeries() is applied to
% each individual stripe. For each stripe the same spline fit is used as
% in the above case of the global drift correction, with the goal to
% remove the artificial evolution of the cummulative drift and regain the
% actual image-to-image displacements.
%
% NOTE: A necessary condition for the successful use of the above proposed
%       stripe-based cross correlation approach is that each stripe 
%       contains a large enoug number of objects to average over
%       individual displacements that might be artificially introduced 
%       due to the shape change of an object as the cross sections
%       traverse the object. Fortunately, past experiments indicate that
%       the images can be rather coarse in nature with few objects to still
%       produce reliable displacement values from the cross correlation.
%
% Stripe heights on order 100-200 pixels have been employed and appear to
% sample the scan distortions on a fine enough grid.
%
% The obtained short range distortions (dx,dy) are collected as function
% of the center of each stripe along the y-direction and interpolated to
% pixel-resolution. The corresponding section produces movies for both x
% and y-displacements that allow to cross check whether the y-dependence
% is smooth. So far, this has been the case for all analyzed serial section
% runs indicating that the above mentioned conditions appear to be valid.
%
% For the actual corrections of the scan distortion it is possible to 
% employ the Image Toolbox function imwarp(), which uses a displacement 
% field that can be constructed from the interpolated local scan 
% distortions.



%% ... Y-DEP DRIFT . Optimum slice thickness

% TASK  Find slice thickness that keeps the residual thickness to a
%       reasonable minimum

% UPDATE -- PROCESS -----------------------------------------------------

% range of 0.5 stripe heights
sHalf   = 45:150;


% -----------------------------------------------------------------------


% optimum stripe thickness

% image height
nY = info.ProcessList.nY( processStep );

% stripe height
sHeight = 2*sHalf + 1;

% number of stripes
nS = floor( nY./sHeight );

% residual
nRes = nY - nS.*sHeight;

figure
plot( sHeight,nRes,'b+-' )
xlabel( 'stripe height' )
ylabel( 'residual pixel count' )



%%   7 . Y-DEP DRIFT . Cross correlate stripes

% TASK  calculate drift for each stripe

% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 7;
processName = "Y-DEP DRIFT . Croos correlate stripes";

% process Note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders (true,false)
useFolders = true;

if useFolders == true

    projectFolderRead      = info.Project.projectFolder;

end


% UPDATE -- PROCESS -----------------------------------------------------

% CONTROL PARAMETERS

% Flags
ccIn.YNGauss           = 'yes';  % smooth with Gaussian
ccIn.YNgrad            = 'yes';  % use gradient image for cc    
ccIn.YNpartial         = 'no';  % use reduced area for cc (not here)
ccIn.YNleastCommonArea = 'no';  % use only common non-0 area for cc    *1
ccIn.YNfit             = 'yes';  % apply polynomial fit to cc peak      *2

% *1 in case images were shifted during the run
%    CAREFUL: stripe height needs to be significantly larger than applied
%             y-shift during run
% *2 two polynomial fits are applied
%    fit 1: general fit, plotted over given sub-pixel grid
%    fit 2: peak postion is parameterized, results of fit 1 are used as
%           starting values

if strcmpi( ccIn.YNGauss,'yes' )

    % Gaussian kernel parameters
    ccIn.Gauss.h = 10;
    ccIn.Gauss.s = 3;
end

if strcmpi( ccIn.YNpartial,'yes' )

    % use cropped region for cross correlation
    % (indeces must fit into the confines of the used images)
    ccIn.Crop.ixS = 0;
    ccIn.Crop.ixE = 0;
    ccIn.Crop.iyS = 0;
    ccIn.Crop.iyE = 0;

end


% reference image for global cummulative drift
ccIn.idxRef = 730;

% image for which peak fit is plotted
% ( value > iStart, using the value 0 suppresses the plot )
% ( here too many calls, set to 0 )
ccIn.idxPlot = 0;


% NOT YET IMPLEMENTED
% % half-width search window for peak finding
% % concentrate the peak search on region near center of cc, in case
% there are spurious peaks further away due to imperfect image conditions
% ccIn.s2 = 50;                       


% CROSS CORRELATION PARAMETERS ( for FastCrossCorr function )

% mode of operation: image-to-image shift (don't change)
ccIn.mode = 'shift'; 

% taper width of image (set to 0, does not apply here)
ccIn.w = 0;             

% sub-pixel grid over which fit 1 is calculated
ccIn.fit.dpf = 0.01;    

% peak fit parameters for fit 2
ccIn.fit.deg = 4;
ccIn.fit.sx2 = 3;
ccIn.fit.sy2 = 3;
ccIn.plotfitnum = 2;    % show either fit 1 or 2 if plotted

% which fit should be used
ccIn.usefit = 2;


% STRIPES

% half height
sHalf = 99;

% NOT IMPLEMENTED (for now no overlap between stripes)
% % increment 
% sInc = 201;


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Y-DEP DRIFT . Cross correlate stripes ...\n', processStep );

% start timer
tic


% PARAMETERS

% update input structure for SerialCrossCorrSeries()

if useFolders == true

    % slice names
    sliceName = info.Project.sliceName;

    % input slice folder
    [sliceFolderRead,~ ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % start/end slice
    iStart = info.ProcessList.iStart( processStep-1 );
    iEnd   = info.ProcessList.iEnd( processStep-1 );

    ccIn.inputType           = 'files';
    ccIn.Input.iStart        = iStart;
    ccIn.Input.iEnd          = iEnd;
    ccIn.Input.projectFolder = projectFolderRead;
    ccIn.Input.sliceFolder   = sliceFolderRead;
    ccIn.Input.sliceName     = sliceName;

else

    ccIn.inputType = 'variable';

end


% stripe height
sHeight = 2*sHalf + 1;

% current image height
if useFolders == true

    ny = info.ProcessList.nY( processStep-1 );

else
    [ny,~,~] = size(IMS);

end

% number of stripes 
nStripes = floor(ny/sHeight);

% stripe centers
si      = 0:nStripes-1;
sCenter = sHalf+1 + sHeight*si;

% for collecting all cross correlation results
ccList = cell(nStripes,1);


% CALCULATE CROSS CORRELATIONS

if useFolders == true

    % generate temporary image folder for each stripe
    
    for i=1:nStripes

        if verbose == true
            PrintIndex(i,nStripes);
        end

        % current stripe folder name
        stripeFolder = sprintf( 'Stripe_%03d',i );

        % generate folder
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderRead,stripeFolder );
        if SUCCESS == 0
            error( '%s (ID %s)',MESSAGE,MESSAGEID )
        end

    end


    % save stripes for each image

    for k=iStart:iEnd

        if verbose == true
            PrintIndex(k,iEnd);
        end

        % read current image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(k) );
        IT = imread( fname );

        for i=1:nStripes

            % extract slice
            S = IT( sCenter(i)-sHalf:sCenter(i)+sHalf,: );

            % corresponding stripe folder
            stripeFolder = sprintf( 'Stripe_%03d',i );

            % save image
            fname = fullfile( projectFolderRead,stripeFolder,sliceName(k) );
            imwrite( S,fname,'tiff', ...
                'Compression','packbits','WriteMode','overwrite')

        end
    end
 

    % cross correlation

    parfor i=1:nStripes
    % for i=1:1

    if verbose == true
        PrintIndex(i,nStripes);
    end

        % set file I/O
        ccInT = ccIn;
        ccInT.Input.sliceFolder = sprintf( 'Stripe_%03d',i );

        % do cross correlation
        ccList{i} = SerialCrossCorrSeries( ccInT );
       
    end


% FOR NOW: DELETE STRIPE FOLDERS IN WINDOWS
%     for i=1:nStripes
% 
%         % delete temporary folders
% 
% 
%     end
    
else

    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    % rearrange image series for parfor loop
    [~,nx,ni] = size( IMS );
    S = zeros( sHeight,nx,ni,nStripes,'like',IMS );

    for i=1:nStripes
        S(:,:,:,i) = IMS( sCenter(i)-sHalf:sCenter(i)+sHalf,:,: );
    end


    parfor i=1:nStripes
        % for i=1:1

        if verbose == true
            PrintIndex(i,nStripes);
        end

        % general cross correlation
        ccList{i} = SerialCrossCorrSeries( ccIn,S(:,:,:,i) );

    end

end

elapsedTime = toc;
% __ sec

% current time
timeFinished = datetime( 'now' );



% SAVE . update info structure

info.YdepDrift.Session          = info.Session;
info.YdepDrift.processStep      = processStep;
info.YdepDrift.processAppliedTo = processStep - 1;
info.YdepDrift.Date             = timeFinished;
info.YdepDrift.Note             = processNote;

info.YdepDrift.sHalf  = sHalf;
info.YdepDrift.ccIn   = ccIn;
info.YdepDrift.ccList = ccList;


% SAVE . add process step

newProcess            = table;
newProcess.Session    = info.Session;
newProcess.Step       = processStep;
newProcess.Process    = processName;
newProcess.Save       = "";

newProcess.iStart     = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd       = info.ProcessList.iEnd( processStep-1 );
newProcess.nX         = info.ProcessList.nX( processStep-1 );
newProcess.nY         = info.ProcessList.nY( processStep-1 );
newProcess.nZ         = info.ProcessList.nZ( processStep-1 );

newProcess.Data       = "YdepDrift";
newProcess.Note       = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

info.ProcessList(processStep,:) = newProcess;


% info structure modified
info.Dates.Modified = timeFinished;

if autoSaveInfo == true

    % update time
    info.Dates.Saved = timeFinished;

    % save to file
    save( infoFile,'info','-v7.3' )

end


fprintf( '\n ... Finished successfully\n' );



%% ... PLOT . Cummulative drift

% TASK  overlay drift values for first and last stripe

% UPDATE ----------------------------------------------------------------
% -----------------------------------------------------------------------

% number of stripes
nStripes = max(size( info.YdepDrift.ccList ));

% read cummulative drifts
cdrift1 = info.YdepDrift.ccList{ 1 }.dflcum;
cdrift2 = info.YdepDrift.ccList{ nStripes }.dflcum;


% PLOT  - cummulative drift for first and last stripe

figure
subplot(1,2,1)
title('Cummulative Drift')
hold on, plot( cdrift1(1,:) )
hold on, plot( cdrift2(1,:) )
xlabel('slice')
ylabel('y')
hold off
subplot(1,2,2)
hold on, plot( cdrift1(2,:) )
hold on, plot( cdrift2(2,:) )
hold off
xlabel('slice')
ylabel('x')
legend('Stripe 1',sprintf('Stripe %d',nStripes))




%% 8 . Y-DEP DRIFT . Spline fitting

% TASK  separate statistical image-to-image shifts from long range
%       cummulative drifts induced by changes of the image content as
%       function of z

% NOTE  To estimate the smoothing factor extract an array from the 
%       variable ccC.dflcum and save into a variable. Use the fitting tool
%       'cftool' and reduce the fitting parameter until the fit is smooth
%       across several image slice numbers


% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 8;
processName = "Y-DEP DRIFT . Spline fit";

% process Note
processNote = "";


% UPDATE -- PROCESS -----------------------------------------------------

% average smoothing parameter
smoothingParameter = 0.001;
 
% plot fit (true/false)
% ( CAREFUL: generates a plot for each stripe )
plotFit = true;


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Y-DEP DRIFT . Spline fitting ...\n', processStep );

% start timer
tic


% set up fit
ft   = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = smoothingParameter;

% copy cross correlation data
ccList = info.YdepDrift.ccList;

% number of stripes
nStripes = max(size( info.YdepDrift.ccList ));

for i=1:nStripes
    % for i=1:1

    if verbose == true
        PrintIndex(i,nStripes);
    end

    % current drift values
    ccC = ccList{i};
    
    % drift along image axes
    d1 = transpose( ccC.dflcum(1,:) );
    d2 = transpose( ccC.dflcum(2,:) );
    
    
    % fit axis 1
    [xData1, yData1] = prepareCurveData( [], d1 );
    [fitresult1, ~]  = fit( xData1, yData1, ft, opts );
    d1Fit            = feval(fitresult1,xData1 );
    
    % fit axis 2
    [xData2, yData2] = prepareCurveData( [], d2 );
    [fitresult2, ~]  = fit( xData2, yData2, ft, opts );
    d2Fit            = feval(fitresult2,xData2 );
    
    
    % remove long-range drift
    d1s = d1 - d1Fit;
    d2s = d2 - d2Fit;
    
    % save fit and image-to-image drift
    ccC.dflcumFit   = [ transpose(d1Fit); transpose(d2Fit) ];
    ccC.dflcumShort = [ transpose(d1s); transpose(d2s) ];
    

    % update cell array
    ccList{i} = ccC;
    
    if plotFit == true

        % PLOT

        figure
        subplot(2,2,1)
        hold on, plot(d1,'b+')
        hold on, plot(d1Fit,'r')
        subplot(2,2,2)
        hold on, plot(d2,'b+')
        hold on, plot(d2Fit,'r')
        subplot(2,2,3)
        hold on, plot(d1-d1Fit,'b+')
        subplot(2,2,4)
        hold on, plot(d2-d2Fit,'b+')

    end

end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

info.YdepDrift.nStripes           = nStripes;
info.YdepDrift.fitType            = 'SmoothingSpline';
info.YdepDrift.smoothingParameter = smoothingParameter;
info.YdepDrift.ccList             = ccList;


% SAVE . add process step

newProcess            = table;
newProcess.Session    = info.Session;
newProcess.Step       = processStep;
newProcess.Process    = processName;
newProcess.Save       = "";

newProcess.iStart     = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd       = info.ProcessList.iEnd( processStep-1 );
newProcess.nX         = info.ProcessList.nX( processStep-1 );
newProcess.nY         = info.ProcessList.nY( processStep-1 );
newProcess.nZ         = info.ProcessList.nZ( processStep-1 );

newProcess.Data       = "YdepDrift";
newProcess.Note       = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

info.ProcessList(processStep,:) = newProcess;


% info structure modified
info.Dates.Modified = timeFinished;

if autoSaveInfo == true

    % update time
    info.Dates.Saved = timeFinished;

    % save to file
    save( infoFile,'info','-v7.3' )

end


fprintf( '\n ... Finished successfully\n' );



%% ... PLOT . Short-range drift statistics

% TASK  plot mean and std of drift along x and y as function of y position

% NEED  data from previous section


% standard deviations
d1sStd = zeros(nStripes,1);
d2sStd = zeros(nStripes,1);

% mean deviation
d1sMean = zeros(nStripes,1);
d2sMean = zeros(nStripes,1);

% mean absolute deviations
d1sMeanAbs = zeros(nStripes,1);
d2sMeanAbs = zeros(nStripes,1);


for i=1:nStripes
    PrintIndex(i,nStripes);
    
    % current drift values
    ccC = ccList{i};

    % standard deviations
    d1sStd(i) = std( ccC.dflcumShort(1,:) );
    d2sStd(i) = std( ccC.dflcumShort(2,:) );
    
    % standard deviations
    d1sMean(i) = mean( ccC.dflcumShort(1,:) );
    d2sMean(i) = mean( ccC.dflcumShort(2,:) );
    
    % mean abs deviations
    d1sMeanAbs(i) = mean(abs( ccC.dflcumShort(1,:) ));
    d2sMeanAbs(i) = mean(abs( ccC.dflcumShort(2,:) ));

end


% PLOT - deviations as function of y position

% y positions
sCenterT = transpose( sCenter );

% add reference point (center of central marker triplet)
sCenterT = cat(1,-225,sCenterT);


d1sStd = cat(1,0,d1sStd);
d2sStd = cat(1,0,d2sStd);
d1sMean = cat(1,0,d1sMean);
d2sMean = cat(1,0,d2sMean);


figure
hold on, plot(sCenterT,d1sStd,'b+-')
hold on, plot(sCenterT,d1sMean,'b+')
hold on, plot(sCenterT,d2sStd,'ro-')
hold on, plot(sCenterT,d2sMean,'ro')
% hold on, plot(sCenterT,d1sMeanAbs,'b')
% hold on, plot(sCenterT,d2sMeanAbs,'r')
legend('y axis - std','y axis - mean','x axis - std', 'x axis - mean')
xlabel( 'slice center position' )
ylabel( 'displacement (pixels) ')



%%  9.  Y-DEP DRIFT . Interpolate sampled drift

% TASK  prepare displacement matrix with pixel resolution

% PROC  1. Collect drift values as function of central y-position for
%          stripes
%       2. Interpolate drift values to pixel resolution
%       3. Plot interpolated x,y displacements for analysis


% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 9;
processName = "Y-DEP DRIFT . Interpolate sampled drift";

% process Note
processNote = "";


% UPDATE -- PROCESS -----------------------------------------------------

interpolationMethod = 'pchip';
 
% -----------------------------------------------------------------------

fprintf( 'Process step %u: Y-DEP DRIFT . Interpolate sampled drift ...\n', processStep );

% start timer
tic


% PARAMETERS

% number of stripes
nStripes = max(size( info.YdepDrift.ccList ));

% number of images
iStart = info.ProcessList.iStart( processStep-1 );
iEnd   = info.ProcessList.iEnd( processStep-1 );
ni     = iEnd - iStart + 1;

% sampled displacement series
dSample = zeros(2,nStripes,ni);

% current series dimension
ny = info.ProcessList.nY( processStep-1 );

% sample points for interpolation
iY = 1:ny;

% interpolated drift (from dSample at sCenter points)
dInt = zeros(2,ny,ni);

% stripe centers
nStripes = info.YdepDrift.nStripes;
sHalf    = info.YdepDrift.sHalf;
s        = 2*sHalf + 1;

si      = 0 : nStripes-1;
sCenter = info.YdepDrift.sHalf+1 + s*si;


% COLLECT SAMPLED DRIFT

for i=1:nStripes

    if verbose == true
        PrintIndex(i,nStripes);
    end
    
    % current strip
    ccC = info.YdepDrift.ccList{i};
    
    % current displacement
    dc = ccC.dflcumShort;
    
    % add to series
    dSample(1,i,:) = dc(1,:);
    dSample(2,i,:) = dc(2,:);
    
end


% INTERPOLATE AXIS 1

% video name
vNameY = sprintf('S%u_P%u_DisplacementY.avi', info.Session,processStep );
vNameY = fullfile( info.Project.projectFolder,'Figures',vNameY );

% prepare video
v1 = VideoWriter( vNameY );
open(v1);
figure

for i=1:ni
    PrintIndex(i,ni);
    
    % sampled displacements of current image
    ds1i = squeeze(dSample(1,:,i));
    
    % interpolate for current image
    di1i = interp1(sCenter,ds1i,iY,interpolationMethod,'extrap');
    
    % save into array
    dInt(1,:,i) = transpose( di1i );
    
    % plot and save frame in video
    plot(sCenter,ds1i,'bo')
    hold on, plot(iY,di1i,'b')
    hold off
    
    ylim([-20,40])
    frame = getframe;
    writeVideo(v1,frame);
    
end

close(v1);


% INTERPOLATE AXIS 2

% video name
vNameX = sprintf('S%u_P%u_DisplacementX.avi', info.Session,processStep );
vNameX = fullfile( info.Project.projectFolder,'Figures',vNameX );

v2 = VideoWriter(  vNameX);
open(v2);
figure

for i=1:ni
    PrintIndex(i,ni);
    
    % sampled displacements of current image
    ds2i = squeeze(dSample(2,:,i));
    
    % interpolate for current image
    di2i = interp1(sCenter,ds2i,iY,interpolationMethod,'extrap');
    
    % save into array
    dInt(2,:,i) = transpose( di2i );
    
    % plot and save frame in video
    plot(sCenter,ds2i,'bo')
    hold on, plot(iY,di2i,'b')
    hold off
    
    ylim([-10,25])
    frame = getframe;
    writeVideo(v2,frame);
    
end

close(v2);

% process finished
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

info.AlignYdep.Session          = info.Session;
info.AlignYdep.processStep      = processStep;
info.AlignYdep.processAppliedTo = processStep - 1;
info.AlignYdep.Date             = timeFinished;
info.AlignYdep.Note             = processNote;

info.AlignYdep.sCenter  = sCenter;
info.AlignYdep.dSample  = dSample;
info.AlignYdep.interpolationMethod  = interpolationMethod;
info.AlignYdep.dInt     = dInt;
info.AlignYdep.iY       = iY;


% SAVE . add process step

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = "";

newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.nX       = info.ProcessList.nX( processStep-1 );
newProcess.nY       = info.ProcessList.nY( processStep-1 );
newProcess.nZ         = info.ProcessList.nZ( processStep-1 );
newProcess.Data     = "AlignYdep";

newProcess.Note     = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

info.ProcessList(processStep,:) = newProcess;


% info structure modified
info.Dates.Modified = timeFinished;

if autoSaveInfo == true

    % update time
    info.Dates.Saved = timeFinished;

    % save to file
    save( infoFile,'info','-v7.3' )

end


fprintf( '\n ... Finished successfully\n' );


%%  10 . Y-DEP DRIFT . Warp images  >> Ydd

% TASK  calculate displacement matrix and apply imwarp()

% NEED  all sections of Y-DEP DRIFT group

% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 10;
processName = "Y-DEP DRIFT . Warp images";

% process name acronym (3 characters, start with capital letter)
processSave = "Ydd";

% process Note ( appears in process list )
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders (true,false)
useFolders = true;

if useFolders == true

    projectFolderRead = info.Project.projectFolder;
    projectFolderSave = projectFolderRead;

end


% UPDATE -- PROCESS -----------------------------------------------------


% -----------------------------------------------------------------------

fprintf( 'Process step %u:Y-DEP DRIFT . Warp images ...\n', processStep );

% start timer
tic


% PARAMETERS

if useFolders == true

    % slice names
    sliceName = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate output folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end

end


% image dimensions of current series
iStart = info.ProcessList.iStart( processStep-1 );
iEnd   = info.ProcessList.iEnd( processStep-1 );
ny     = info.ProcessList.nY( processStep-1 );
nx     = info.ProcessList.nX( processStep-1 );

% 1D displacement components
dy = squeeze( info.AlignYdep.dInt(1,:,:) );
dx = squeeze( info.AlignYdep.dInt(2,:,:) );


if useFolders == false

    if exist( 'IMSP','var' )

        % reset image series variable from previous process
        IMS = IMSP;
        clearvars IMSP

    end

    % check consistency
    ni = iEnd - iStart + 1;
    [nyV,nxV,niV] = size( IMS );

    if nyV~=ny || nxV~=nx || niV~=ni
        error('At least on dimension of loaded image series does not fit to chosen parameters')
    end

    % corrected image series
    IMSD = zeros( size(IMS),'like',IMS );

end

% for i=1:5
for i=iStart:iEnd

    if verbose == true
        PrintIndex(i,iEnd);
    end
    
    % index in array
    k = i - iStart + 1;

    if useFolders == true

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

    else

        IT = IMS(:,:,k);

    end

    % 2D displacement field components
    DX = repmat(dx(:,k),1,nx);
    DY = repmat(dy(:,k),1,nx);
    
    % 2D displacement field
    D = zeros(ny,nx,2);
    D(:,:,1) = DX;
    D(:,:,2) = DY;
    
    % apply imwarp
    IW = imwarp(IT,-D,'cubic');


    if useFolders == true

        % save to file
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( IW,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')
    else

        IMSD(:,:,k) = IW;

    end
    
end


% stop timer
elapsedTime = toc;
% __ sec

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

info.AlignYdep.sCenter  = sCenter;
info.AlignYdep.dSample  = dSample;
info.AlignYdep.interpolationMethod  = interpolationMethod;
info.AlignYdep.dInt     = dInt;
info.AlignYdep.iY       = iY;


% SAVE . add process step

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.nX       = info.ProcessList.nX( processStep-1 );
newProcess.nY       = info.ProcessList.nY( processStep-1 );
newProcess.nZ       = info.ProcessList.nZ( processStep-1 );
newProcess.Data     = "AlignYdep";

newProcess.Note     = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

info.ProcessList(processStep,:) = newProcess;


% info structure modified
info.Dates.Modified = timeFinished;

if autoSaveInfo == true

    % update time
    info.Dates.Saved = timeFinished;

    % save to file
    save( infoFile,'info','-v7.3' )

end


fprintf( '\n ... Finished successfully\n' );


%% ... PLOT . Compare before and after

% TASK  Illustrate effect of drift correction on representative volume

% HERE  - Assume that folders are used
%       - use same image region that was used for denoising

% UPDATE -- PROCESS -----------------------------------------------------

%i n/out
loadFolder        = info.Project.projectFolder;
sliceName         = info.Project.sliceName;

sliceFolderBefore = 'ESB_ExtCroDen';
sliceFolderAfter  = 'ESB_ExtCroDenYdd';

% image region
xS = 1288;
xE = 1898;
yS = 1663;
yE = 2231;

% z range for test volume
zS = info.ProcessList.iStart( processStep ) + 305;
zE = info.ProcessList.iStart( processStep ) +406;

% yz-slice in collage for comparison
yzOrig = 378;
yzDC   = 986;

% -----------------------------------------------------------------------

% sub volume dimensions
nx = xE - xS + 1;
ny = yE - yS + 1;
nz = zE - zS + 1;

% initialize sub volume
IDCS = zeros(ny,2*nx,nz,info.Extract.ImageClass );

% initialize running index
k = 0;

for i=zS:zE
    PrintIndex( i,zE );

    % read images
    fname = fullfile( loadFolder,sliceFolderBefore,sliceName(i) );
    IORIG = imread( fname );

    fname = fullfile( loadFolder,sliceFolderAfter,sliceName(i) );
    IDC   = imread( fname );

    % extract regions
    IORIG = IORIG(yS:yE,xS:xE);
    IDC   = IDC(yS:yE,xS:xE);

    % concatenate
    ICOL = cat(2,IORIG,IDC);

    % running index
    k = k + 1;

    % add to subvolume
    IDCS(:,:,k) = ICOL;

end

% SAVE . sub volume
fImg = sprintf( 'S%u_P%u_DriftCorrection.mrc',...
    info.Session,processStep);
fImg = fullfile(info.Project.projectFolder,'Figures',fImg );
TomoSaveMRC( IDCS,fImg,6 );


% compare yz-slize through bacterium
ICOMP = cat( 2,squeeze( IDCS(:,yzOrig,:) ),squeeze( IDCS(:,yzDC,:) ) );

% PLOT . yz-slice comparison
plim( ICOMP )

% SAVE . comparison file
fComp = sprintf( 'S%u_P%u_DriftCorr_Compare_YZslice.tif',...
    info.Session,processStep);
fComp = fullfile(info.Project.projectFolder,'Figures',fComp );
SaveImage( ICOMP,fComp,'TIF' );


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% ADD MORE PROCESS SECTIONS HERE




