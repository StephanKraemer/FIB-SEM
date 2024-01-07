%
% NOTE
%
% This template was extracted from the corresponding Matlab script in the
% project demo. At the top it contains all sections that are necessary to 
% initiallize a second session as well as input/output sections. The last
% section illustrates how to use y-dependent drift values from the first
% session and apply it to the un-denoised data set.
%
% CAREFULLY ADJUST THE UPDATE SECTIONS TO YOUR GIVEN PROJECT AND ADD MORE
% SECTIONS AS NEEDED
% 


%% INFO . Load

% TASK  reload project processing data to continue the current or copy
%       data for a new session

% NEED  go to Align folder

% UPDATE -- WORKFLOW ----------------------------------------------------

% info structure of current session (true/false)
% (otherwise assumed to be a previous session that will be built on)
currentSession = false;

% session number
session = 1;


% PROCESS CONTROL PARAMETERES

% print image indeces during calculations (true/false)
verbose = true;

% automatically save info file at end of each process (true/false)
autoSaveInfo = true;


% -----------------------------------------------------------------------

% info file name for loading
iFile = sprintf( 'info-Session%d.mat',session );

% load info structure
ldata = load( iFile );


if currentSession == true

    % extract info structure
    info = ldata.info;

    % save loading time (to compare with saving time)
    info.Dates.Loaded = datetime( 'now' );

    % set global path for saving info file
    infoFile = fullfile( info.Project.projectFolder,'Align',iFile );

else

    % extract previous info structure
    infoPrev = ldata.info;

    % save loading time
    infoPrev.Dates.Loaded = datetime( 'now' );

end


% clean up
clearvars ldata


%% Session 2

% TASK  start a new session that builds on a previous workflow
%       1. generate info structure copying processes and corresponding
%          data from previous session(s)
%       2. define script behavior

% NEED  1. load the previous session using section 'INFO . Load' into
%          variable 'infoPrev' !
%       2. make sure there is no variable 'info' in workspace


if exist( 'info','var')

    % warn if info variable exists
    error( 'Workspace contains info structure. Make sure that all relevant data are saved and delete the variable')

end


% UPDATE -- INFO --------------------------------------------------------

% project identifier
info.Name = infoPrev.Name; 

% new session number
info.Session = 2;

% description
info.Description = "Apply alignment to original data";

% detector signal used in this script
info.Detector = infoPrev.Detector;

% last process number in previous session before new workflow
procLast = 3;


% UPDATE -- WORKFLOW ----------------------------------------------------

% print image indeces during calculations?
verbose = true;

% automatically save info file at end of each process (true/false)
autoSaveInfo = true;


% -----------------------------------------------------------------------

% initialize new dates
timeNow             = datetime( 'now' );
info.Dates.Created  = timeNow;
info.Dates.Modified = timeNow;
info.Dates.Loaded   = 0;
info.Dates.Saved    = 0;

% COPY . Run Notes ProcessList
info.Run         = infoPrev.Run;
info.Notes       = infoPrev.Notes;
info.ProcessList = infoPrev.ProcessList(1:procLast,:);

% copy data
for i=1:procLast
    PrintIndex( i,procLast );

    % current sub structure
    subStruct = infoPrev.ProcessList.Data(i);

    % copy structure
    info.(subStruct) = infoPrev.(subStruct);

end


% info file name for saving
infoFile = sprintf( 'info-Session%d.mat',info.Session );

% set global path for saving info file
infoFile = fullfile( info.Project.projectFolder,'Align',infoFile );



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% PROJECT NOTES

% TASK  Note particularities of the run and important results throughout
%       the analysis

% PROC  Comment out previous entry to omit multiple entries. Add new note.


% UPDATE ----------------------------------------------------------------

% TASK  ( 'add note','save to Excel', copy from here! )
noteTask = 'add note';

% process number for note ( if applicable, set 0 for general note )
noteProcess = 7;

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
saveToExcel = true;

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
session = 2;


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

% save to info structure
saveToInfo = false;

% type of load folder ( Project,Storage )
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
pStep = 10;


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

% dynamic range of TIFF image ('uint8,'uint16' ) 
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
%%  4 . Y-DEP DRIFT . Warp images  >> Ydd

% TASK  calculate displacement matrix and apply imwarp()

% NEED  all sections of Y-DEP DRIFT group

% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 4;
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

% transfer relevant data from previous session
info.YdepDrift = infoPrev.YdepDrift;
info.AlignYdep = infoPrev.AlignYdep;


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

info.AlignYdep.D  = D;


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


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% ADD MORE PROCESS SECTIONS HERE


