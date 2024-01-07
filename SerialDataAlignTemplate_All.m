%% ABOUT THIS SCRIPT

% Last modified:    2023-12-31

% This script performs image processing and alignment tasks on serial
% section images in preparation for data segmentation. It handles images
% recorded on
% 
%   - a Helios FIB (via Slice & View )
%   - a Crossbeam FIB (via Atlas tomography widget)
% 
% The goal of this script is to provide a flexible and well documented
% workflow for the data analysis depending on the processing needs of the
% given serial section dataset. It affords some patience during the
% setup but will hopefully help to keep track of the applied squence of 
% image processing tasks and relevant parameters.
% 
% Version 1.0   Most commonly needed processes, systematically structured
%               process sections
% 
% 
% For questions, bug reports and general comments please contact
% 
%       Stephan Kraemer (skraemer@fas.harvard.edu)
%       Center for Nanoscale Systems, Harvard University



%% INTRODUCTION 

% SCRIPT NAVIGATION AND EXECUTION
% 
% The script makes use of a series of _Matlab sections_ (denoted by the
% double-percent symbol '%%' in the  text). It can be best viewed in the
% Matlab Editor and easily navigated with the 'Go To' button in the
% Editor control tab.
% 
% Every section represents one process step which follows strictly from the
% previous ones. The individual section can be executed either by clicking 
% on the blue ribbon on the left edge of the text editor that appears when
% the cursor is placed inside of it or via the 'Run Section' command in
% the SECTION panel of the Editor tab.
% 
% Each section has the same structure containing a text block that briefly
% describes the given process, an UPDATE section that allows the user to
% define the values of parameters relevant for the given process and a 
% program body which then executes the process, adds a line to a process
% list and saves relevant results into an _info_ structure.
% 
% Most sections use their input from the _info_ structure that contains
% relevant results from previous sections and the parameters in the UPDATE
% section. In rare cases a section relies on local variables that were
% generated in the immediately preceding section(s). Those are marked with
% three dots '...' at the end of the section title.
% 
% A new section can be built using the TEMPLATE at the end this template
% file. The template represents a scaffold that handles the basic
% organization of the script. Modify the code following the instructions 1
% to 10 and add the process-specific code at the marked location.
% 
% Simple tasks are directly programmed in the script, programming-intensive
% steps are executed via _Matlab functions_ which are needed to run the
% script. For time intensive calculations (such as denoising) the functions
% are compiled into MEX files.
% 
% At this point the functions do not make explicit use of GPUs to speed up
% image processing tasks. Where possible, parfor loops are used to
% distribute the work over multiple Matlab clones.
% 
% 
% BUILDING A SCRIPT FROM SCRATCH
% 
% 1. Build a script by copying Matlab sections from this template file
%    in the sequence as needed by the given data set.( Each scripts
%    represents a processing _session_ with ONE given series of processes. 
%    Build separate scripts for alternate sequences.)
% 
% 2. Change the process step number in the title of the section and update
%    the corresponding variable in the UPDATE-WORKFLOW section ( THIS STEP
%    IS IMPORTANT TO INSURE THE PROPER FLOW OF THE SCRIPT !).
% 
% 3. Modify all the values for the variables defined in the UPDATE section
%    as needed by the given problem.
% 
% Two example scripts are provided. Each contain a typical set of sections:
% reading and saving of the info structure, reading and saving of images,
% initialize the session and info structure, parsing the image folders, 
% denoising the images and aligning the images. The two scripts distinguish
% in their alignment process. The first uses a basic rigid body 
% displacement via cross correlation of entire images (see GLOBAL DRIFT
% section group in this template). The second one can be used in case 
% images show y-dependent scan distortions. In this case the images are
% split into horizontal stripes. The drift values obtained from cross
% correlation in each stripe is used to calculate a y-dependent displacment
% field that is used by the imwarp() function to correct for those
% distortions (see Y-DEP DRIFT section group).
% 
% 
% STARTING A SESSION
% 
% A new script is started using the 'SESSION X' section. Each session
% represents one sequence of process steps. Replace X with 1, or any higher 
% number if different process sequences are tested. Session 2 could be also
% used if multiple signals are recorded during the serial section run and
% the processes developed in session 1 are applied to the images of the
% other signal in session 2.
% 
% The 'SESSION X' script allows to define two global workflow parameters.
% The variable 'verbose' defines if each image number is shown as each
% section is executing the process. The variable 'autoSaveInfo' defines
% whether the info structure is autmatically saved after the execution of
% each process. It is recommended to keep it turned on since each process 
% depends on the results of the previous one and a sesion is to be
% continued on another day.
% 
% 
% The following sections initialize the info structure:
% 
%   The 'RUN DATA' setion allows to store parameters particular to the
%   execution of the given serial section run, such as microscope type,
%   beam parameters, ect.
% 
%   The 'PROJECT FOLDER' defines the location of a given project on the
%   computer or computer system. The script can distinguish between
%   various locations. See more explanation below. For saving and reading
%   data _absolute_ paths are used.
% 
%   The 'PROJECT NOTES' allows to save information that is pertinent to the
%   properties of the given data set or its analysis. The user can add a
%   new line to the  table at any point during the execution of the script.
% 
%   The 'PROJECT LIST' section initializes the table info.ProcessList. At
%   the end of each section that changes the dataset or generates
%   information that is used later in the script, a new line is added to
%   the ProcessList table. This table documents the sequence of processes
%   applied to the dataset and stores the dimenions of the dataset! THIS
%   TABLE ALLOWS THE USER TO SEE IN A GLANCE WHAT PROCESSES HAVE BEEN
%   APPLIED AND HOW THE IMAGE DIMENSIONS ARE CHANGED IN EACH STEP.
% 
% 
% CONTINUE A SESSION
% 
% The analysis of a dataset can be resumed at any time. Just use the 'INFO
% . Load' section to load the info structure back into memory. 
% 
% How to continue next depends on how serial section data are handled. 
% If the entire dataset is loaded into memory, do so with the 'IMAGES . 
% Load from folder' section. If each section handles the loading of
% individual images, one can immediately continue with the pending
% processing step. See next info block.
% 
% As in the 'SESSION X' section, it is possible to define the two global
% workflow paramters.
% 
% 
% STORAGE AND HANDLING OF SERIAL SECTION DATA
% 
% Serial section data are currently solely stored as a series of TIF files
% in a folder. Two sections ('IMAGES . Load from folder' and 'IMAGES . 
% Save to folder') handle loading and saving of entire data cubes. 
% 
% For inspection purposes, it is possible to load a subset of images by
% skipping a number of slices. To further save memory, images can be
% reduced to 8-bit.
% 
% Image data can be handled in two different ways, depending on the amount
% of memory available on the computer. In the first approach images are
% successively loaded from the hard drive at each processing step, and
% resulting images are stored back onto the hard drive after each
% task. In order to speed up the process one can set up a work folder on
% a SSD hard drive. For the second approach the entire data set is loaded
% into memory. The latter case is intended for speedy evaluation. In that
% case, the intermediate result HAS to be saved via section 'IMAGES
% . Save to folder' if the analysis is stopped and continued at another
% time.
% 
% There is currently no memory management. If the dataset is loaded into 
% memory, the memory has to be large enough to hold the given number of 
% Matlab clones, two representations of the data, before and after
% processing, plus overhead for the actual calculations.
% 
% 
% STORAGE OF INTERMEDIATE RESULTS
% 
% The image series can be saved after significant process steps allowing to
% reanalyze the effect of a given processing task. Each series is stored
% in a folder which starts with the acronym of the currently used signal, 
% e.g. 'ESB','SE2', 'InLens' for Crossbeam data, 'TLD','ICD','CBS' for 
% Helios data. Each process step is described by a three-letter acronym.
% The acronym is saved in the ProcessList. The folder name contains the
% sequence of applied processing steps.
% 
% As example, the folder 'ESB_ExtDenCroYdd' contains an image series that
% has undergone the following sequence of processing steps
% 
%   1. Ext: extract the 'XROI' from a Crossbeam dataset
%   2. Den: denoise via non-local means averaging
%   3. Cro: crop series to region of interest
%   4. Ydd: correct image drift, taking into account y-dependent drift
%      values in case images suffer from drift-induced local scan
%      distortions.
% 
% Each section that saves an intermediate result is marked in the title by 
% '>> <ACRONYM>'
% 
% 
% STORAGE OF PROCESSING PARAMETERS
% 
% Processing parameters and results are saved in the _info_ structure,
% which can be stored and reloaded via section 'INFO . Load','INFO . Save'. 
% Each section saves parameters into a corresponding sub structures. The
% info structure is the core element of this script that stores and
% analysis results and contains parameters determining the workflow of this
% script.
% 
% Whenever a process (apart from intermediate steps) is finished, a line
% in the _ProcessList_ table is added, documenting the processing steps so
% far applied to the dataset. The table also contains information about
% the data volume size at each step, the name of sub structures in the info
% structure which contains input paramters and results for the given
% process step. Notes can be added and the duration of each process is
% stored.
% 
% 
% STRUCTURAL ELEMENTS OF A SECTION
% 
% Each section represents a specific processing task (such as cropping
% the images series, denoising, etc.). Every script is built from template
% sections and is tailored to the processing needs of the specific data
% series. The analysis process is a linear workflow where each section
% follows _strictly_ from the previous one.
% 
% - Sections that start with a _number_ represent a task that extracts
%   necessary information or modifies the image series. The numbering
%   allows to guide the reader through the sequence of analysis steps.
% 
%  - Sections marked with 'PLOT' visualize a certain aspect of an
%    individual process step.
% 
% - Sections marked with '...' denote intermediate processing steps that
%   calculate parameters needed in a following section.
% 
% 
% All sections have the same layout of structural elements:
% 
% - The 'TASK' line briefly describes the process.
% 
% - A 'NEED' line indicates if the analysis needs specific variables,
%   for example the image series, if it is loaded into memory. In general
%   sections were designed to be as autonomous as possible where most
%   variables can be retrieved from the _info_ structure and _ProcessList_
%   table (see above). This allows to interrup the analysis and resume at
%   another time. It also allows to apply the same script to a second
%   signal with minimal changes, which can be useful for example in cases
%   where both secondary and backscatter images were recorded during the
%   run.
% 
% - A 'PROC' line describes the intended handling of the sections if more
%   specific instructions are necessary.
% 
% - An 'UPDATE' text block contains the input parameters needed for the
%   given task. It is distinguished betwen WORKLFOW, INPUT/OUTPUT and
%   PROCESS parameters. The user only needs to interact with this part of
%   the section. All following code should run without intervention.
% 
% - At the end of the section the _info_ structure and _ProcessList_ table
%   are updated with input paramters, results of the process and potential
%   changes to the size of the data set.
% 
% - The sections use tic/toc command pairs to monitor the time spent for
%   each processing step which is saved into the ProcessList table.
% 
% - A 'Dates' substructure keeps track of when the script was created,
%   modified, laded and saved
% 
% 
% As mentioned above, a user can easily generate a new process using the 
% 'TEMPLATE' section, located at the end of this template file. 
% It contains the main building blocks. Code that needs to be adjusted is
% marked with instructions % TODO X:...'.



%% CONVENTIONS
 
% COORDINATE SYSTEM
% 
% Axes are named following the conventions of the Image Processing
% Toolkit (and most other image processing software):
%           x: horizontal cross section image axis
%           y: vertical cross section image axis
%           z: slice direction
%  
% NOTE: The _indeces_ of the corresponding 3D data cube correspondingly
%       represent the y,x, and z axis, e.g.
%           [ny,nx,nz] = size( IMS );
% 
% SLICE INDEXING
% 
% The script reads and saves the filenames of a data folder in the
% 'sliceName' string array. The initial set of files contains 
% 'nSlicesRecorded' slices, with indeces 
%           iStart = 1;
%           iEnd   = nSlicesRecorded;
% 
% NOTE: The above slice indexing only follows the order of files, given
%       by the formatting of their names. The numbers can differ from the
%       ones in the filenames if those start with a number different from 1.



%% PROJECT ORGANIZATION

% PROJECT FOLDERS
%
% All data related to one serial section run are stored in a _Project_ 
% folder. The script allows to distribute the data over several platforms. 
% Their absolute locations are stored in the info.Project structure and
% each section provides the option to set the desired read and save
% locations. Currently three locations are considered:
%
%   1. external _storage_ drive: used to store the raw data including
%      information about the run provided by the respective instrument. It
%      is also intended for longterm storage once the analysis of the serial
%      section run is finished to free up space on the local active drive.
%
%   2. local _project_ drive: contains currently active projects. It stores
%      all intermediate results. Once the analysis of a run is finished the
%      user can decide whether to delete intermediate results or move them
%      to the storage drive to free up space.
%
%   3. _work_ drive: a drive that allows fast file reading and saving such
%      as a SSD drive. The entire project folder can be moved/copied there
%      temporarily for the time-intensive part of the analysis. This
%      can be useful if the local memory is small and images are
%      continuously read from the hard drive (See section 'INTRODUCTION').
%
% Examples
%
%   info.Project.storageFolder = 'F:\User\Stephan\Serial\Projects\2022-12-09 Copper single particle I'; 
%   info.Project.projectFolder = 'D:\User\Stephan\Serial\Projects\2022-12-09 Copper singel particle I';
%   info.Project.workFolder    = 'C:\WORK\Current Project';


% PROJECT FOLDER STRUCTURE
%
% The project folder is divided into a few sub folders:
%
% 'Run'             all supplementary data that are produced during the 
%                   preparation and execution of the run, such as log 
%                   files, Atlas 3D setup files, etc. Typically stored on
%                   storage drive.
%
% 'Align'           Matlab script(s) 'SerialDataEval_Session1.m', etc
%                   Corresponding 'info' structure stored as *.mat file
%                   'info-Session1.mat', etc. See below definition of a
%                   _session_.
%
% 'Figures'         any images and plots that describe analysis results
%
% 'Segmentation'    further Matlab scripts that contain
%                   segmentation-related processes such as morphological
%                   operations, maybe at some point Deep Learning based
%                   segmentation
% 
% data folders      as described in the 'INTRODUCTION' section,
%                   intermediate serial section data, marked by the name
%                   of the detector and applied processes, symbolized by
%                   3-letter acronyms. For example
%                   
%                   ESB_Ext, ESB_ExtCroDen, ESB_ExtCroDen, ...
%
%                   The above naming convention allows to infer the applied
%                   process sequence on the operating system level without
%                   having to open the info file.

% A SESSION
%
% A session contains one specific sequence of processes applied to the
% image series. Once a session is complete, it is possible to make use of
% the resulting data in the 'info' structure to build further sequences:
%
% One potential scenario is to use the same sequence of processes on
% another signal, if a multi-signal run was performed. In this case, run
% a second session with a copy of the first script and change the input
% data from sigal one to two. In order to apply the same image alignments
% that were used on the first data set, skip all pre-processing steps and
% immediately apply the beam shift values from the first session.
%
% Another scenario is to evaluate alternate processing sequencies, for
% example to test the influence of denoising on segmentation results. 
% In this case it is possible to build a starting info structure using
% section 'COPY INFO FROM PREVIOUS SESSION' and start a new script with the
% following processes that build on the copied ones.



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
sliceFolder = 'ESB_ExtDen';

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
sliceFolderSave = 'ESB_ExtDenYddInv_8bit';


% UPDATE -- PROCESS -----------------------------------------------------

% slice increment when skipping images
nInc = 1;

% dynamic range of TIFF image ('uint8,'uint16' )
dynRange = 'uint16';

% image format ( 'TIFF','JPG' )
iFormat = 'TIFF';

% image quality of JPG
jpgQuality = 80;


% -----------------------------------------------------------------------


% start and end slice
iStart = info.ProcessList.iStart( pStep );
iEnd   = info.ProcessList.iEnd( pStep );

% slice names
sliceName = info.Project.sliceName;

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



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% SESSION 1

% TASK  Set up info structure and start first session

% NEED  - go to 'Align' folder in project folder
%       - make sure that all data from previous sessions are save 
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
info.Name = "Veera Panova . Run 2022-12-09 . Copper Single Particle I";

% current session number
info.Session = 1;

% description
info.Description = "Align images of TLD series";

% Detector signal used in this session
info.Detector = "TLD";


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

info.Run.Instrument         = "Helios";
info.Run.Software           = "Slice&View 4";
info.Run.Date               = datetime("2022-12-09");

info.Run.Project            = "Veera Panova";
info.Run.Name               = "Copper Single Particle I";

info.Run.Units.Lenght       = "nm";
info.Run.Units.Energy       = "kV";
info.Run.Units.Current      = "pA";

info.Run.pixelSize          = 13;
info.Run.sliceThickness     = 50;

info.Run.Detectors          = ["TLD","ICD"];   % multiple: string array

info.Run.Ebeam.Voltage      = 3;
info.Run.Ebeam.Current      = 1600;
info.Run.Ibeam.Voltage      = 30;
info.Run.Ibeam.Current      = 800;

info.Run.Projection         = "Corrected";


% number of digits of slice number in file string
nDigits                     = 3;

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
info.Project.workFolder    = 'C:\WORK\Stephan\Current Project';
info.Project.projectFolder = 'D:\User\Stephan\Serial\Project . Veera Panova\2022-12-09 (Veera) Copper singel particle I';
info.Project.storageFolder = 'F:\User\Stephan\Serial\Projects . Veera Panova\2022-12-09 (Veera) Copper single particle I'; 

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
noteProcess = 0;

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

        % save to Excel spread sheet

        % filename
        fXls = sprintf( 'ProcessList_Session%u.xlsx', info.Session );
        fXls = fullfile( info.Project.projectFolder,'Align',fXls );

        % save to Excel
        writetable( info.ProcessList,fXls );

    end
end


        
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Each sequence of processes is saved as a _session_. All relevant data are
% stored in the corresponding _info_ structure. The _ProcessList_ table 
% provides the size of the dataset after each process.
% 
% Use the section below in case an alternate sequence of processes is to be
% used that builds on part of the processes applied in the previous
% session.


%% Session X

% TASK  start a new session that builds on a previous workflow
%       1. generate info structure copying processes and corresponding
%          data from previous session(s)
%       2. define script behavior

% NEED  1. load the previous session using section 'INFO . Load' into
%          variable 'infoPrev' !
%       2. make sure there is no variable 'info' in workspace


if exist( 'info','var')

    % check if info variable exists
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
%%  0 . PREPARE . Crossbeam . Parse simple file list

% TASK  extract filenames and get image dimensions from first image

% NOTE  Assumes that - all images have the same size
%                    - name of the folder corresponds to info.Detector

% UPDATE -- WORKFLOW ----------------------------------------------------

processStep = 0;
processName = "PREPARE . Crossbeam . Parse simple file list";

processNote = "";


% UPDATE -- PROCESS -----------------------------------------------------

% data location
projectFolderRead = info.Project.storageFolder;


% -----------------------------------------------------------------------

fprintf( 'ProcessStep %u: Crossbeam . Parse simple file list ...\n',processStep );

% start timer
tic


% FILE NAMES AND NUMBER OF SLICES

% input folder
sliceFolder = convertStringsToChars( info.Detector );

% files in slice directory
dList = dir(fullfile( projectFolderRead,sliceFolder ));

if isempty( dList )
    error('Files not found. Please check if path is correct.')
end

% exclude folders '.', '..' , ...
isDir = cat( 1,dList.isdir );
dList( isDir==1 ) = [];

% number of slices
nSlicesRecorded = max(size( dList ));

% file names sliceName = strings(nSlicesRecorded,1);
sliceName = strings( nSlicesRecorded,1 );
for i=1:nSlicesRecorded
    sliceName(i) = convertCharsToStrings( dList(i).name );
end


% IMAGE INFORMATION

% read first image of series
fname  = fullfile( projectFolderRead,sliceFolder,sliceName(1) );
IT     = imread( fname );

% image dimensions
[ny,nx] = size( IT );


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

info.Project.nSlicesRecorded    = nSlicesRecorded;
info.Project.sliceName          = sliceName;

% beginning start and end slices
info.Project.iStart = 1;
info.Project.iEnd   = nSlicesRecorded;
info.Project.iRef   = round( 0.5*(info.Project.iStart + info.Project.iEnd) );


% SAVE  . generate first entry in process list table

info.ProcessList.Session    = info.Session;
info.ProcessList.Step       = processStep;
info.ProcessList.Process    = processName;
info.ProcessList.Save       = "";

info.ProcessList.iStart     = 1;
info.ProcessList.iEnd       = nSlicesRecorded;
info.ProcessList.nX         = nx;
info.ProcessList.nY         = ny;
info.ProcessList.nZ         = nSlicesRecorded;

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

fprintf( '\n ... Finished successfully\n' );



%%  0 . PREPARE . Crossbeam . Parse Atlas file list

% TASK  . extract filenams and slice thicknesses from filenames

% UPDATE -- WORKFLOW ----------------------------------------------------

processStep = 0;
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



%%  0 . PREPARE . Crossbeam . Extract XROI  >> Ext

% TASK  Extract largest rectangle in images

% NOTE  . in case images contain also tracking markers
%       . assuming only one XROI

% PROC  algorithm finds horizontal and vertical edges of patches, chooses
%       rectangle with largest edge length


% UPDATE -- WORKFLOW ----------------------------------------------------

processStep = 0;
processName = "PREPARE . Crossbeam . Extract XROI";

% process name acronym (3 characters, start with capital letter)
processSave = "Ext";

% process note
processNote = "";


% plot figures (true/false)
plotResults = true;


% UPDATE -- IN/OUT ------------------------------------------------------

% location of raw data
projectFolderRead = info.Project.storageFolder;

% use folders
useFolders = true;

if useFolders == true

    % output
    projectFolderSave = info.Project.projectFolder;

end


% UPDATE -- PROCES ------------------------------------------------------

% extracted area: envelope of all images or smallest common area
area = 'Envelope';
% area = 'Intersection';


% -----------------------------------------------------------------------

fprintf( 'ProcessStep %u: Crossbeam . Extract XROI ...\n',processStep );

% start timer
tic


% slice folder names
[sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
    info.Detector,processStep,info.ProcessList,processSave );

% slice names
sliceName = info.Project.sliceName;

% generate save folder
[SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
if SUCCESS == 0
    error( '%s (ID %s)',MESSAGE,MESSAGEID )
end

% parameters
pPrev   = processStep - 1;
iStart  = info.ProcessList.iStart( pPrev );
iEnd    = info.ProcessList.iEnd( pPrev );
nSlices = iEnd - iStart + 1;

% boundbox
BoundingBox = zeros( nSlices,4 );

% image class list
ImgClass = string;


% FIND BOUNDING BOX

parfor i=iStart:iEnd

    if verbose == true
        PrintIndex( i,iEnd );
    end

    % read current slice image
    fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
    IT    = imread( fname );

    % image class (bit depth)
    ImgClass(i) = convertCharsToStrings(class( IT ));

    % non-zero pixels
    BW = int8( IT>0 );

    % internal edges of regions
    E  = uint16( ...
        ( BW - circshift(BW,[0,1])  )>0 | ...
        ( BW - circshift(BW,[0,-1]) )>0 | ...
        ( BW - circshift(BW,[1,0])  )>0 | ...
        ( BW - circshift(BW,[-1,0]) )>0 );

    % length of horizontal and vertical lines
    lh = sum(E,2);
    lv = sum(E,1);

    % figure
    % subplot(1,2,1)
    % plot(lh,'b+-')
    % subplot(1,2,2)
    % plot(lv,'b+-')

    % find largest rectangle from two longest edges in both directions
    [~,ihs] = sort( lh,"descend" );
    [~,ivs] = sort( lv,"descend" );
    ihs      = ihs(1:2);
    ivs      = ivs(1:2);

    % origin (pixel diagonal in front of rectangle) [ox,oy]
    ox = min( ivs ) - 1;
    oy = min( ihs ) - 1;

    % size of rectangle
    nx = max( ivs ) - ox;
    ny = max( ihs ) - oy;

    % save bounding box
    BoundingBox(i,:) = [ ox,oy,nx,ny ];

end


% IMAGE CLASS

if     all( ImgClass == "uint8" ),  imgClass = 'uint8';
elseif all( ImgClass == "uint16" ), imgClass = 'uint16';
else,  error( 'Folder contains images with varying dynamic range' )
end


% GLOBAL BOUNDING BOX

switch area
    case 'Envelope'
        oxg = min( BoundingBox(:,1) );
        oyg = min( BoundingBox(:,2) );

        exmax = max( BoundingBox(:,1) + BoundingBox(:,3) );
        eymax = max( BoundingBox(:,2) + BoundingBox(:,4) );

        nxg = exmax - oxg;
        nyg = eymax - oyg;

    otherwise
        oxg = max( BoundingBox(:,1) );
        oyg = max( BoundingBox(:,2) );

        exmin = min( BoundingBox(:,1) + BoundingBox(:,3) );
        eymin = min( BoundingBox(:,2) + BoundingBox(:,4) );

        nxg = exmin - oxg;
        nyg = eymin - oyg;

end

BoundingBoxGlobal = [ oxg,oyg,nxg,nyg ];


% EXTRACT DESIRED AREA

if useFolders == true

    parfor i=iStart:iEnd

        if verbose == true
            PrintIndex( i,iEnd );
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT    = imread( fname );

        % extract XROI
        ITE = IT( oyg+1:oyg+nyg,oxg+1:oxg+nxg );

        % save to file
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( ITE,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')
    end

else

    % allocate memory for image series
    IMS = zeros( nyg,nxg,nslices,imgClass );

    parfor i=1:nSlices

        if verbose == true
            PrintIndex( i,nSlices );
        end

        % image index
        imgIdx = iStart + k - 1;

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(imgIdx) ); %#ok<PFBNS> 
        IT    = imread( fname );

        % extract XROI
        ITE = IT( oyg+1:oyg+nyg,oxg+1:oxg+nxg );

        % save to variable
        IMS(:,:,i) = ITE;

    end    

end


% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

info.Extract.Session            = info.Session;
info.Extract.processStep        = processStep;
info.Extract.processAppliedTo   = pPrev;
info.Extract.Date               = timeFinished;   
info.Extract.Note               = processNote;

info.Extract.BoundingBox        = BoundingBox;
info.Extract.BoundingBoxGlobal  = BoundingBoxGlobal;
info.Extract.GlobalArea         = area;
info.Extract.ImageClass         = imgClass;


% SAVE  . update process list

newProcess            = table;
newProcess.Session    = info.Session;
newProcess.Step       = processStep;
newProcess.Process    = processName;
newProcess.Save       = processSave;

newProcess.iStart     = info.Project.iStart;
newProcess.iEnd       = info.Project.iEnd;
newProcess.nX         = nxg;
newProcess.nY         = nyg;
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


% PLOT . Bounding box

if plotResults == true

    figure

    subplot(2,2,1)
    hold on, plot( BoundingBox(:,1))
    xlabel( 'slice' )
    ylabel( 'ox' )
    title( ' Global Bounding Box')

    subplot(2,2,2)
    hold on, plot( BoundingBox(:,2))
    xlabel( 'slice' )
    ylabel( 'oy' )

    subplot(2,2,3)
    hold on, plot( BoundingBox(:,3))
    xlabel( 'slice' )
    ylabel( 'nx' )

    subplot(2,2,4)
    hold on, plot( BoundingBox(:,4))
    hold off
    xlabel( 'slice' )
    ylabel( 'ny' )

end

fprintf( '\n ... Finished successfully\n' );



%%  0 . PREPARE . Crossbeam . Extract Manually  >> Ext

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
processStep = 0;
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
ox = 2664;
oy = 4572;
wx = 3144 + 4993 - ox;
wy = 4584 + 4364 - oy;


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



%%  0 . PREPARE . Helios . Parse file list

% TASK  extract filenames and image information, for now:
%           <ScanArea>:<Width>,<Height>
%           <PixelSize>:<X>,<Y>

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number / name
processStep = 0;
processName = "PREPARE . Helios . Parse file list";

processNote = "";


% UPDATE -- INPUT/OUTPUT ------------------------------------------------

% collect input data
projectFolderRead = info.Project.storageFolder;


% UPDATE -- PROCESS -----------------------------------------------------

% pixel units
Image.PixelUnits = 'm';

% meta key strings to find pixel size
pixStart = '<PixelSize>';
pixEnd   = '</PixelSize>';
xStart   = '<X unit="m" unitPrefixPower="1">';
xEnd     = '</X>';
yStart   = '<Y unit="m" unitPrefixPower="1">';
yEnd     = '</Y>';


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Helios . Parse file list ...\n', processStep );

% start timer
tic


% set workflow parameters
Image.ProcessStep = processStep;


% FILE NAMES AND NUMBER OF SLICES

% input folder
sliceFolder = convertStringsToChars( info.Detector );

% files in slice directory
dList = dir(fullfile( projectFolderRead,sliceFolder ));

if isempty( dList )
    error('Files not found. Please check if path is correct.')
end

% exclude folders '.', '..' , ...
isDir = cat( 1,dList.isdir );
dList( isDir==1 ) = [];

% number of slices
nSlicesRecorded = max(size( dList ));

% file names sliceName = strings(nSlicesRecorded,1);
sliceName        = strings( nSlicesRecorded,1 );
for i=1:nSlicesRecorded
    sliceName(i) = convertCharsToStrings( dList(i).name );
end


% IMAGE INFORMATION

% initialize image dimensions (x,y) and pixel size
Image.Width      = zeros( nSlicesRecorded,1 );
Image.Height     = zeros( nSlicesRecorded,1 );
Image.PixelSizeX = zeros( nSlicesRecorded,1 );
Image.PixelSizeY = zeros( nSlicesRecorded,1 );
Image.BitDepth   = zeros( nSlicesRecorded,1 );

for i=1:nSlicesRecorded
    % for i=100:100

    if verbose == true
        PrintIndex( i,nSlicesRecorded );
    end

    % read image information
    fname  = fullfile( projectFolderRead,sliceFolder,sliceName(i) );
    imInfo = imfinfo( fname );


    % image width and height
    Image.Width(i)  = imInfo.Width;
    Image.Height(i) = imInfo.Height;

    % bit depth
    Image.BitDepth(i) = imInfo.BitDepth;


    % xml string of meta data (somewhat obscure)
    xmlString = imInfo.UnknownTags(3).Value;

    % extract pixel size
    kpb = strfind( xmlString, pixStart ) + length( pixStart );
    kpe = strfind( xmlString, pixEnd ) - 1;
    pixString = xmlString(kpb:kpe);

    kxb = strfind( pixString, xStart ) + length( xStart );
    kxe = strfind( pixString, xEnd ) - 1;
    Image.PixelSizeX(i) = str2double( pixString(kxb:kxe) );

    kyb = strfind( pixString, yStart ) + length( yStart );
    kye = strfind( pixString, yEnd ) - 1;
    Image.PixelSizeY(i) = str2double( pixString(kyb:kye) );

end


% check consistency of image bit depth

if ~all( Image.BitDepth == 8 ) &&  ~all( Image.BitDepth == 16 )

    error( 'Folder contains images with different dynamic range' )

end


% end timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% RANGE OF DIMENSIONS

Image.Range.Width      = MiMa( Image.Width );
Image.Range.Height     = MiMa( Image.Height );
Image.Range.PixelSizeX = MiMa( Image.PixelSizeX );
Image.Range.PixelSizeY = MiMa( Image.PixelSizeY );


% SAVE . update info structure

info.Project.Session          = info.Session; 
info.Project.processStep      = processStep;
info.Project.processAppliedTo = 0;
info.Project.Date             = timeFinished;

info.Project.nSlicesRecorded  = nSlicesRecorded;
info.Project.sliceFolder      = sliceFolder;
info.Project.sliceName        = sliceName;


info.Image.Session          = info.Session;
info.Image.processAppliedTo = 0;
info.Image.Date             = timeFinished;

info.Image.Width            = Image.Width;
info.Image.Height           = Image.Height;
info.Image.PixelSizeX       = Image.PixelSizeX;
info.Image.PixelSizeY       = Image.PixelSizeY;
info.Image.BitDepth         = Image.BitDepth;
info.Image.Range            = Image.Range;


% SAVE . update process list

newProcess            = table;
newProcess.Session    = info.Session;
newProcess.Step       = processStep;
newProcess.Process    = processName;
newProcess.Save       = "";

newProcess.iStart     = 1;
newProcess.iEnd       = nSlicesRecorded;
newProcess.nX         = 0;
newProcess.nY         = 0;
newProcess.nZ         = nSlicesRecorded;

newProcess.Data       = "Project, Image";
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

figure
subplot(2,2,1)
plot(Image.Width)
ylabel('Width')
subplot(2,2,2)
plot(Image.Height)
ylabel('Height')
subplot(2,2,3)
plot(Image.PixelSizeX)
xlabel('Slice')
ylabel('Pixel Size X')
subplot(2,2,4)
plot(Image.PixelSizeY)
xlabel('Slice')
ylabel('Pixel Size Y')

fprintf( '\n ... Finished successfully\n' );



%%  0 . PREPARE . Helios . Scale Images >> Scl

% TASK  resample images, correcting for varying pixel sizes induced by
%       adjusting the selected area location and projection effect along y

% UPDATE -- WORKFLOW ------------------------------------------------------

% process number and name
processStep = 0;
processName = "PREPARE . Helios . Scale Images";

% process name acronym (3 characters, start with capital letter)
processSave = "Scl";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% input/output
projectFolderRead = info.Project.storageFolder;
projectFolderSave = info.Project.projectFolder;


% UPDATE -- PROCESS -----------------------------------------------------

% scale from projection to actual pixel size on cross section ?
% ( true: if not corrected during recording of the data )
correctProjection = false;

% interpolation method
method = 'bicubic';

% new pixel size
% ( set pNew to PixelSizeX in case projection has to be coorected along y )
pNew = max( info.Image.PixelSizeX );


% -----------------------------------------------------------------------

fprintf( 'ProcessStep %u: Helios . Scale Images ...\n',processStep );

% start timer
tic


% PARAMETERS

% input/output folders
[sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
    info.Detector,processStep,info.ProcessList,processSave );

% generate folder
[SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
if SUCCESS == 0
    error( '%s (ID %s)',MESSAGE,MESSAGEID )
end

% slice name
sliceName = info.Project.sliceName;

% scale factors for image width based on provided pixel size
% ( correct changes during run )
scaleX = info.Image.PixelSizeX ./ pNew;
scaleY = info.Image.PixelSizeY ./ pNew;

if correctProjection == true

    % additionally correct for projection along y
    scaleY = scaleY ./ sind( 52 );
end


% current image dimensions
imageWidth  = info.Image.Width;
imageHeight = info.Image.Height;

% new image dimensions (x,y) and pixel size
nSlices               = info.Project.nSlicesRecorded;
imageScaledWidth      = zeros( nSlices,1 );
imageScaledHeight     = zeros( nSlices,1 );
imageScaledPixelSizeX = repmat( pNew,[nSlices,1] );
imageScaledPixelSizeY = repmat( pNew,[nSlices,1] );


parfor i=1:nSlicesRecorded
    % for i=596:597

    if verbose == true
        PrintIndex( i,nSlicesRecorded );
    end

    % read image
    fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
    IT = imread( fname );

    % image clase
    imgClass = class( IT ); %#ok<PFTUSW> 


    % scale image
    nx = round( scaleX(i) * imageWidth(i) );
    ny = round( scaleY(i) * imageHeight(i) );
    IS = imresize( double( IT ), [ny nx], method);

    % save dimensions to structure
    imageScaledWidth(i)  = nx;
    imageScaledHeight(i) = ny;


    % save scaled image
    fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );

    switch imgClass
        case 'uint8'
            imwrite( uint8( IS ),fname,'tiff', ...
                'Compression','packbits','WriteMode','overwrite' )
        case 'uint16'
            imwrite( uint16( IS ),fname,'tiff', ...
                'Compression','packbits','WriteMode','overwrite' )
    end

end


% end timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

info.ImageScaled.Session          = info.Session;
info.ImageScaled.processStep      = processStep;
info.ImageScaled.processAppliedTo = processStep - 1;
info.ImageScaled.Date             = timeFinished;
info.ImageScaled.Note             = processNote;

info.ImageScaled.pNew             = pNew;
info.ImageScaled.scaleX           = scaleX;
info.ImageScaled.scaleY           = scaleY;

info.ImageScaled.Width            = imageScaledWidth;
info.ImageScaled.Height           = imageScaledHeight;
info.ImageScaled.PixelSizeX       = imageScaledPixelSizeX;
info.ImageScaled.PixelSizeY       = imageScaledPixelSizeY;

info.ImageScaled.Range.Width      = MiMa( imageScaledWidth );
info.ImageScaled.Range.Height     = MiMa( imageScaledHeight );
info.ImageScaled.Range.PixelSizeX = MiMa( imageScaledPixelSizeX );
info.ImageScaled.Range.PixelSizeY = MiMa( imageScaledPixelSizeY );

info.ImageScaled.BitDepth         = info.Image.BitDepth;


% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

% update dimensions
newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.nX       = info.ProcessList.nX( processStep-1 );
newProcess.nY       = info.ProcessList.nY( processStep-1 );
newProcess.nZ       = info.ProcessList.nZ( processStep-1 );

newProcess.Data     = "ImageScaled";
newProcess.Note     = processNote;     
newProcess.Date     = timeFinished;
newProcess.Duration = elapsedTime;

% info.ProcessList(pStep,:) = struct2table( newProcess );
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

figure
subplot(2,2,1)
plot(info.ImageScaled.Width)
ylabel('Width')
subplot(2,2,2)
plot(info.ImageScaled.Height)
ylabel('Height')
subplot(2,2,3)
plot(info.ImageScaled.PixelSizeX)
xlabel('Slice')
ylabel('Pixel Size X')
subplot(2,2,4)
plot(info.ImageScaled.PixelSizeY)
xlabel('Slice')
ylabel('Pixel Size Y')


fprintf( '\n ... Finished successfully\n' );



%%  0 . PREPARE . Helios . Extract images  >> Ext

% TASK  Equalize potential image size differences and save into datacube
%       with fixed width/height

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number and name
processStep = 0;
processName = "PREPARE . Helios . Extract images";

% process name acronym (3 characters, start with capital letter)
processSave = "Ext";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% input
projectFolderRead = info.Project.projectFolder;

% use folders for output
useFolders = true;

if useFolders == true

    projectFolderSave = info.Project.projectFolder;

end


% UPDATE . PROCESS ------------------------------------------------------

% desired size of area ('smallest','largest')
area = 'smallest';


% centering of smaller image in larger or cropping of larger image to fit
% into smaller frame

% x: 'left','center','right'
% y: 'top', 'center','bottom'

centerX = 'center';
centerY = 'bottom';

% extract as image class 'uint8' or 'uint16'
% ( see info.Image.BitDepth for dynamic range of original images )
imageClass = 'uint16';


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Helios . Extract images\n',processStep );

% start timer
tic


% slice folder names
[sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
    info.Detector,processStep,info.ProcessList,processSave );

% generate folder
[SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
if SUCCESS == 0
    error( '%s (ID %s)',MESSAGE,MESSAGEID )
end

% slice name
sliceName = info.Project.sliceName;


% start and end slices
iStart = 1;
iEnd   = info.Project.nSlicesRecorded;


% dynamic range scale factors
s8  = double(intmax( 'uint8' ));
s16 = double(intmax( 'uint16' ));


% image sizes
if isfield( info,'ImageScaled' )
    imageInfo = info.ImageScaled;
else
    imageInfo = info.Image;
end


switch area
    case 'smallest'
        % smallest common area
        nx = min( imageInfo.Width );
        ny = min( imageInfo.Height );

    otherwise
        % largest common area
        nx = max( imageInfo.Width );
        ny = max( imageInfo.Height );
end


if useFolders ~= true

    % allocate image series memory
    IMS = zeros( ny,nx,info.Project.nSlicesRecorded,imageClass );

end


parfor i=iStart:iEnd

    if verbose == true
        PrintIndex( i,iEnd );
    end

    % default positions
    ox = 0;
    oy = 0;


    % load image
    fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
    IT = imread( fname );

    % current image size
    [nyt,nxt] = size( IT );


    if ~isa( IT,imageClass )

        % adjust dyanmice range
        switch imageClass
            case 'uint8'
                IT = uint8( s8/s16 .* double( IT ));

            case 'uint16'
                IT = uint16( s16/s8 .* double( IT ));

        end
    end


    % initialize extracted image frame
    IE = zeros( ny,nx,imageClass );


    if strcmpi( area,'smallest' )

        if nxt == nx
            % current image has smallest width
            ox = 0;

        elseif nxt > nx
            % crop to smallest width
            switch centerX
                case 'left'
                    ox = 0;

                case 'center'
                    ox = floor( 0.5*(nxt - nx ));

                case 'right'
                    ox = nxt - nx;

            end
        end

        if nyt == ny
            % current image has smallest height
            oy = 0;
        
        elseif nyt > ny
            % crop to smallest height
            switch centerY
                case 'left'
                    oy = 0;

                case 'center'
                    oy = floor( 0.5*(nyt - ny ));

                case 'right'
                    oy = nyt - ny;

            end
        end

        % shrink larger areas to smallest common area
        IE  = IT(oy+1:oy+ny,ox+1:ox+nx);

    end


    if strcmpi( area,'largest' )

        if nxt == nx
            % current image has largest width
            ox = 0;
        
        elseif nxt < nx
            % placement in larger frame
            switch centerX
                case 'left'
                    ox = 0;

                case 'center'
                    ox = floor( 0.5*(nx - nxt ));

                case 'right'
                    ox = nx - nxt;

            end
        end


        if nyt == ny
            % current image has largest height
            oy = 0;
        
        elseif nyt < ny
            % placement in larger frame
            switch centerY
                case 'left'
                    oy = 0;

                case 'center'
                    oy = floor( 0.5*(ny - nyt ));

                case 'right'
                    oy = ny - nyt;

            end
        end

        % place image in larger frame
        IE(oy+1:oy+ny,ox+1:ox+nx) = IT;

    end
    

    if useFolders == true

        % save image
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( IE ,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    else

        IMS(:,:,i) = IE;

    end
 
end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE .  update  info structure

info.Extract.Session           = info.Session;
info.Extract.processStep       = processStep;
info.Extract.processAppliedTo  = processStep - 1;
info.Extract.Date              = timeFinished;
info.Extract.Note              = processNote;

info.Extract.projectFolderRead = projectFolderRead;
info.Extract.sliceFolderRead   = sliceFolderRead;

if useFolders == true
info.Extract.projectFolderSave = projectFolderSave;
info.Extract.sliceFolderSave   = sliceFolderSave;
end

info.Extract.sliceName         = sliceName;
info.Extract.iStart            = iStart;
info.Extract.iEnd              = iEnd;
info.Extract.area              = area;
info.Extract.centerX           = centerX;
info.Extract.centerY           = centerY;
info.Extract.ImageClass        = imageClass;


% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = "Ext";

% update dimensions
newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.nX       = nx;
newProcess.nY       = ny;
newProcess.nZ       = info.ProcessList.nZ( processStep-1 );

newProcess.Data     = "Extract";
newProcess.Note     = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

% info.ProcessList(pStep,:) = struct2table( newProcess );
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
%%  0 . HISTOGRAM . Evaluate series

% TASK  Calculate histogram for full series

% NOTE  It is possible to
%       - calculate histogram for any given processing stage for which
%         data were saved. Either load the image series into 'IMS' 
%         variable or provide 'sliceFolderRead' link.
%       - skip saving the results, if only data exploration is the goal


% UPDATE -- WORKFLOW ----------------------------------------------------

% process step / name
processStep = 0;
processName = "HISTOGRAM . Histogram series";

% process note
processNote = "";


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



%% 0 . HISTOGRAM . Global Stretch  >>  Hst

% TASK  stretch histogram in case experimental dynamic range narrower than
%       the one given by image class

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number and name
processStep = 0;
processName = "HISTOGRAM . Global Stretch";

% process name acronym ( 3 characters, start with capital letter )
processSave = "Hst";


% process note ( appears in ProcessList table )
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders
useFolders = true;

if useFolders == true

    % input, output
    projectFolderRead = info.Project.projectFolder;
    projectFolderSave = projectFolderRead;

end


% UPDATE -- PROCESS -----------------------------------------------------

% min/max intensity range
iMin = 10000;
iMax = 60000;


% -----------------------------------------------------------------------

fprintf( 'ProcessStep %u: HISTOGRAM . Global Stretch ...\n',processStep );

% start timer
tic


% <PARAMETERS>

% dynamic range of image class

switch info.Extract.ImageClass

    case 'uint8'
        scl  = double(intmax( 'uint8' ));

    case 'uint16'
        scl = double(intmax( 'uint16' ));

    otherwise
        error( 'Image class not implemented' )
end



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

    iStart = info.ProcessList.iStart(processStep-1);
    iEnd   = info.ProcessList.iEnd(processStep-1);


    % NOTE: THE FOLLOWING FOR LOOP CAN BE REPLACED BY A PARFOR LOOP IF 
    %       VARIABLES FOLLOW THE NECESSARY RULES ( SEE MATLAB DOCU )

    for i=iStart:iEnd

        if verbose == true
            PrintIndex(i,iEnd);
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % remember image class (uint8 or uint16)
        % NOTE: IN CASE IMAGES ARE TRANSFORMED TO CLASS DOUBLE
        % TO ENSURE PRECISION OF CALCULATIONS
        imgClass = class( IT );


        % crop intensities
        IT( IT<iMin ) = iMin;
        IT( IT>iMax ) = iMax;

        IP = scl/(iMax-iMin) * ( double(IT) - iMin );
        
       
        switch imgClass
            case 'uint8'
                IP = uint8( IP );

            case 'uint16'
                IP =  uint16( IP );

            otherwise
                error('image class not implemented')
        end


        % save
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( IP,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    end

else

    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    % allocate memory for resulting image series of this setion
    IMSP = zeros(size(IMS),'like',IMS);

    % remember image class
    imgClass = class( IMS );

    % number of images
    [~,~,ni] = size( IMS );


    for i=1:ni

        if verbose == true
            PrintIndex( i,ni );
        end

        % read current image
        % 
        IT = IMS(:,:,i);

        % crop intensities
        IT( IT<iMin ) = iMin;
        IT( IT>iMax ) = iMax;

        IP = scl/(iMax-iMin) * ( double(IT) - iMin );


        % save
        switch imgClass
            case 'uint8'
                IMSP(:,:,i) = uint8( IP );

            case 'uint16'
                 IMSP(:,:,i) = uint16( IP );

            otherwise
                error('image class not implemented')
        end

    end
end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

nInstance = 0;
if isfield( info,'HistStretch' )

    % number of process events so far
    nInstance = max(size( info.HistStretch ));

end

info.HistStretch( nInstance+1 ).Session          = info.Session;
info.HistStretch( nInstance+1 ).processStep      = processStep;
info.HistStretch( nInstance+1 ).processAppliedTo = processStep - 1;
info.HistStretch( nInstance+1 ).Date             = timeFinished;
info.HistStretch( nInstance+1 ).Note             = processNote;

if useFolders == true
    info.HistStretch( nInstance+1 ).sliceFolderRead     = sliceFolderRead;
    info.HistStretch( nInstance+1 ).sliceFolderSave     = sliceFolderSave;
end

info.HistStretch( nInstance+1 ).iMin = iMin;
info.HistStretch( nInstance+1 ).iMax = iMax;
    

% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;


iProcPrev = processStep-1;
newProcess.iStart   = info.ProcessList.iStart( iProcPrev );
newProcess.iEnd     = info.ProcessList.iEnd( iProcPrev );
newProcess.nX       = info.ProcessList.nX( iProcPrev );
newProcess.nY       = info.ProcessList.nY( iProcPrev );
newProcess.nZ       = info.ProcessList.nZ( iProcPrev );

newProcess.Data     = "HistStretch";
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



%% ... HISTOGRAM . Fit peaks ...

% TASK  fit n Gaussian peaks to histogram to evaluate intensity drift 
%       ( and dynamic range ? )

% NOTE  Selecting index range in info.Histogram(iHist).HS:
%          - If several peaks are clearly separated throughout the series 
%            concentrate the range that contains a peak in the histogram.
%          - If the histogram contains one broad intensity distribution
%            that might be affected by global intensity changes try to use
%            a broad intensity range covering the entire intenity band.


% UPDATE -- Workflow ----------------------------------------------------

% UPDATE -- Process -----------------------------------------------------

% histogram number (default 1, other if histograms at different processing
% steps were calculated)
iHist = 1;

% fit range in histogram image (see NOTE)
% (index values, intensities calculated using ival)
idxStart = 1;
idxEnd   = 2048;

% number of Gaussian peaks
% !!! so far only 1 !!!
nPeaks = 1;

% function: g = a1* exp( -( (x-b1)/c1^2 ))


% test images for plotting
iTest = [ 1,260,548,554 ];

% save movie (true,false)
saveMovie = true;


% -----------------------------------------------------------------------

fprintf( 'Process step %u: Fit histogram peaks ...\n', processStep );


% cropped histogram
HS = info.Histogram(iHist).HS(:,idxStart:idxEnd);

% corresponding intensity bins
ival = info.Histogram(iHist).ival(idxStart:idxEnd);

% corresponding fit
HF = zeros( size(HS) );

% number of images
[nSlices,~] = size( HS );


% maximal intensity
iMax = max( ival );

% set up fit
ftGauss = sprintf( 'gauss%u',nPeaks );
ft      = fittype( ftGauss );

% fit options
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';


% saved peak postion and peak width
pPosition = zeros( nSlices,1 );
pWidth    = zeros( nSlices,1 );


for i=1:nSlices
    PrintIndex( i,nSlices );

    % prepare input variables
    [x,y] = prepareCurveData( ival,HS(i,:) );

    % set starting values (from histogram peak)
    % (peak height, position, width)
    [hMax,mPos] = max(y);
    opts.StartPoint = [ hMax,ival(mPos),0.25*iMax ];

    % fit Gaussian
    [hFit,~] = fit( x,y,ft,opts );

    % fit curve
    yf   = feval( hFit,x );

    % peak position
    coeff = coeffvalues( hFit );
    pPosition(i) = coeff(2);
    pWidth(i)    = coeff(3);


    % save fit
    HF(i,:) = yf';

end


% PLOT

% peak positions and widths
figure
subplot(1,2,1)
plot( pPosition )
xlabel('slice')
ylabel('peak intensity')
subplot(1,2,2)
plot( pWidth )
xlabel('slice')
ylabel('peak width')


% profile for test images
nTest = max(size( iTest ));

figure
for i=1:nTest
    hold on, plot( ival,HS(iTest(i),:),'b+')
    hold on, plot( ival,HF(iTest(i),:),'r')
end
hold off


if saveMovie == true

    % movie file
    fname = sprintf( 'HistogramFit.avi' );
    fname = fullfile( info.Project.projectFolder,'Figures',fname );


    % start video
    v1 = VideoWriter( fname );
    open(v1);
    figure

    for i=1:nSlices
        PrintIndex(i,nSlices);

        % plot data
        plot( ival,HS(i,:),'b+')
        hold on, plot( ival,HF(i,:),'r')
        hold off

        % add frame to video
        frame = getframe;
        writeVideo(v1,frame);

    end

    % stop video
    close(v1);

end

fprintf( '\n ... Finished successfully\n' );



%%  0 . HISTOGRAM . Scale  >> Hsc

% TASK  Scale intensities based on fitted peak positions in histogram
%       ( optional: scale based on peak position or fitted width )

% NEED  results from previous section

% UPDATE -- WORKFLOW ----------------------------------------------------

% process step / name
processStep = 0;
processName = "HISTOGRAM . Scale";

% process name acronym (3 characters, start with capital letter)
processSave = "Hsc";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders
useFolders = true;

if useFolders == true

    % input/output
    projectFolderRead      = info.Project.projectFolder;
    projectFolderSave      = info.Project.projectFolder;
    
end


% UPDATE -- PROCESS -----------------------------------------------------

% mean peak intensity ( manual value, mean(pPosition) )
% pMean = 0;
pMean = mean( pPosition );

% scale mode
% 'none'  do not scale intensities (just shift)
% 'peak'  scale with peak position (relative to lower bound 'iLow')
% 'width' use pWidth to scale relative to reference width
scaleMode = 'peak';

switch scaleMode

    case 'peak'
        % lower limit of histogram
        iLow = info.Histogram(iHist).ival(1);

    case 'width'
        % reference width of intensity distribution 
        % (manual value, width at a particular slice, ...)
        refWidth = 6000;

end


% -----------------------------------------------------------------------

fprintf( 'Process Step %u: Scale histogram ...\n', processStep );

% start timer
tic


if useFolders == true

    % slice names
    sliceName = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end

    iStart = info.ProcessList.iStart( processStep-1 );
    iEnd   = info.ProcessList.iEnd( processStep-1 );


    for i=iStart:iEnd
        % for i=2:2

        if verbose == true
            PrintIndex(i,iEnd);
        end

        % array index
        k = i - iStart + 1;

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT    = imread( fname );

        % remember image class
        imgClass = class( IT );


        switch scaleMode

            case 'peak'

                % scale factor
                sP = (pMean - iLow) / (pPosition(k) - iLow);

                % scale image
                IS = iLow + sP*( double(IT) - iLow );

            case 'width'

                % scale factor
                sW = refWidth / pWidth(k);

                % center on pMean and scale by measured peak width (???)
                IS = pMean + sW*( double(IT) - pPosition(k) );


            otherwise

                % just shift
                IS = double(IT) + pMean - pPosition(k);


        end


        % save
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );

        switch imgClass
            case 'uint8'
                imwrite( uint8(IS),fname,'tiff', ...
                    'Compression','packbits','WriteMode','overwrite')

            case 'uint16'
                imwrite( uint16(IS),fname,'tiff', ...
                    'Compression','packbits','WriteMode','overwrite')

            otherwise
                error('image class not implemented')
        end

    end


else

    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    % allocate memory for resulting image series of this setion
    IMSP = zeros(size(IMS),'like',IMS);

    % remember image class
    imgClass = class( IMS );

    % number of images
    [~,~,ni] = size( IMS );

    % consistency check
    if( ni ~= nSlices )

        error( 'IMS dimensions do not fit to given number of slices of histogram' );

    end


    for i=1:ni

        if verbose == true
            PrintIndex( i,ni );
        end


        switch scaleMode

            case 'peak'

                % scale factor
                sP = (pMean - iLow) / (pPosition(i) - iLow);

                % scale image
                IS = iLow + sP*( double(IMS(:,:,i)) - iLow );

            case 'width'

                % scale factor
                sW = refWidth / pWidth(i);

                % center on pMean and scale by measured peak width (???)
                IS = pMean + sW*( double(IMS(:,:,i)) - pMean );


            otherwise

                % just shift
                IS = double(IMS(:,:,i)) + pMean - pPosition(i);


        end

        % save
        switch imgClass
            case 'uint8'
                IMSP(:,:,i) = uint8( IP );

            case 'uint16'
                IMSP(:,:,i) = uint16( IP );

            otherwise
                error('image class not implemented')
        end

    end

end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

nInstance = 0;
if isfield( info,'HistScale' )

    % number of process events so far
    nInstance = max(size( info.HistScale ));

end

info.HistScale( nInstance+1 ).Session             = info.Session;
info.HistScale( nInstance+1 ).processStep         = processStep;
info.HistScale( nInstance+1 ).processAppliedTo    = processStep-1;
info.HistScale( nInstance+1 ).Updated             = timeFinished;
info.HistScale( nInstance+1 ).processNote         = processNote;

info.HistScale( nInstance+1 ).type                = 'Histogram.Peak';

if useFolders == true
    info.HistScale( nInstance+1 ).sliceFolderRead = 'SE2_ExtCro';
    info.HistScale( nInstance+1 ).sliceFolderSave = 'SE2_ExtCroHeq';
end

info.HistScale( nInstance+1 ).iHist               = iHist;
info.HistScale( nInstance+1 ).idxStart            = idxStart;
info.HistScale( nInstance+1 ).idxEnd              = idxEnd;
info.HistScale( nInstance+1 ).iTest               = iTest;
info.HistScale( nInstance+1 ).pPosition           = pPosition;
info.HistScale( nInstance+1 ).pWidth              = pWidth;

info.HistScale( nInstance+1 ).pMean               = pMean;
info.HistScale( nInstance+1 ).ScaleMode           = scaleMode;
switch scaleMode
    case 'peak'
        info.HistScale( nInstance+1 ).iLow        = iLow;
        info.HistScale( nInstance+1 ).refWidth    = 0;
    case 'width'
        info.HistScale( nInstance+1 ).iLow        = 0;
        info.HistScale( nInstance+1 ).refWidth    = refWidth;
    otherwise
        info.HistScale( nInstance+1 ).iLow        = 0;
        info.HistScale( nInstance+1 ).refWidth    = 0;        
end


% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

iProcPrev = processStep-1;
newProcess.iStart   = info.ProcessList.iStart( iProcPrev );
newProcess.iEnd     = info.ProcessList.iEnd( iProcPrev );
newProcess.nX       = info.ProcessList.nX( iProcPrev );
newProcess.nY       = info.ProcessList.nY( iProcPrev );
newProcess.nZ       = info.ProcessList.nZ( iProcPrev );

newProcess.Data     = "HistScale";
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
%%  0 . DENOISE - Estimate noise level ...

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
processStep = 0;
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
idx = [100,200,300,400];

% uncomment if above indeces are given in cropped coordinates: 
% transform back to indices of raw dataset
% idx = idx + info.ProcessList.iStart( processStep-1 ) - 1;

% fit degree of polynomial to remove long-range intensity variation
degree = 4;

if subProcess == 2

    % (index 1: test image, index 2: region with near-constant intensity)
    p(1,1,:) = [927 1442 80 80] ;
    p(2,1,:) = [310 1113 42 189];
%     p(2,2,:) = [1130 1778 71 59];
%     p(2,3,:) = [1227 1975 99 90];
    p(3,1,:) = [1519 867 58 76];
    p(4,1,:) = [1027 1106 33 35];

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
p = [5442 329 1158 876];

% non-local means parameter space
% noise level, probe window, averaging window half-width
nlmNoise = [2000,4000];
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


% SAVE .  IDS AND IRS to mrc files

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

% table of image indeces and corresponding parameters
nlmTable = array2table( nlmShow, ...
    'VariableNames',{'Image','Noise','Probe','Avg'});
info.Denoise.Test.ParameterList = nlmTable;


if autoSaveInfo == true

    % update time
    info.Dates.Saved = datetime( 'now' );

    % save to file
    save( infoFile,'info','-v7.3' )

end



%%  0 . DENOISE . All images  >> Den

% TASK  Remove noise (using non-local means, Buades et. al.)

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number / name
processStep = 0;
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
roiStd = 430;

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

% TASK show degree of denoising

% UPDATE ----------------------------------------------------------------

%in/out
loadFolder        = info.Project.storageFolder;
sliceName         = info.Project.sliceName;

sliceFolderBefore = 'ESB_Ext';
sliceFolderAfter  = 'ESB_ExtDen';

% test slice
iTest = 312;


% -----------------------------------------------------------------------

% read images
fname = fullfile( loadFolder,sliceFolderBefore,sliceName(iTest) );
IB = double( imread( fname ));

fname = fullfile( loadFolder,sliceFolderAfter,sliceName(iTest) );
IA = double( imread( fname ));

% difference
ID = IB - IA;

plim( IB )
% plim( IA )
% plim( ID )



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% GLOBAL DRIFT - Including separation of jumps

% Use this set of sections in case images do not suffer from scan
% distortions and a simple image-to-image drift correction is sufficient
%
% It allows to locate isolated larger jumps that are treated as a separate
% entity. This step can be skipped if the cummulative drift contains only a
% slowly varying continuous function with additive noise.

% The analysis contains the following steps. Key variables are noted
%
% 1. Cross correlation
%    The function 'SerialCrossCorrSeries' uses a Fourier-based cross
%    correlation function to determine the image-to-image drift. A
%    polynomial fif of the cross correlation peak allows sub-pixel
%    resolution.
%
%    Two sub-pixel results are provided. The first one does a general
%    polynomial fit and interpolates the peak on a user-defined grid. A
%    second fit uses a polynomial function that contains the peak center as
%    a parameter. The first fit is used to calculate starting values for
%    the fit of the latter function. The second fit allows to estimate the
%    error in finding the peak position. Due to the inherent smoothness of
%    the cross correlation function the error is typically better than a
%    tenth of a pixel.
%
%    The function also provides the cummulative drift based on values of
%    either fit 1 or fit 2, depending on the user's choice. The
%    corresponding variables are 
%       dfl1        (sub-pixel image-to-image drift on fixed grid)
%       dfl2        (sub-pixel image-to-image drift parameterized)
%       dflcum      (cummulative sum of either fit 1 or fit 2)
%
% 2. Separate single jumps from continuous drift
%    Plotting a histogram of the image-to-image drifts (for example dfl2)
%    leads to a Gaussian distribution resulting from the statistical nature
%    of the drifts occuring during recording of the series. Jumps show up
%    as outliers.
%
%    The user provides a threshold value to separate drifts that are
%    considered jumps versus noise given as multiples of the standard
%    deviation of a Normal distribution, e.g. 3sigma or 5sigma. Jumps
%    should only show up sporadically. The resulting variables are stored
%    in the substructure 'info.ContDrift'
%       d1Step, d2Step     singular image jumps (1:y-axis, 2:x-axis)
%       d1Cont, d2Cont     continuous drift part with drift noise
%
% 3. Integrate along image series
%    The image-to-image drifts are combined via cummulative sum. Individual
%    jumps produce step functions, the continous part a slowly undulating
%    function. The long-range undulations are a consequence of the
%    sensitivity of the cross correlation to spatial changes of the overall 
%    intensity distribution.
%       d1StepCum, d1ContCum
%       d2StepCum, d2ContCum
%       d1Comb, d2Comb          (sum of both cummulative drifts)
%
% 4. Remove undulations of continuous cummulative drift
%    A spline fit is used to separate the slowly varying drift component
%    described above from the actual image-to-image drift which is more
%    statistical in nature. The fit result are stored under
%       info.SplineFit.d1Fit, info.SplineFit.d2Fit
%
%    NOTE 1: This correction ultimately assumes that there is no systematic
%            shift along the z direction, i.e. that the sample surface is
%            perfectly flat and perpendicular to the ion beam. A deviation
%            of this condition will lead to a distortion of the data cube.
%            Nevertheless, such a systematic variation, if measureable by
%            other means, can be easily added to the shift correction
%            introducing a seperate variable.
%    NOTE 2: There is a degree of arbitrariness in the choice of the
%            smoothing factor. ...
%
% 5. Correct shift via circshift
%    The ultimate drift correction is calculated by subtracting the spline
%    fit from the total cummulative drift. In this implementation a basic
%    circshift is used for the shift. This approach leaves the intensities
%    unaltered due to the pixel-wise shift. The bounds of the corresponding
%    precision of the drift corretion is +-0.5 pixels which appears to be
%    perfectly adequate.



%%  0 . GLOBAL DRIFT . Cross correlation

% TASK  measure image-to-image drift

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number / name
processStep = 0;
processName = "GLOBAL DRIFT . Cross correlation";

% process note
processNote = "";


% UPDATE -- Input/Output ------------------------------------------------

% use folders (true,false)
useFolders = true;

if useFolders == true

    projectFolderRead = info.Project.projectFolder;

end


% UPDATE -- Process -----------------------------------------------------

% CONTROL PARAMETERS

% Flags
ccIn.YNGauss           = 'yes';  % smooth with Gaussian
ccIn.YNgrad            = 'yes';  % use gradient image for cc    
ccIn.YNpartial         = 'no';  % use reduced area for cc (not here)
ccIn.YNleastCommonArea = 'yes';  % use only common non-0 area for cc    *1
ccIn.YNfit             = 'yes';  % apply polynomial fit to cc peak      *2

% *1 in case images were shifted during the run
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
ccIn.idxRef = 332;

% image for which peak fit is plotted
% ( value > iStart, using the value 0 suppresses the plot )
% ( here too many calls, set to 0 )
ccIn.idxPlot = 0;


% NOT YET IMPLEMENTED
% % half-width search window for peak finding
% % concentrate on region near center of cc, in case there are spurious peaks
% % further away due to imperfect image conditions
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


% -----------------------------------------------------------------------

fprintf( 'Process step %u: GLOBAL DRIFT . Cross correlation ...\n', processStep );

% start timer
tic

% update input structure

if useFolders == true

    % slice names
    sliceName = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,~ ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

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


if useFolders == true

    % calculate cross correlations
    ccOut = SerialCrossCorrSeries( ccIn );

else

    if exist( 'IMSP','var' )

        % reset image series variable from previous process
        IMS = IMSP;
        clearvars IMSP

    end

    ccOut = SerialCrossCorrSeries( ccIn,IMS );

end

% stop timer
elapsedTime = toc;
% 1025 sec (40 cores)

% current time
timeFinished = datetime( 'now' );


% SAVE .  update info structure

info.GlobalDrift.Session          = info.Session;
info.GlobalDrift.processStep      = processStep;
info.GlobalDrift.processAppliedTo = processStep - 1;
info.GlobalDrift.Date             = timeFinished;
info.GlobalDrift.Note             = processNote;

info.GlobalDrift.Parameters  = ccIn;
info.GlobalDrift.Results     = ccOut;


% SAVE . update process list

iProc                 = processStep-1;

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

newProcess.Data       = "GlobalDrift";
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

figure
subplot(1,2,1)
plot(ccOut.dflcum(1,:),'b+-')
xlabel ('slice')
ylabel('Y drift')
subplot(1,2,2)
plot(ccOut.dflcum(2,:),'b+-')
xlabel ('slice')
ylabel('X drift')


fprintf( '\n ... Finished successfully\n' );



%%  0 . GLOBAL DRIFT . Separate isolated jumps

% TASK  separate image-too-image jumps from artificial continuous
%       long-range drift induced by cross correlation

% NEED  results from previous session, stored in info structure

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number / name
processStep = 0;
processName = "GLOBAL DRIFT . Separate isolated jumps";

% process note
processNote = "";


% UPDATE -- Input/Output ------------------------------------------------

% UPDATE -- Process -----------------------------------------------------

% image-to-image drift fit used (1 or 2, 2 preferred)
useFit = 2;

% standard deviation scale factor (default 3)
sclStd = 3.0;

% histogram parameters
hpar.nbin   = 100;
hpar.mima   = [-200,200];
hpar.YNplot = 'yes';
fitIn.ng    = 1;


% -----------------------------------------------------------------------

fprintf( 'Process step %u: GLOBAL DRIFT . Separate isolated jumps ...\n', processStep );

% start timer
tic


% ESTIMATE NOISE OF IMAGE-TO-IMAGE DRIFT FROM HISTOGRAM

% image drift values used
if useFit == 1
    d1 = info.GlobalDrift.Results.dfl1(1,:);
    d2 = info.GlobalDrift.Results.dfl1(2,:);

else
    d1 = info.GlobalDrift.Results.dfl2(1,:);
    d2 = info.GlobalDrift.Results.dfl2(2,:);

end


% calculate histogarm
[ival1,nval1,fitOut1] = TomoHist(d1,hpar,fitIn);
[ival2,nval2,fitOut2] = TomoHist(d2,hpar,fitIn);

% standard deviation from Gaussian fit
d1Sig = fitOut1.coeff(3)/sqrt(2);
d2Sig = fitOut2.coeff(3)/sqrt(2);


% FIND OUTLIERS

% outliers
iOut1 = find( abs(d1) > sclStd*d1Sig );
iOut2 = find( abs(d2) > sclStd*d2Sig );

% fill array with outlier values
d1Step = zeros(size(d1));
d2Step = zeros(size(d2));

d1Step(iOut1) = d1(iOut1);
d2Step(iOut2) = d2(iOut2);


% fill 2nd array with all other drift values
d1Cont = d1;
d2Cont = d2;

% remove outliers
d1Cont(iOut1) = 0;
d2Cont(iOut2) = 0;

% cummulative drift
d1ContCum = cumsum(d1Cont);
d1StepCum = cumsum(d1Step);
d2ContCum = cumsum(d2Cont);
d2StepCum = cumsum(d2Step);

% combined cummulative drift
d1Comb = d1ContCum + d1StepCum;
d2Comb = d2ContCum + d2StepCum;



% PLOT . image-to-image drift separated, and corresponding cummulative
%        drift

% plot mutually excluded
nSlices = max(size( d1 ));
iSlice  = 1:nSlices;

iSliceCont1          = iSlice;
iSliceCont1( iOut1 ) = [];
d1ContP              = d1Cont;
d1ContP( iOut1 )     = [];
iSliceStep1          = iSlice( iOut1 );
d1StepP              = d1Step( iOut1 );

iSliceCont2          = iSlice;
iSliceCont2( iOut2 ) = [];
d2ContP              = d2Cont;
d2ContP( iOut2 )     = [];
iSliceStep2          = iSlice( iOut2 );
d2StepP              = d2Step( iOut2 );


figure
subplot(2,2,1)
hold on, plot( d1,'k+-' )
hold on, plot( iSliceCont1,d1ContP,'bs' )
hold on, plot( iSliceStep1,d1StepP,'ro' )
ylabel( 'image-to-image drift' )
title( 'Y Axis')
subplot(2,2,2)
hold on, plot( d2,'k+-' )
hold on, plot( iSliceCont2,d2ContP,'bs')
hold on, plot( iSliceStep2,d2StepP,'ro')
title( 'X Axis')
subplot(2,2,3)
hold on, plot( d1ContCum,'b+-' )
hold on, plot( d1StepCum,'r+-' )
xlabel( 'slice' )
ylabel( 'Cummulative drift' )
subplot(2,2,4)
hold on, plot( d2ContCum,'b+-' )
hold on, plot( d2StepCum,'r+-' )
xlabel( 'slice' )


% % PLOT . check consistency
% 
% % original cummulative drift
% d1Cum = cumsum(d1);
% d2Cum = cumsum(d2);
% 
% figure
% subplot(1,2,1)
% hold on, plot(d1Cum,'b+-')
% hold on, plot(d1Comb,'ro')
% subplot(1,2,2)
% hold on, plot(d2Cum,'b+-')
% hold on, plot(d2Comb,'ro')



% stop timer
elapsedTime = toc;
% 1025 sec (40 cores)

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

info.ContDrift.Session          = info.Session;
ifo.ContDrift.processStep       = processStep;
info.ContDrift.processAppliedTo = processStep - 1;
info.ContDrift.Date             = timeFinished;
info.ContDrift.Note             = processNote;

info.ContDrift.FitNoise.hpar    = hpar;
info.ContDrift.FitNoise.ival1   = ival1;
info.ContDrift.FitNoise.nval1   = nval1;
info.ContDrift.FitNoise.fitOut1 = fitOut1;
info.ContDrift.FitNoise.d1Sig   = d1Sig;

info.ContDrift.d1Step      = d1Step;
info.ContDrift.d1Cont      = d1Cont;
info.ContDrift.d1StepCum   = d1StepCum;
info.ContDrift.d1ContCum   = d1ContCum;
info.ContDrift.d1Comb      = d1Comb;

info.ContDrift.FitNoise.ival2    = ival2;
info.ContDrift.FitNoise.nval2    = nval2;
info.ContDrift.FitNoise.fitOut2  = fitOut2;
info.ContDrift.FitNoise.d2Sig    = d2Sig;

info.ContDrift.d2Step      = d2Step;
info.ContDrift.d2Cont      = d2Cont;
info.ContDrift.d2StepCum   = d2StepCum;
info.ContDrift.d2ContCum   = d2ContCum;
info.ContDrift.d2Comp      = d2Comb;


% SAVE . update process list

iProc                 = processStep-1;

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

newProcess.Data       = "ContDrift";
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



%%  0 . GLOBAL DRIFT - Spline fit of cummulative drift

% TASK  Remove long range drift with predefined smoothing parameter

% NEED  results from previous section

% NOTE  To estimate the smoothing factor generate a variable from the 
%       arrays 'info.ContDrift.d1ContCum' and/or 'info.ContDrift.d2ContCum'
%       Use the fitting tool 'cftool' and reduce the fitting parameter
%       until the fit is smooth across several image slice numbers

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number / name
processStep = 0;
processName = "GLOBAL DRIFT - Spline fit of cummulative drift";

% process note
processNote = "";


% UPDATE -- Input/Output ------------------------------------------------

% UPDATE -- Process -----------------------------------------------------

% smoothing parameter axis 1 (y), axis 2 (x)
pSmooth1 = 0.001;
pSmooth2 = 0.001;

% -----------------------------------------------------------------------

fprintf( 'Process step %u: GLOBAL DRIFT - Spline fit of cummulative drift ...\n', processStep );

% start timer
tic


% drift, jumps removed
d1 = transpose( info.ContDrift.d1ContCum );
d2 = transpose( info.ContDrift.d2ContCum );

% set up fit
ft = fittype( 'smoothingspline' );


% fit axis 1

opts1 = fitoptions( 'Method', 'SmoothingSpline' );
opts1.SmoothingParam = pSmooth1;

[xData1, yData1] = prepareCurveData( [], d1 );
[fitresult1, ~]  = fit( xData1, yData1, ft, opts1 );
d1Fit            = feval(fitresult1,xData1 );


% fit axis 2

opts2 = fitoptions( 'Method', 'SmoothingSpline' );
opts2.SmoothingParam = pSmooth2;

[xData2, yData2] = prepareCurveData( [], d2 );
[fitresult2, ~]  = fit( xData2, yData2, ft, opts2 );
d2Fit            = feval(fitresult2,xData2 );



% stop timer
elapsedTime = toc;
% 1025 sec (40 cores)

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

info.SplineFit.Session          = info.Session;
info.SplineFit.processStep      = processStep;
info.SplineFit.processAppliedTo = processStep - 1;
info.SplineFit.Date             = timeFinished;
info.SplineFit.Note             = processNote;

info.SplineFit.FitType  = 'smoothingspline';
info.SplineFit.opts1    = opts1;
info.SplineFit.opts2    = opts2;

info.SplineFit.xData1   = xData1;
info.SplineFit.yData1   = yData1;
info.SplineFit.d1Fit    = d1Fit;

info.SplineFit.xData2   = xData2;
info.SplineFit.yData2   = yData2;
info.SplineFit.d2Fit    = d2Fit;


% SAVE . update process list

iProc                 = processStep-1;

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

newProcess.Data       = "SplineFit";
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


% PLOT . fit results for both axes

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


fprintf( '\n ... Finished successfully\n' );



%%  0 . GLOBAL DRIFT . Correct drift  >> Gdc

% TASK  remove long-range drift and correct only image-to-image shift

% NEED  results from previous section

% UPDATE -- WORKFLOW ----------------------------------------------------

% process step / name
processStep = 0;
processName = "GLOBAL DRIFT . Correct drift";

% process name acronym (3 characters, start with capital letter)
processSave = "Gdc";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders (true,false)
useFolders = true;

if useFolders == true

    projectFolderRead = info.Project.projectFolder;
    projectFolderSave = projectFolderRead;

end

% UPDATE -- Process -----------------------------------------------------

% crop to smallest commmon area (true/false)
cropFrame = false;


% -----------------------------------------------------------------------


fprintf( 'Process step %u: GLOBAL DRIFT . Correct drift ...\n', processStep );

% start timer
tic


% % remove long-range drift
ds1 = d1Comb - transpose( info.SplineFit.d1Fit );
ds2 = d2Comb - transpose( info.SplineFit.d2Fit );

% shift vector
shft  = [ds1; ds2 ];
shftr = round( shft );

if cropFrame == true

    nx = info.ProcessList.nX( processStep-1 );
    ny = info.ProcessList.nY( processStep-1 );

    cxS = - min( shftr(2,:) ) + 1;
    cxE = nx - max( shftr(2,:) );

    cyS = - min( shftr(1,:) ) + 1;
    cyE = ny - max( shftr(1,:) );

end


if useFolders == true

    % slice names
    sliceName = info.Project.sliceName;

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
        % for i=2:2

        if verbose == true
            PrintIndex(i,iEnd);
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % image index in array
        k = i - iStart + 1;

        % apply drift
        ITS = circshift( IT,-shftr(:,k) ); %#ok<PFBNS> 

        % save back
        % ftest = fullfile( projectFolder,sliceFolder,'test.tif' );
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( ITS,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    end

else

    if exist( 'IMSP','var' )

        % reset image series variable from previous process
        IMS = IMSP;
        clearvars IMSP

    end

    % shifted image series
    IMSP = zeros(size(IMS),'like',IMS);

    % number of slices
    [~,~,nSlices] = size( IMS );


    parfor i=1:nSlices

        if verbose == true
            PrintIndex( i,nSlices );
        end

        IMSP(:,:,i) = circshift( IMS(:,:,i),-shftr(:,i) );

    end
end

% stop timer
elapsedTime = toc;
% 1025 sec (40 cores)

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

info.ContDriftValues.Shift = shft;
info.ContDriftValues.CircShift = -shftr;


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

newProcess.Data       = "ContDriftValues";
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


if useFolders == false

    % PLOT . Compare before and after

    sz = size( IMS );
    ix = round( 0.5*sz(2) );
    plim(cat(2,squeeze(IMSP(:,ix,:)),squeeze(IMSP(:,ix,:))))

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

% last process step
processStep = 0;

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



%%   0 . Y-DEP DRIFT . Cross correlate stripes

% TASK  calculate drift for each stripe

% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 0;
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
ccIn.YNleastCommonArea = 'yes';  % use only common non-0 area for cc    *1
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
ccIn.idxRef = 500;

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
sHalf = 87;

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




%% 0 . Y-DEP DRIFT . Spline fitting

% TASK  separate statistical image-to-image shifts from long range
%       cummulative drifts induced by changes of the image content as
%       function of z

% NOTE  To estimate the smoothing factor extract an array from the 
%       variable ccC.dflcum and save into a variable. Use the fitting tool
%       'cftool' and reduce the fitting parameter until the fit is smooth
%       across several image slice numbers


% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 0;
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



%%  0.  Y-DEP DRIFT . Interpolate sampled drift

% TASK  prepare displacement matrix with pixel resolution

% PROC  1. Collect drift values as function of central y-position for
%          stripes
%       2. Interpolate drift values to pixel resolution
%       3. Plot interpolated x,y displacements for analysis


% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 0;
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


%%  0 . Y-DEP DRIFT . Warp images  >> Ydd

% TASK  calculate displacement matrix and apply imwarp()

% NEED  all sections of Y-DEP DRIFT group

% UPDATE -- WORKFLOW ----------------------------------------------------

% process step number / name
processStep = 0;
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
%%  0 . MISC . Convert 16 to 8 bit  >> Dyn

% TASK  reduce image size by reducing dynamic range from 16 to 8 bit

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number and name
processStep = 0;
processName = "MISC . Convert 16 to 8 bit";

% process name acronym (3 characters, start with capital letter)
processSave = "Dyn";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------


% use folders (true/false)
useFolders = true;

if useFolders == true

    % input/output
    projectFolderRead      = info.Project.projectFolder;
    projectFolderSave      = info.Project.projectFolder;

end


% -----------------------------------------------------------------------

fprintf( 'Process Step %u: Convert 16 to 8 bit ...\n', processStep );

% start timer
tic


% PARAMETERS

% scale factors
s8  = double(intmax( 'uint8' ));
s16 = double(intmax( 'uint16' ));


if useFolders == true
    
    % slice names
    sliceName = info.Project.sliceName;

    % slice folder names
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end


    % start and end slices after previous process
    iStart  = info.ProcessList.iStart(processStep-1);
    iEnd    = info.ProcessList.iEnd(processStep-1);


    for i=iStart:iEnd

        if verbose == true
            PrintIndex( i,iEnd );
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        if i==iStart
            if isa( IT,'uint8' )
                error(' Image series already 8 bit' )
            end
        end

        % convert 16 to 8 bit
        IC = uint8( s8/s16 .* double( IT ) );

        % save
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( IC,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    end

else
    
    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    if isa( IMS,'uint8' )
        error(' Image series already 8 bit' )
    end

    % allocate memory for resulting image series of this setion
    IMSP = zeros( size(IMS),'uint8' );

    % number of images
    [~,~,ni] = size( IMS );


    for i=1:ni

        if verbose == true
            PrintIndex( i,ni );
        end

        % convert image
        IMSP(:,:,i) = uint8( s8/s16 .* double( IMS(:,:,i) ));
    
    end

end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

% copy dimensions (that don't change in this process)
newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.nX       = info.ProcessList.nX( processStep-1 );
newProcess.nY       = info.ProcessList.nY( processStep-1 );
newProcess.nZ         = info.ProcessList.nZ( processStep-1 );

newProcess.Data     = "";
newProcess.Note     = processNote;     
newProcess.Date     = timeFinished;
newProcess.Duration = elapsedTime;

% info.ProcessList(pStep,:) = struct2table( newProcess );
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



%%  0 . MISC . Flip images along y  >> Fly

% TASK  Reverse image along y-axis in case sample sits upside down
%       (e.g. high-pressure frozen and freeze-substituted samples)

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number and name
processStep = 0;
processName = "MISC . Flip images along y";

% process name acronym (3 characters, start with capital letter)
processSave = "Fly";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders (true/false)
useFolders = true;

if useFolders == true

    % input/output
    projectFolderRead = info.Project.projectFolder;
    projectFolderSave = info.Project.projectFolder;
    
end


% -----------------------------------------------------------------------

fprintf( 'Process Step %u: Flip images along y ...\n', processStep );

% start timer
tic
    

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

    iStart = info.ProcessList.iStart(processStep-1);
    iEnd   = info.ProcessList.iEnd(processStep-1);


    for i=iStart:iEnd

        if verbose == true
            PrintIndex( i,iEnd );
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % flip first dimension of image
        IF = flipud( IT );
            
        % save
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( IF,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    end

else

    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    % flip first index (y-axis)
    IMSP = flip( IMS,1 );

end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

% copy dimensions (that don't change in this process)
newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.iRef     = info.ProcessList.iRef( processStep-1 );
newProcess.nX       = info.ProcessList.nX( processStep-1 );
newProcess.nY       = info.ProcessList.nY( processStep-1 );
newProcess.nZ         = info.ProcessList.nZ( processStep-1 );

newProcess.Data     = "";
newProcess.Note     = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

% info.ProcessList(pStep,:) = struct2table( newProcess );
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



%%  0 . MISC . Invert intensity  >> Inv

% TASK  Change from back-scatter to TEM bright-field like contrast

% NOTE  Section allows to stretch at the same time the dynamic range by
%       providing min/max intensities

% UPDATE -- WORKFLOW ------------------------------------------------------

% process number and name
processStep = 0;
processName = "MISC . Invert intensity";

% process name acronym (3 characters, start with capital letter)
processSave = "Inv";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders
useFolders = true;

if useFolders == true

    % input/output
    projectFolderRead = info.Project.projectFolder;
    projectFolderSave = info.Project.projectFolder;

end


% UPDATE -- PROCESS -----------------------------------------------------

% total min/max estimated from series histogram
iMin = 5000;
iMax = 57000;


% -----------------------------------------------------------------------

fprintf( 'Process Step %u: Invert images ...\n', processStep );

% start timer
tic

if useFolders == true

    % slice names
    sliceName       = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end

    iStart = info.ProcessList.iStart(processStep-1);
    iEnd   = info.ProcessList.iEnd(processStep-1);

    for i=iStart:iEnd

        if verbose == true
            PrintIndex(i,iEnd);
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        if i==iStart

            % image class (uint8 or uint16)
            imgClass = class( IT );

            % scale factors
            classMax = double(intmax( imgClass ));
            scl      = classMax/(iMax-iMin);

        end

        % crop min/max
        IT(IT<iMin) = iMin;
        IT(IT>iMax) = iMax;

        % invert
        switch imgClass

            case 'uint8'
                II = uint8( classMax - scl*(double( IT ) - iMin) );
    
            case 'uint16'
                II = uint16( classMax - scl*(double( IT ) - iMin) );

            otherwise
                error('Image class not implemented')
        end

        % save
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( II,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    end

else

    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    % allocate memory for resulting image series of this setion
    IMSP = zeros(size(IMS),'like',IMS);

    % image clase
    imgClass = class( IMS );
   
    % number of images
    [~,~,ni] = size( IMS );


    for i=1:ni

        if verbose == true
            PrintIndex( i,ni );
        end

        % current image
        IT = IMS(:,:,i);

         % crop min/max
        IT(IT<iMin) = iMin;
        IT(IT>iMax) = iMax;
       

        % invert
        switch imgClass

            case 'uint8'
                II = uint8( classMax - scl*(double( IT ) - iMin) );
    
            case 'uint16'
                II = uint16( classMax - scl*(double( IT ) - iMin) );

            otherwise
                error('image class not implemented')
        end
            
        % save
        IMSP(:,:,i) = II;

    end
end


% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

ninv = 0;
if isfield( info,'Invert' )

    % number of rotation events so far
    ninv = max(size( info.Invert ));

end

info.Invert( ninv+1 ).Session          = info.Session;
info.Invert( ninv+1 ).processStep      = processStep;
info.Invert( ninv+1 ).processAppliedTo = processStep - 1;
info.Invert( ninv+1 ).Date             = timeFinished;
info.Invert( ninv+1 ).Note             = processNote;

if useFolders == true
    info.Invert( ninv+1 ).sliceFolderRead  = sliceFolderRead;
    info.Invert( ninv+1 ).sliceFolderSave  = sliceFolderSave;
end

info.Invert( ninv+1 ).iMin = iMin;
info.Invert( ninv+1 ).iMax = iMax;


% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

% copy dimensions (that don't change in this process)
newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.nX       = info.ProcessList.nX( processStep-1 );
newProcess.nY       = info.ProcessList.nY( processStep-1 );
newProcess.nZ         = info.ProcessList.nZ( processStep-1 );

newProcess.Data     = "Invert";
newProcess.Note     = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

% info.ProcessList(pStep,:) = struct2table( newProcess );
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



%%  0 . MISC . Crop  >> Cro

% TASK  Crop to region of interest

% UPDATE -- WORKFLOW ------------------------------------------------------

% process number name
processStep = 0;
processName = "MISC . Crop";

% process name acronym (3 characters, start with capital letter)
processSave = "Cro";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders
useFolders = true;

if useFolders == true

    % input/output
    projectFolderRead = info.Project.projectFolder;
    projectFolderSave = info.Project.projectFolder;

end


% UPDATE -- PROCESS -----------------------------------------------------

% NOTE  Use values from ProcessList of previous process for those values
%       that stay uncahnged

% start and end pixel indeces
ixS = 1;
ixE = 5473;
iyS = 16;
iyE = 4358;

% new start and end image slice
iStart = info.ProcessList.iStart( processStep-1 );
iEnd   = info.ProcessList.iEnd( processStep-1 );


% -----------------------------------------------------------------------

fprintf( 'Process Step %u: Crop images ...\n', processStep );

% start timer
tic

% updated image series dimensions
nX = ixE - ixS + 1;
nY = iyE - iyS + 1;
nZ = iEnd - iStart + 1;

% start and end image in index space
izS = iStart - info.ProcessList.iStart( processStep-1 ) + 1;
izE = iEnd   - info.ProcessList.iStart( processStep-1 ) + 1;


if useFolders == true

    % slice names
    sliceName = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end


    for i=iStart:iEnd

        if verbose == true
            PrintIndex( i,iEnd );
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % crop
        ITC = IT( iyS:iyE, ixS:ixE );
            
        % save
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( ITC,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    end

else

    % !!! loaded image series must belong to process (processStep-1)

    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    % crop series
    IMSP = IMS( iyS:iyE, ixS:ixE, izS:izE );

    % updated image dimensions
    [nY,nX,nZ] = size( IMSP );

end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

ncrop = 0;
if isfield( info,'Crop' )

    % number of crop events so far
    ncrop = max(size( info.Crop ));

end

info.Crop( ncrop+1 ).Session          = info.Session;
info.Crop( ncrop+1 ).processStep      = processStep;
info.Crop( ncrop+1 ).processAppliedTo = processStep - 1;
info.Crop( ncrop+1 ).Date             = timeFinished;
info.Crop( ncrop+1 ).Note             = processNote;

if useFolders == true
    info.Crop( ncrop+1 ).sliceFolderRead  = sliceFolderRead;
    info.Crop( ncrop+1 ).sliceFolderSave  = sliceFolderSave;
end

info.Crop( ncrop+1 ).ixS = ixS;
info.Crop( ncrop+1 ).ixE = ixE;
info.Crop( ncrop+1 ).iyS = iyS;
info.Crop( ncrop+1 ).iyE = iyE;
info.Crop( ncrop+1 ).izS = izS;
info.Crop( ncrop+1 ).izE = izE;


% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

% update dimensions
newProcess.iStart   = iStart;
newProcess.iEnd     = iEnd;
newProcess.nX       = nX;
newProcess.nY       = nY;
newProcess.nZ       = nZ;

newProcess.Data     = "Crop";
newProcess.Note     = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

% info.ProcessList(pStep,:) = struct2table( newProcess );
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



%%  0 . MISC . Rotate  >> Rot

% TASK  Rotate images in case features (e.g. interfaces) are inclined

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number and name
processStep = 0;
processName = "MISC . Rotate";

% process name acronym (3 characters, start with capital letter)
processSave = "Rot";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders
useFolders = true;

if useFolders == true

    % input, output
    projectFolderRead      = info.Project.projectFolder;
    projectFolderSave      = projectFolderRead;

end


% UPDATE -- PROCESS -----------------------------------------------------

% angle
angle = -1.53;


% -----------------------------------------------------------------------

fprintf( 'ProcessStep %u: Rotate images by %3.2f degrees ...\n',...
    processStep,angle );

% start timer
tic


if useFolders == true

    % slice names
    sliceName = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end

    iStart = info.ProcessList.iStart(processStep-1);
    iEnd   = info.ProcessList.iEnd(processStep-1);


    for i=iStart:iEnd

        if verbose == true
            PrintIndex(i,iEnd);
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % save image class (uint8 or uint16)
        imgClass = class( IT );


        % rotate image
        IR = imrotate( double( IT ),angle,'bicubic','crop' );


        % save back
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );

        switch imgClass
            case 'uint8'
                imwrite( uint8(IR),fname,'tiff', ...
                    'Compression','packbits','WriteMode','overwrite')

            case 'uint16'
                imwrite( uint16(IR),fname,'tiff', ...
                    'Compression','packbits','WriteMode','overwrite')

            otherwise
                error('image class not implemented')
        end

    end

else

    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    % allocate memory for resulting image series of this setion
    IMSP = zeros(size(IMS),'like',IMS);

    % remember image class
    imgClass = class( IMS );

    % number of images
    [~,~,ni] = size( IMS );


    for i=1:ni

        if verbose == true
            PrintIndex( i,ni );
        end

        % current image
        IT = double( IMS(:,:,i) );

        % rotate image
        IR = imrotate( IT,angle,'bicubic','crop' );

        % save
        switch imgClass
            case 'uint8'
                IMSP(:,:,i) = uint8( IR );

            case 'uint16'
                 IMSP(:,:,i) = uint16( IR );

            otherwise
                error('image class not implemented')
        end

    end
end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

nrot = 0;
if isfield( info,'Rotate' )

    % number of rotation events so far
    nrot = max(size( info.Rotate ));

end

info.Rotate( nrot+1 ).Session          = info.Session;
info.Rotate( nrot+1 ).processStep      = processStep;
info.Rotate( nrot+1 ).processAppliedTo = processStep - 1;
info.Rotate( nrot+1 ).Date             = timeFinished;
info.Rotate( nrot+1 ).Note             = processNote;

if useFolders == true
    info.Rotate( nrot+1 ).sliceFolderRead  = sliceFolderRead;
    info.Rotate( nrot+1 ).sliceFolderSave  = sliceFolderSave;
end

info.Rotate( nrot+1 ).angle            = angle;
    

% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

% copy parameter that don't change
newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.nX       = info.ProcessList.nX( processStep-1 );
newProcess.nY       = info.ProcessList.nY( processStep-1 );
newProcess.nZ         = info.ProcessList.nZ( processStep-1 );

newProcess.Data     = "Rotate";
newProcess.Note     = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

% info.ProcessList(pStep,:) = struct2table( newProcess );
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



%%   0 . MISC . Bin  >> Bin

% TASK  Bin images

% UPDATE -- WORKFLOW ----------------------------------------------------

% process number and name
processStep = 0;
processName = "MISC . Bin";

% process name acronym (3 characters, start with capital letter)
processSave = "Bin";

% process note
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders
useFolders = true;

if useFolders == true

    % input, output
    projectFolderRead      = info.Project.projectFolder;
    projectFolderSave      = projectFolderRead;

end


% UPDATE -- PROCESS -----------------------------------------------------

% bin factor
bin = 2;


% -----------------------------------------------------------------------

fprintf( 'ProcessStep %u: Bin images by a factor of  %d \n',...
    processStep,bin );

% start timer
tic


if useFolders == true

    % slice names 
    sliceName = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end

    iStart = info.ProcessList.iStart(processStep-1);
    iEnd   = info.ProcessList.iEnd(processStep-1);


    for i=iStart:iEnd

        if verbose == true
            PrintIndex(i,iEnd);
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % bin image
        IB = TomoBinSeries( IT,bin,false );

        % save back
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( IB,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    end

    % new image dimensions
    [nY,nX] = size( IB );

else

    % overwrite IMS with previous result and remove previous result from
    % memory
    if exist( 'IMSP','var' )

        IMS = IMSP;
        clearvars IMSP

    end

    % TODO compile a mex file
    IMSP = TomoBinSeries( IMS,bin,verbose );

    % updated image size
    [nY,nX,~] = size( IMSP );

end

% end timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

nbin = 0;
if isfield( info,'Bin' )

    % number of rotation events so far
    nbin = max(size( info.Bin ));

end

info.Bin( nbin+1 ).Session          = info.Session;
info.Bin( nbin+1 ).processStep      = processStep;
info.Bin( nbin+1 ).processAppliedTo = processStep - 1;
info.Bin( nbin+1 ).Date             = timeFinished;
info.Bin( nbin+1 ).Note             = processNote;

if useFolders == true
    info.Bin( nbin+1 ).sliceFolderRead  = sliceFolderRead;
    info.Bin( nbin+1 ).sliceFolderSave  = sliceFolderSave;
end

info.Bin( nbin+1 ).factor = bin;


% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;

% update dimensions
newProcess.iStart   = info.ProcessList.iStart( processStep-1 );
newProcess.iEnd     = info.ProcessList.iEnd( processStep-1 );
newProcess.nX       = nX;
newProcess.nY       = nY;
newProcess.nZ       = info.ProcessList.nZ( processStep-1 );

newProcess.Data     = "Bin";
newProcess.Note     = processNote;     
newProcess.Date       = timeFinished;
newProcess.Duration   = elapsedTime;

% info.ProcessList(pStep,:) = struct2table( newProcess );
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

% Below template section contains the overall structure handling file
% input/ouput and organization of the workflow. In order to adapt the
% section do the following
%
% 1. replace place holders <...>
% 2. follow instructions '% TODO 1: ....' through '% TODO x: ...'



%% 0 . < PROCESS GROUP . Process name >  >>  < Acronym >

% TASK  < Longer description of given task >

% NEED  < list all conditions that need to be fulfilled for the section
%         to work. Use sparingly if at all. Preferably all necessary input
%         parameters are stored by the previous processes in the 'info'
%         structure. Remove if not needed. >


% UPDATE -- WORKFLOW ----------------------------------------------------

% process number and name
processStep = 0;
processName = "< process name >";

% preceding process step ( the process that this section builts on, 
% typically the directly preceding one in the process list. The rule can
% be broken if multiple data need to be prepared in intermediate process
% steps in order for this section to work, for example if a mask needs to
% loaded or prepared that is then applied to current image series.)
processPrec = processStep - 1;


% TODO 1: : DEFINE THREE-CHARACTER ACRONYM THAT IS USED FOR SAVING
%           INTERMITTENT RESULTS AND DESCRIBES THE UNDERLYING PROCESS
%           ( IF THIS PROCESS MODIFIES THE IMAGE SERIES, OTHERWISE LEAVE 
%           EMPTY )

% process name acronym ( 3 characters, start with capital letter )
% ( place the acronym in the section title to indicate that this section
% is modifying the image series )
processSave = "";


% process note ( appears in ProcessList table )
processNote = "";


% UPDATE -- IN/OUT ------------------------------------------------------

% use folders ( true/false )
useFolders = true;

if useFolders == true

    % TODO 2: CHANGE IN/OUT VARIABLES AS NEEDED,  IF IMAGE FOLDERS ARE 
    %         USED

    % input, output
    projectFolderRead = info.Project.projectFolder;
    projectFolderSave = projectFolderRead;

end


% UPDATE -- PROCESS -----------------------------------------------------

% TODO 3: DEFINE HERE INPUT PARAMETERS NEEDED FOR PROCESS, COLLECT DATA
%         FROM PRECEDING PROCESSES IN INFO STRUCTURE


% -----------------------------------------------------------------------

fprintf( 'ProcessStep %u: %s ...\n',processStep,processName );

% start timer
tic


% <PARAMETERS>

% TODO 4: DEFINE ADDITIONAL LOCAL PARAMETERS, IF NEEDED

% </PARAMETERS>


if useFolders == true

    % slice names
    sliceName = info.Project.sliceName;

    % input/output folders
    [sliceFolderRead,sliceFolderSave ] = SerialSliceFolderNames( ...
        info.Detector,processStep,info.ProcessList,processSave );

    % generate folder
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir( projectFolderSave,sliceFolderSave );
    if SUCCESS == 0
        error( '%s (ID %s)',MESSAGE,MESSAGEID )
    end

    iStart = info.ProcessList.iStart(processStep-1);
    iEnd   = info.ProcessList.iEnd(processStep-1);


    % NOTE: THE FOLLOWING FOR LOOP CAN BE REPLACED BY A PARFOR LOOP IF 
    %       VARIABLES FOLLOW THE NECESSARY RULES ( SEE MATLAB DOCU )

    for i=iStart:iEnd

        if verbose == true
            PrintIndex(i,iEnd);
        end

        % read current slice image
        fname = fullfile( projectFolderRead,sliceFolderRead,sliceName(i) );
        IT = imread( fname );

        % remember image class (uint8 or uint16)
        % NOTE: IN CASE IMAGES ARE TRANSFORMED TO CLASS DOUBLE
        % TO ENSURE PRECISION OF CALCULATIONS
        imgClass = class( IT );


        % <CODE>

        % TODO 5: PLACE FUNCTION OR SCRIPT BLOCK HERE

        % </CODE>
        
        
        % TODO 6: 1. CHANGE VARIABLE NAME OF PROCESSED IMAGE (IP here) AS
        %            PRODUCED IN THE ABOVE CODE BLOCK
        %         2. 'switch' BLOCK IS NOT NECESSARY IF DYNAMIC RANGE OF
        %            IMAGE IS NOT ALTERED BY THE CODE

        switch imgClass
            case 'uint8'
                IP = uint8( IP );

            case 'uint16'
                IP =  uint16( IP );

            otherwise
                error('image class not implemented')
        end


        % save
        fname = fullfile( projectFolderSave,sliceFolderSave,sliceName(i) );
        imwrite( IP,fname,'tiff', ...
            'Compression','packbits','WriteMode','overwrite')

    end

else

    if exist( 'IMSP','var' )

        % reset image series variable
        IMS = IMSP;
        clearvars IMSP

    end

    % allocate memory for resulting image series of this setion
    IMSP = zeros(size(IMS),'like',IMS);

    % remember image class
    imgClass = class( IMS );

    % number of images
    [~,~,ni] = size( IMS );


    for i=1:ni

        if verbose == true
            PrintIndex( i,ni );
        end

        % read current image
        % 
        IT = IMS(:,:,i);


        % <CODE>

        % TODO 6: PLACE FUNCTION OR SCRIPT BLOCK HERE
        %         CHANGE PROCESSED IMAGE NAME (HERE 'IP') AS NEEDED.
        %
        %         ALTERNATIVELY THE ENTIRE FOR LOOP CAN BE REPLACED BY A
        %         FUNCTION THAT HANDLES THE ENTIRE SERIES IMS

        % </CODE>


        % save
        switch imgClass
            case 'uint8'
                IMSP(:,:,i) = uint8( IP );

            case 'uint16'
                 IMSP(:,:,i) = uint16( IP );

            otherwise
                error('image class not implemented')
        end

    end
end

% stop timer
elapsedTime = toc;

% current time
timeFinished = datetime( 'now' );


% SAVE . update info structure

% TODO 8: CHANGE THE SUB-STRUCTURE NAME ( Here ProcessStruct ) TO FIT THE
%         DESIRED PROCESS

nInstance = 0;
if isfield( info,'ProcessStruct' )

    % number of process events so far
    nInstance = max(size( info.ProcessStruct ));

end

info.ProcessStruct( nInstance+1 ).Session          = info.Session;
info.ProcessStruct( nInstance+1 ).processStep      = processStep;
info.ProcessStruct( nInstance+1 ).processAppliedTo = processStep - 1;
info.ProcessStruct( nInstance+1 ).Date             = timeFinished;
info.ProcessStruct( nInstance+1 ).Note             = processNote;

if useFolders == true
    info.ProcessStruct( nInstance+1 ).sliceFolderRead     = sliceFolderRead;
    info.ProcessStruct( nInstance+1 ).sliceFolderSave     = sliceFolderSave;
end

% TODO 9: ADD CORRESPDONING INPUT VARIABLES AND RESULTS THAT ARE
%         IMPORTANT TO REMEMBER, E.G.
%
% info.ProcessStruct( nInstance+1 ).angle            = angle;
    

% SAVE . add process to list

newProcess          = table;
newProcess.Session  = info.Session;
newProcess.Step     = processStep;
newProcess.Process  = processName;
newProcess.Save     = processSave;


% TODO 10: COPY DIMENSIONAL PARAMETERS THAT DON'T CHANGE FROM PREVIOUS
%          PROCESS AND MODIFY THOSE THAT HAVE CHANGED

iProcPrev = processStep-1;
newProcess.iStart   = info.ProcessList.iStart( iProcPrev );
newProcess.iEnd     = info.ProcessList.iEnd( iProcPrev );
newProcess.nX       = info.ProcessList.nX( iProcPrev );
newProcess.nY       = info.ProcessList.nY( iProcPrev );
newProcess.nZ       = info.ProcessList.nZ( iProcPrev );

newProcess.Data     = "ProcessStruct";
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




