function [sliceFolderIn,sliceFolderOut] = SerialSliceFolderNames(detector,process,list,acronym )
%  prepare folder names for consecutive processes
% -----------------------------------------------------------------------
%
%  Generate variables 'sliceFoderRead' and 'slcieFolderSave' for
%  SerialDataEval script. Use in case current process follows strictly
%  from previous one.
% 
%  SYNTAX   [sliceFolderIn,sliceFolderOut] = ...
%               SerialSliceFolderNames(detector,process,list,extension )
%
%  INPUT    detector    variable 'info.Detector' (string)
%           process     variable 'processStep' 
%           list        variable 'info.ProcessList' (table)
%           extension   variable 'processSave' (string)
% 
%  OUTPUT   sliceFolderIn   variable 'sliceFolderRead'
%           sliceFolderOut  variable 'sliceFolderSave'
%
% -----------------------------------------------------------------------

% convert to character character
detChar  = convertStringsToChars( detector );
acroChar = convertStringsToChars( acronym );

% concatenate process acronyms for previous step
procList = '';

for i=1:process-1
    procList = cat( 2,procList,...
        convertStringsToChars( list.Save(i) ));
end


% input folder name
if isempty( procList )

    % read from raw data folder
    sliceFolderIn = sprintf( '%s', detChar );

else

    sliceFolderIn = sprintf( '%s_%s', detChar,procList );

end


% output folder nmae
sliceFolderOut = sprintf( '%s_%s%s', detChar,procList,acroChar );


end
