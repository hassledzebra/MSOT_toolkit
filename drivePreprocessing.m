% DRIVEPREPROCESSING.M Format pipeline inputs and configuration.
%   This step should precede all others and should generate a configuration
%   file which dictates the operation of the rest of the pipeline. 
%

function [completedTask] = drivePreprocessing(inputPath,outputPath)
    
    try
        jsonInputs = util.JSON2struct(fullfile(inputPath,'params.json'));
    catch
        % TODO: Gracefully handle this by loading default pipeline params.
        error('Something went wrong when loading pipeline parameters');
    end
    
    
    % Add java paths. 
    util.addJavaLibraries;
    
    % Resolve any variable/struct assignments, and create a struct without
    % redundant variable names. 
    [resolvedNames,overwriteStructs] = cleanInputs(jsonInputs); 
    resolvedNames.inputDirectory = inputPath;
    resolvedNames.outputDirectory = outputPath;
    
    disp('Initial pipeline meta creation');
    pipelineMeta = formatInputs(resolvedNames,overwriteStructs); % Assign variables to where they belong.
    save('debug_dump.mat','pipelineMeta'); % Save the created pipelineMeta before any further fill-in so we can debug later.
    
    % Create a temporary loader. This serves the purpose of helping the overall
    % pipeline fail quickly if there's something wrong with the files as
    % provided, and also automatically links up an msotData object which we can
    % use to populate settings.
    loader = util.MSOTSignalLoader(pipelineMeta.loaderSettings);
    meta = loader.Meta;
    calculateStructure(meta);
    
    
    % Fill in any settings which inherit from meta (e.g. sampling rate)
    disp('Populating from data meta');
    pipelineMeta = util.settingsFromMeta(pipelineMeta,meta);
    
    % Degenerate means of specifying which frames to reconstruct, but a
    % good placeholder for now. 
    if ~isfield(pipelineMeta.reconSettings,'FramesToReconstruct')
        pipelineMeta.reconSettings.FramesToReconstruct = 1:numel(meta.ScanFrames);
    else
        if ischar(pipelineMeta.reconSettings.FramesToReconstruct)
            pipelineMeta.reconSettings.FramesToReconstruct = jsondecode(pipelineMeta.reconSettings.FramesToReconstruct);
        end
        pipelineMeta.reconSettings.FramesToReconstruct = pipelineMeta.reconSettings.FramesToReconstruct(:).';
    end
    
    
    switch pipelineMeta.reconSettings.TuneSpeedOfSound % TODO: Factor this out into another function which can just work off the pipelineMeta.
        
        case 'auto' % Automatically tune the speed of sound.
            [pipelineMeta.reconSettings.SpeedOfSound,tuningResults] = recon.backProjectTune(...
                loader , meta , max(300,pipelineMeta.reconSettings.N_x) , pipelineMeta.reconSettings.FieldOfView_X);
            pipelineMeta.debugSettings.DIAGNOSTICS.SpeedOfSoundTuning = tuningResults;
        case 'settings' % Use the speed of sound in the .nod file, which reflects user settings during the scan.
            try
            loadedNod = util.nodImport(pipelineMeta.loaderSettings.NodFile);
            
            pipelineMeta.reconSettings.SpeedOfSound = 1520 + loadedNod.TrimSpeedOfSound;
            catch
                 pipelineMeta.reconSettings.SpeedOfSound = 1520 + meta.MeasurementDesc.TrimSOS;
                if meta.MeasurementDesc.TrimSOS == 0
                    warning('Speed of sound at default of 1520m/s! Was it tuned during acquisition?');
                end
            end
            
        case 'user' % Defined by the user input.
            if ~isnumeric(pipelineMeta.reconSettings.UserSpeedOfSound)
                pipelineMeta.reconSettings.SpeedOfSound = str2double(pipelineMeta.reconSettings.UserSpeedOfSound);
            else
                pipelineMeta.reconSettings.SpeedOfSound = pipelineMeta.reconSettings.UserSpeedOfSound;
            end
                
            
        otherwise
            % It should be set to a number or an image of the speed of
            % sound distribution. If not, we'll run into an error pretty
            % fast.
            error('Invalid speed of sound configuration');
    end
    
    
    
    % Import spectra, overwriting the current unmixing settings. 
    pipelineMeta.unmixSettings = util.importSpectra(pipelineMeta.unmixSettings);
 
    
    save(fullfile(outputPath,'pipeline_settings.mat'),'pipelineMeta','-v7.3');
    
    drawnow;pause(1);drawnow; % Flushes buffers to avoid crashes due to incomplete writing of mat files.
    
    completedTask = fullfile(outputPath,'pipeline_settings.mat');
    pause(2); % Give any file writing a fair chance of finishing.
    
    
end


% Function which cleans up the import from JSON files created by Nextflow. 
% Generally just makes sure that the Nextflow parameters (which are doubled
% up using camelcase and underscore names) are reduced, making sure that we
% take care of any params which we want to use for struct assignment later.

function [cleanedStruct,dirtyStructs] = cleanInputs(structIn)
    
    structFields = fieldnames(structIn);
    structValues = struct2cell(structIn);
    
    structCat = [structFields structValues];
    
    hasUnderscores_mask = cellfun(@(x) contains(x,'_'),structFields);
    
    hasUnderscores = structCat(hasUnderscores_mask,:);
    camelcaseGroup = structCat(~hasUnderscores_mask,:);
    
    
    
    
    % Camelcase form will have no underscores UNLESS there's a _DOT_ flag.
    hasDotFlag_mask = cellfun(@(x) numel(x)>1,regexp(hasUnderscores(:,1),'_DOT_','split'));
    
    hasDotFlag = hasUnderscores(hasDotFlag_mask,:);
    
    % Underscore form will have no caps UNLESS there's a dot, in which case
    % the immediately following character will be capitalized.
    hasCapsFlag = regexp(hasUnderscores(:,1),'[A-Z]','split');
    hasCapsFlag_mask = cellfun(@(x) numel(x)>1,hasCapsFlag);
    hasUnderscoreWithCaps_mask = (hasCapsFlag_mask & ~hasDotFlag_mask);
    
    hasUnderscoreWithCaps = hasUnderscores(hasUnderscoreWithCaps_mask,:);
    
    % Clean up. All underscores and camelcase without issues should be in
    % one group now.
    underscoreGroup = hasUnderscores((~hasCapsFlag_mask & ~hasDotFlag_mask),:);
    recameledGroup = underscoreGroup;
    recameledGroup(:,1) = regexprep(underscoreGroup(:,1),'_([a-z])','${upper($1)}');
    
    
     
    underscoreWithCaps_renamed = hasUnderscoreWithCaps;
    underscoreWithCaps_repdot(:,1) = regexprep(hasUnderscoreWithCaps(:,1),'_([A-Z])','.$1');
    underscoreWithCaps_renamed(:,1) = regexprep(underscoreWithCaps_repdot(:,1),'_([a-z])','${upper($1)}');
    
    hasDot_renamed = hasDotFlag;
    hasDot_renamed(:,1) = regexprep(hasDotFlag(:,1),'(_DOT_)','.');
    
    cleanedStruct = cell2struct(recameledGroup(:,2),recameledGroup(:,1),1);
   
    dirtyStructs = makeStructuresFromParameters(hasDot_renamed(:,1),hasDot_renamed(:,2));
    if isempty(fieldnames(cleanedStruct)) && ~isempty(structCat)
        structCat = structCat.';
        for q = 1:numel(structCat)
            if iscell(structCat{q})
                structCat{q} = {structCat{q}};
            end
        end
        cleanedStruct = struct(structCat{:});
    end
end



% Format the inputs from the inputParser to make sure that all necessary
% parameters are either set to defaults or to whatever the user entered.
function formattedStructure = formatInputs(parsedInputs,overwriteStructs)
    
    addpath(pwd)
    
    %     unmatched = parsedInputs.Unmatched; %TODO: More extensive handling of unmatched entries.
    %     results = parsedInputs.Results;
    
    makeFullName = @(x) fullfile(x.folder,x.name);
    
    %% Handle inputs.
    msotDir = dir(fullfile(parsedInputs.inputDirectory,'*.msot'));
    binDir = dir(fullfile(parsedInputs.inputDirectory,'*.bin'));
    irfDir = dir(fullfile(parsedInputs.inputDirectory,'*.irf'));
    nodDir = dir(fullfile(parsedInputs.inputDirectory,'*.nod'));
    
    
    loaderSettings.MsotFile = fullfile(parsedInputs.inputDirectory,msotDir.name);
    loaderSettings.BinFile = fullfile(parsedInputs.inputDirectory,binDir.name);
    loaderSettings.IrfFile = fullfile(parsedInputs.inputDirectory,irfDir.name);
    
    try
        loaderSettings.NodFile = fullfile(parsedInputs.inputDirectory,nodDir.name);
    catch
        loaderSettings.NodFile = '';
    end
    
    
    filterSettings = struct();
    
    reconSettings.ModelType = parsedInputs.reconModel;
    reconSettings.ModelSolver = parsedInputs.reconSolver;
    
    reconSettings.FieldOfView_X = getStrNum(parsedInputs.fieldOfView)./1000; % Convert mm to m
    reconSettings.FieldOfView_Y = getStrNum(parsedInputs.fieldOfView)./1000;
    
    reconSettings.N_x = getStrNum(parsedInputs.pixelResolution);
    reconSettings.N_y = getStrNum(parsedInputs.pixelResolution);
    
    reconSettings.TuneSpeedOfSound = parsedInputs.tuneSpeedOfSound;
    reconSettings.UserSpeedOfSound = parsedInputs.userSpeedOfSound;
    
    
    % Move all endmember names into a standard form. 
    unmixSettings.EndmemberNames = util.fixEndmembers(parsedInputs.unmixingEndmembers);
    
    try
    	unmixSettings.SpectraFiles = regexp(parsedInputs.spectraFiles,'(\w+?)Spectrum\.','tokens');
        if iscell(unmixSettings.SpectraFiles{1})
                newUnmixSpectraFiles = {};
            for k = 1:numel(unmixSettings.SpectraFiles)
                newUnmixSpectraFiles{k} = unmixSettings.SpectraFiles{k}{1};
            end
            unmixSettings.SpectraFiles = newUnmixSpectraFiles;
        end
    catch
        unmixSettings.SpectraFiles = {};
    end
    
    
    
    % If we do have any user spectra files, check to see if we need to
    % override any default spectra.
    if ~isempty(unmixSettings.SpectraFiles)
        util.struct2JSON(fullfile('.','spectral'),'.','USER_SPECTRA_FOLDER');
        spectraFile_parsedNames = cat(2,unmixSettings.SpectraFiles{:});
        
        unmixSettings.EndmemberNames = union(spectraFile_parsedNames,unmixSettings.EndmemberNames);
        
    end
    unmixSettings.UnmixSolverType = parsedInputs.unmixingSolver;
    unmixSettings.MultispectralStateFilterType = parsedInputs.multispectralPrefilter;
    parallelSettings.parallelType = parsedInputs.parallelizationSettings;
    
    
    
    outputSettings = struct();
    
    debugSettings.doDebug = parsedInputs.doDebug;
    
    
    try parallelSettings = unmatched.parallelSettings; catch end
    try debugSettings = unmatched.debugSettings; catch end
    
    
    
    
    %% Load the default structure.
    formattedStructure = util.loadDefault('PipelineSettings');
    
    loadSet = util.loadDefault('LoaderSettings');
    formattedStructure.loaderSettings = fuseStructures(loadSet,loaderSettings);
    
    filterSet = util.loadDefault('FilterSettings');
    formattedStructure.filterSettings = fuseStructures(filterSet,filterSettings);
    
    reconSet = util.loadDefault('ReconSettings');
    formattedStructure.reconSettings = fuseStructures(reconSet,reconSettings);
    
    unmixSet = util.loadDefault('UnmixSettings');
    formattedStructure.unmixSettings = fuseStructures(unmixSet,unmixSettings);
    
    parallelSet = util.loadDefault('ParallelSettings');
    formattedStructure.parallelSettings = fuseStructures(parallelSet,parallelSettings);
    
    debugSet = util.loadDefault('DebugSettings');
    formattedStructure.debugSettings = fuseStructures(debugSet,debugSettings);
    
    outputSet = util.loadDefault('OutputSettings');
    formattedStructure.outputSettings = fuseStructures(outputSet,outputSettings);
    
    for k = 1:numel(overwriteStructs)
        formattedStructure.(overwriteStructs{k}.name) = fuseStructures(util.loadDefault(firstupper(overwriteStructs{k}.name)),overwriteStructs{k}.contents);
    end
    %% Handle the unmatched inputs.loadedNod.Trim
    %%TODO%%
    
    function str = firstupper(str)
        str(1)=upper(str(1));
    end
    function number = getStrNum(entry)
       if isnumeric(entry)
           number = entry;
       elseif ischar(entry) || isstring(entry)
               number = str2double(entry);
       else
           number = NaN;
       end
    end
end


% Given a list of input names and input values, create holder structures
% so that further parsing can take place if necessary.
function outputStructures = makeStructuresFromParameters(inputNames,inputValues)
    
    % Make sure all the input names are signified correctly.
    
    % Analyze the input names to determine how many structures to make.
    % Skip over the first character (as it's a dash);
    [splitNameList,splitDotList] = cellfun( @(x) split(x,'.'),inputNames,'UniformOutput',false);
    
    
    combinedLists = fuseStructureLists(splitNameList,splitDotList);
    structureNames = cellfun( @(x) x{1},combinedLists,'UniformOutput',false);
    
    [U,A,C] = unique(structureNames);
    outputStructures = cell(numel(U),1);
    
    for structureIndex = 1 : numel(U)
        outputStructures{structureIndex}.name = U{structureIndex};
        sameHeadInds = find(C == structureIndex);
        holdStruct = struct();
        
        for structureCoordinateIndex = 1 : numel(sameHeadInds)
            assignmentBundle = substruct(combinedLists{sameHeadInds(structureCoordinateIndex)}{2:end});
            holdStruct = subsasgn(holdStruct,assignmentBundle,inputValues{sameHeadInds(structureCoordinateIndex)});
        end
        outputStructures{structureIndex}.contents = holdStruct;
    end
    
    
    
    % Helper function that joins together reference lists for structures.
    function fusedLists = fuseStructureLists(names,delimiters)
        fusedLists = cell(size(names));
        for strucInd = 1:numel(fusedLists)
            holdList = cell(1,numel(names{strucInd})+numel(delimiters{strucInd}));
            holdList(1:2:(numel(names{strucInd})*2 - 1)) = names{strucInd};
            holdList(2:2:(numel(delimiters{strucInd})*2 + 1)) = delimiters{strucInd};
            fusedLists(strucInd) = {holdList};
        end
    end
    
    
end

% Joins two structures together to make another structure. structB
% dominates structA, so if structA and structB have the same fields,
% structB's value will be the one in fusedStructs.
function fusedStructs = fuseStructures(structA,structB)
    
    fusedStructs = structA;
    
    Bfields = fieldnames(structB);
    NBfields = numel(Bfields);
    
    for k = 1:NBfields
        try

            if isnumeric(structB.(Bfields{k})) || islogical(structB.(Bfields{k}))
                fusedStructs.(Bfields{k}) = structB.(Bfields{k});
                continue;
            end
            
            convVal = str2double(structB.(Bfields{k}));
            if isnan(convVal)
                switch structB.(Bfields{k})
                    case true
                    case false
                    case {'true',"true"}
                        convVal = true;
                    case {'false',"false"}
                        convVal = false;
                    otherwise
                        convVal = structB.(Bfields{k});
                end
            else
            end
        catch
            
            convVal = structB.(Bfields{k});
        end
        fusedStructs.(Bfields{k}) = convVal;
    end
    
end



