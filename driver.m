
% Main entrypoint for pipeline processing which passes through arguments from
% command line to the specified drive* code.

function [runCompleted] = driver(stepToRun,varargin)
    
    save('driverInput.mat','varargin');
    % Format output style and such.
    format compact;
    
    % Deployment-specific tasks.
    if isdeployed
        fprintf('CTF Root: %s\n',ctfroot);
    end
    
    % Status announcement.
    fprintf('Present Working Directory: %s\n',pwd);
    warning('off','MATLAB:structOnObject')
    
    
    % Add the path for the toolbox. In deployed sessions this will be in the
    % unpackaged archive directory. In Nextflow, who knows!
    driverPath = fileparts(mfilename('fullpath'));
    fprintf('Location of driver m-file: %s\n', driverPath);
    fprintf('\tAdding path: %s\n',driverPath);
    addpath(driverPath);
    
    testExternal = strsplit(genpath(fullfile(driverPath,'external')),':');
    
    if ~isempty(testExternal)
        disp('External paths:');
        for k = 1:numel(testExternal)
            fprintf('\tAdding path: %s\n',testExternal{k});
            addpath(testExternal{k});
        end
    end
    
    nVargs = numel(varargin);
    switch nVargs 
        case 0 % If no inputs, just use the current working directory.
            inputPath = pwd;
            outputPath = pwd;
        case 1
            if ~isempty(varargin{1}) % If one input and it isn't empty, assume it's a path and use it.
                inputPath = varargin{1};
                outputPath = varargin{1};
            else                     % Otherwise it's empty and we should just use the pwd
                inputPath = pwd;
                outputPath = pwd;
            end
        case 2 % Two inputs, first is the input and second is the output. 
            inputPath = varargin{1};
            outputPath = varargin{2};
        otherwise
            error("Incorrect number of inputs provided to driver.");
    end
    
    if inputPath=='\'
        warning(['Root was passed in as input directory. Reverting to ',pwd]);
        inputPath = pwd;
    end
    
    if outputPath=='\'
        warning(['Root was passed in as output directory. Reverting to ',pwd]);
        outputPath = pwd;
    end
    
    % Make sure the output path exists before any processing occurs.
    if ~isdir(outputPath)
        mkdir(outputPath);
    end
    
    fprintf('Input directory: %s',inputPath);
    fprintf('Output directory: %s\n',outputPath);
    
    fprintf('Java Heap Memory: %d bytes\n ',java.lang.Runtime.getRuntime.maxMemory);
    runCompleted = false;
    
    % General step should be some combination of str2func and feval.
    % TODO: Add check that if the varargin is the string 'test', the step
    % loads up some default data.
    
    switch stepToRun
        case 'drivePreprocessing'
            drivePreprocessing(inputPath,outputPath);
        case 'driveRecons'
            driveRecons(inputPath,outputPath);
        case 'collatePipelineOutput'
            collatePipelineOutput(inputPath,outputPath);
        case 'driveUnmixing'
            driveUnmixing(inputPath,outputPath);
        case 'driveHDF5Writing'
            driveHDF5Writing(inputPath,outputPath);
        case 'debug'
        otherwise
            error(strcat('Unrecognized pipeline step. Valid steps are ',...
                'drivePreprocessing, driveRecons, collatePipelineOutput,',...
                'driveUnmixing, or driveHDF5Writing'));
    end
    
    runCompleted = true;
    
end