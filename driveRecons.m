
% DRIVERECONS.M Perform reconstruction process on raw MSOT data.
%   driveRecons takes configuration inputs, including the MSOT file location
%   and target output directory, and reconstructs the entire dataset according
%   to the predefined configuration.
%
%

% The pipeline expects a settings file in the given output directory.

function [completedTask] = driveRecons(inputPath,outputPath)
    
    
    
    pipelineMeta = load(fullfile(inputPath,'pipeline_settings.mat'));
    assert(numel(fieldnames(pipelineMeta))==1); % Only allow a single variable in the file.
    if isfield(pipelineMeta,'pipelineMeta')
        pipelineMeta = pipelineMeta.pipelineMeta;
    end

    
    %% Add the JARs for the MSOT schemas
    util.addJavaLibraries;
    
    
    %% Debug stuff
    disp(['Working in ' pwd]);
    disp(pipelineMeta.loaderSettings)
    DEBUG = pipelineMeta.debugSettings.doDebug;
    
    %% Initial metadata preprocessing
    % Create the loader for the data and the metadata.
    
    loader = util.MSOTSignalLoader(pipelineMeta.loaderSettings);
    meta = loader.Meta;
    calculateStructure(meta);
    
    
    
    %% Set up prefilter.
    preFilter = recon.MSOTPreFilter(pipelineMeta.filterSettings);
    
    % Index list of which frames to reconstruct
    reconIndex = pipelineMeta.reconSettings.FramesToReconstruct; 
    NRecon = numel(reconIndex);
    
    
    % Connect to cluster, set up pool(s)
    if ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'));
    end
    
    % TODO: Put this in a separate configuration utility.
    clust = parcluster
    SPMD_RANGE = [1 max(clust.NumWorkers)]
%     SPMD_RANGE = [1 max(clust.NumWorkers)]
    
    % Create storage variables for reconstructed data.
%     reconImageStructure = zeros( pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x , NRecon);
%     reconLogStructure = cell(NRecon,1);
    
    
    %% Perform reconstruction.
    clientJavaPath = javaclasspath('-dynamic');
    
    switch pipelineMeta.parallelSettings.parallelType
        case 'single'
            %% Set up reconstruction processing for one worker.
            reconstructor = recon.ReconSystem(pipelineMeta.reconSettings);
            reconstructor.constructModel;
            
            % Output reconstruction files.
            
            dataWriter = filter.BinWriter(...
                outputPath,...
                'recon_images_1',...
                [pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x],...
                'double',...
                NRecon);
            

            
            % Loop over reconstructions.
            for p = reconIndex
                
                if DEBUG
%                     disp(['Iteration ',num2str(p-indexStart+1)]);
                    disp(['Reconstructing frame ',num2str(p)])
                end
                
                % Load data, get the raw signal.
                dataChunk = loader(p,{});
                
                % Filter the signal
                filteredSignal = preFilter(dataChunk);
                
                % Perform reconstruction.
                [reconVector,reconMessage] = reconstructor(filteredSignal);
                reconLog.message = reconMessage;
                
                % Format the reconstruction log.
                metaStructure = struct(dataChunk.Meta);
                reconLog.header = metaStructure;
                
                % Put in the storage structure.
                
                recFrame.Data = reconVector;
                recFrame.Meta = reconLog;
                dataWriter(recFrame,p);
            end
            release(dataWriter);
            % Save data.
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'local'
            %% Set up processing for local workers
            % TODO: Put this in configureParallel.
            workerPool = parpool(clust,SPMD_RANGE(end),'AutoAddClientPath',true)
            % Create output file.
            %% SPMD
            
            output = Composite(); % This is a SUPER dirty way to do this.
            spmd(SPMD_RANGE(end))
                warning('off','MATLAB:structOnObject')
                disp(['Entered parallel session in folder ' pwd]);
                
                
                try
                    util.addJavaLibraries;
                catch
                    javaaddpath(clientJavaPath);
                end
                disp('Added Java Beans path');
                
                
                % Set up debug monitors if set to do so.
                if DEBUG
                    mpiprofile off; % Initialize profiling.
                    mpiprofile on -detail builtin;
                end
                
                % Give each worker a portion of the reconstruction icDMM(ndices
                indexList = codistributed(reconIndex);
                codistr = getCodistributor(indexList);
                
                if DEBUG
                    parallelWorkerList = codistributed(ones(1,numlabs));
                    codistrParallelWorkers = getCodistributor(parallelWorkerList);
                end
                
                % Construct a codist array similar to the index list, but a cell
                outputList = codistributed(cell(size(indexList)));
                
                % Get the local part of the recon index and output lists.
                localIndices = getLocalPart(indexList);
                
                localLoader = util.MSOTSignalLoader(pipelineMeta.loaderSettings);
                
                localReconstructor = recon.ReconSystem(pipelineMeta.reconSettings);
                localReconstructor.constructModel;
                
                % Determine the location of each local index in the overall
                % global list. Load reconstructor into local space.
                NLocalRecon = numel(localIndices);
                globalIndex = globalIndices(indexList,2);
                
%                 localImages = zeros( pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x, NLocalRecon);
%                 localLogs = cell(NLocalRecon,1);
                
                workerReconName = ['recon_images_',num2str(labindex)];
                workerReconFilename = fullfile(outputPath,[workerReconName '.mat'])
                workerLogFilename = fullfile(outputPath,['recon_logs_',num2str(labindex),'.mat'])
                workBinName = strrep(workerReconFilename,'.mat','.bin')
                
                
                % Create .mat file for each worker in the host directory. Since
                % this is local, network transfer shouldn't be an issue, but
                % hard disk writing might be.
                workReconFile = matfile(workerReconFilename,'Writable',true);
               
                dataWriter = filter.BinWriter(...
                        outputPath,...
                        workerReconName,...
                        [pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x],...
                        'double',...
                        NLocalRecon);
                
                
                %% Perform reconstructions on each worker.
                %
                for p = 1:NLocalRecon
                    
                    if DEBUG
                        disp(['reconstructing frame ',num2str(localIndices(p))]);
                    end
                    
                    % Load the signal. TODO: Have a buffer set up so that
                    % the workers can load big chunks of data at once to
                    % avoid network thrashing. This will probably require
                    % switching from a for to a while loop.
                    dataChunk = localLoader(localIndices(p),{});
                    rawSignal = double(dataChunk.Data);
                    
                    % Filter the signal. If in chunk, can do this with all.
                    filteredSignal = preFilter(dataChunk);
                    
                    % Reconstruct the image. If in chunk, loop over chunk.
                    [reconVector,reconMessage] = localReconstructor(filteredSignal);
                    
                    
                    reconLog.message = reconMessage;
                    metaStructure = struct(dataChunk.Meta);
                    % Append additional information about the parallelization.
                    % TODO: How to generalize this?
                    reconLog.header = metaStructure;
                    reconLog.time = dataChunk.Meta.RelTime;
                    reconLog.WL = dataChunk.Meta.Wavelength;
                    reconLog.workerIterationIndex = p;
                    reconLog.reconstructionIndex = localIndices(p);
                    reconLog.globalIndex = globalIndex(p);
                    
                    % Store. If in chunk, write as block.
                    
%                     localImages(:,:,p) = reconVector;
%                     localLogs{p} = reconLog;
                    
                    recFrame.Data = reconVector;
                    recFrame.Meta = reconLog;
                    
                    dataWriter(recFrame,p);
                   if mod(p,20)==0
			disp(p);
		   end 
                end
                
                release(dataWriter);
                disp('Reconstruction Complete');
                
%                 workReconFile.reconImages = localImages;
%                 workReconFile_bin.Data.reconImages = localImages;
%                 workLogFile.reconLogs = localLogs;
                
                
                if DEBUG
                    pInfo = mpiprofile('info');
                    pInfoOutput = codistributed.build(pInfo,codistrParallelWorkers);
                end
                output{labindex} = {workerReconFilename,workerLogFilename}; 
                
                
                
                
            end
            
            
            clientCollect = gather(output);
            disp('Data pulled from workers');
            if DEBUG
                profileData = gather(pInfoOutput);
            end
            
            
            
            delete(workerPool);
            
            
            %%%%%%
        case 'distributed'
            % Set up main client for submitting jobs.
            
        otherwise
    end
    
    
    
    %% Save to the local work directory, or the original host if being called from not-nextflow.
    
    if ~isdir(outputPath)
        mkdir(outputPath);
    end
    
    whos % Post out the variables currently in the workspace so that we can debug later.
    
    if DEBUG
        try
            save(fullfile(outputPath,'recon_parallel_debug.mat'),'profileData','-v7.3');
        catch
        end
    end
    
    metaStruct = struct(meta);
    pause(5);
    save(fullfile(outputPath,'recon_metadata.mat'),'metaStruct','-v7.3');
    pause(5);
    save(fullfile(outputPath,'recon_pipeline_settings.mat'),'pipelineMeta','-v7.3');
    pause(5);
    
    completedTask = fullfile(outputPath);
    
    pause(5); % Give any file writing a fair chance of finishing.
    
    
end





























