function [collatedReconFile,collatedLogFile] = collatePipelineOutput(inputPath,outputPath)
    % collatePipelineOutput organizes the multiple output image and meta files from a
    % parallelized reconstruction job into monolithic image and meta files
    % suitable for use by driveUnmixing and driveHDF5Writing.
    
    
    disp('visible files are')
    dir(inputPath);
    load(fullfile(inputPath,'pipeline_settings.mat'));
    % Check for all files which are of the form 'recon_images_1.meta' or
    % 'recon_logs_2.meta' or similar.
    containedFiles = [dir([inputPath,filesep,'recon_images_*.bin']);dir([inputPath,filesep,'recon_images_*.meta'])];
    
    % Every image should have a matching .meta 
    assert(mod(numel(containedFiles),2)==0);
    nSigs = numel(containedFiles)/2;
    
    % If there's more than one signature, we used parallel reconstruction.
    if nSigs ~= 1
        Nims_overall = numel(pipelineMeta.reconSettings.FramesToReconstruct);
        
        dataWriter = filter.BinWriter(...
                        outputPath,...
                        'recon_images',...
                        [pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x],...
                        'double',...
                        Nims_overall);
        
        % For each of the worker recons, get a loader into it.
        for k = 1:nSigs
            thisLogList = util.JSON2struct(fullfile(inputPath,['recon_images_',num2str(k),'.meta']));
            Nims_part = thisLogList(end).workerIterationIndex;
            reconPartLoader = filter.BinLoader(...
                        inputPath,...
                        ['recon_images_',num2str(k)],...
                        [pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x],...
                        'double',...
                        Nims_part);
                    
            % Write the loaded data into the corresponding location in the
            % output file.
            for p = 1:Nims_part
                locFrame = reconPartLoader(p);
                dataWriter(locFrame,locFrame.Meta.globalIndex); 
                
                if mod(p,100)==0
                    fprintf('Writing of frame %d of part %d complete - global index %d\n',p,k,locFrame.Meta.globalIndex);
                end
            end
            
        end
    else % We used serial reconstruction and the file is already in the form we need, just with the wrong name. 
        
    
        copyfile(fullfile(inputPath,'recon_images_1.bin'),fullfile(outputPath,'recon_images.bin'));
        copyfile(fullfile(inputPath,'recon_images_1.meta'),fullfile(outputPath,'recon_images.meta'));
    end
    
    collatedReconFile = fullfile(outputPath,'recon_images.mat');
    collatedLogFile = fullfile(outputPath,'recon_logs.mat');
end


