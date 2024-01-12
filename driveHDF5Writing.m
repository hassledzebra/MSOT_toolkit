function [completedTask] = driveHDF5Writing(inputPath,outputPath)
    %driveHDF5Writing converts reconstructed and unmixed images, stored as
    %binary files with associated metadata, into corresponding HDF5 / H5
    %files. The InputPath should contain the collated reconstructions
    %and/or the unmixed images; the outputPath should point to a valid path
    %on the invoking system. 
    %
    % It is assumed that the outputPath does not have 'unmixed.h5' or
    % 'recons.h5' in it, otherwise this function will fail. 
    
    completedTask = 0; % Sentinel
    
    % Load the pipeline settings for the data structure.
    pipelineSettingsFile = fullfile(inputPath,'pipeline_settings.mat');
    metadataFile = fullfile(inputPath,'recon_metadata.mat');
    
    %% Get the links into each of the data repos.
    pipelineSettings = load(pipelineSettingsFile);
    pipelineSettings = pipelineSettings.pipelineMeta;
    meta = load(metadataFile);
    
    DEBUG = pipelineSettings.debugSettings.doDebug;
    
    if DEBUG && ~isdeployed
        profile off;
        profile on;
    end
    
    % Get the sizes for each of the output files.
    Ny = pipelineSettings.reconSettings.N_y ;
    Nx = pipelineSettings.reconSettings.N_x;
    Nframes = numel(pipelineSettings.reconSettings.FramesToReconstruct);
    unmixingArgs = pipelineSettings.unmixSettings;
    Ncomps = numel(unmixingArgs.EndmemberNames);
    
    % Loaders for each of recons and unmixed.
    reconLoader = filter.BinLoader(inputPath,...
        'recon_images',...
        [Ny Nx],...
        'double',...
        Nframes);
    
    unmixLoader = filter.BinLoader(inputPath,...
        'unmix_images',...
        [Ny Nx Ncomps],...
        'double',...
        Nframes);
    
    % Get data the HDF5 writer needs for invocation.
    reconVarSize = [Ny Nx Nframes];
    reconVarClass = reconLoader.FrameFormat;
    
    unmixVarClass = unmixLoader.FrameFormat;
    
    
    scanFrameInfo = meta.metaStruct.ScanFrames;
    metaFieldnames = fieldnames(scanFrameInfo);
    metaFieldvalues = struct2cell(scanFrameInfo);
    
    reconInfo = pipelineSettings.reconSettings;
    reconFieldnames = fieldnames(reconInfo);
    reconFieldvalues = struct2cell(reconInfo);
    
    unmixInfo = pipelineSettings.unmixSettings;
    
    unmixFieldnames = fieldnames(unmixInfo);
    unmixFieldvalues = struct2cell(unmixInfo);
    
    %% Write out the reconstructions.
    
    reconFid = h5filecreate(fullfile(outputPath,'recons.h5'),'userblock',1024);
    
        % Attach the original frame metadata arrays.
            for k = 1:numel(metaFieldnames)

                storeName = ['FRAME_META/',metaFieldnames{k}];
                h5datacreate(reconFid,storeName,'size',size(scanFrameInfo),'type',class(metaFieldvalues{k}));

                h5varput(reconFid,storeName,[scanFrameInfo.(metaFieldnames{k})]);
            end
    
        % Create the recon dataset.
            storeName = ['/RECON'];
            h5datacreate(reconFid,storeName,'type',reconVarClass,'size',reconVarSize,'max_size',reconVarSize,'chunk_size',[Ny Nx 1]);
    
        % Populate the dataset.
            for k = 1:reconVarSize(3)
                recFrame = reconLoader(k);
                h5write(fullfile(outputPath,'recons.h5'),'/RECON',recFrame.Data,[1 1 k],[Ny Nx 1],[1 1 1]);
            end
    
    % Attach the processing description to the dataset.
    for k = 1:numel(reconFieldnames)
        storeName = '/RECON';
        if strcmp(reconFieldnames{k},'FramesToReconstruct') % TODO: This is done because of a complaint about a too-large header object. Is this a result of the 'userblock' being a particular size?
            h5attput(fullfile(outputPath,'recons.h5'),storeName,'NFramesReconstructed',numel(reconFieldvalues{k}));
        else
            h5attput(fullfile(outputPath,'recons.h5'),storeName,reconFieldnames{k},reconFieldvalues{k});
        end
    end
    
    %% Write out the unmixed images. Give each endmember its own path so that users can more easily access the time-courses.
    unmixFid = h5filecreate(fullfile(outputPath,'unmixed.h5'),'userblock',1024);
    names = unmixInfo.EndmemberNames;
    
    for p = 1:numel(names)
        storeName = ['/UNMIX/',names{p}];
        h5datacreate(unmixFid,storeName,'type',unmixVarClass,'size',[Ny Nx Nframes],'max_size',[Ny Nx Nframes],'chunk_size',[Ny Nx 1]);
        endmemberOffset = p;
        for k = 1:Nframes
            imageIndex = (k-1).*numel(names) + endmemberOffset;
            unmixFrame = unmixLoader(k).Data;
            h5write(fullfile(outputPath,'unmixed.h5'),storeName,unmixFrame(:,:,endmemberOffset),[1 1 k],[Ny Nx 1],[1 1 1]);
        end
        
    end
    
    
    % Debug stuff. 
    if DEBUG && ~isdeployed
        hdf5ProfileInfo = profile('info');
        try
            save(fullfile(outputPath,'hdf5_debug.mat'),'clientCollectOperation','-v7.3');
        catch
            save(fullfile(outputPath,'hdf5_debug.mat'),'hdf5ProfileInfo','-v7.3');
        end
    end
    
    
    completedTask = 1;
    
    
end

