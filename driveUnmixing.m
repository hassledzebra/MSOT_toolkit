% driveUnmixing performs the unmixing process on a reconstructed MSOT dataset.
% 
function [completedTask] = driveUnmixing(inputPath,outputPath)
    

    
    
    % Format output style and such.
    format compact;
    % Deactivate the warning for invoking struct() on an object. 
    warning('off','MATLAB:structOnObject')
    
    
    
    pipelineMeta = load(fullfile(inputPath,'pipeline_settings.mat'));
    assert(numel(fieldnames(pipelineMeta))==1); % Only allow a single variable in the file.
    if isfield(pipelineMeta,'pipelineMeta')
        pipelineMeta = pipelineMeta.pipelineMeta;
    end
    
    DEBUG = pipelineMeta.debugSettings.doDebug;
    if DEBUG
       profile off; 
       profile on; 
    end
    
    Ny = pipelineMeta.reconSettings.N_y ;
    Nx = pipelineMeta.reconSettings.N_x;
    Nframes = numel(pipelineMeta.reconSettings.FramesToReconstruct);
        unmixingArgs = pipelineMeta.unmixSettings;
        Ncomps = numel(unmixingArgs.EndmemberNames);
    
    disp('Initiating file link');
    % If everything checks out, create the file links.
    dataLoader = filter.BinLoader(inputPath,...
        'recon_images',...
        [ Ny, Nx],...
        'double',...
        Nframes);
    
    

    dataWriter = filter.BinWriter(outputPath,...
        'unmix_images',...
        [Ny,Nx,Ncomps],...
        'double',...
        Nframes);
    

    
    
    % instantiate the unmixing system.
    disp('Setting up unmixing system');
    unmixer = unmix.UnmixSystem(unmixingArgs);

    
    disp('Running unmixing');
    for k = 1:Nframes
        k
        recFrame = dataLoader(k);
        %% Construct DataFrame for input into the UnmixSystem.
        inputImage.Data = recFrame.Data;
        inputImage.Meta = recFrame.Meta.header;
        
        [testOutput,testLog] = unmixer(inputImage);
        
        dataWriter(testOutput,k);
    end

    disp('Completed unmixing');
    %% Save to the local work directory, or the original host if being called from not-nextflow.
    
    
    whos
    
    disp('Starting file writing');
    if DEBUG
        unmixProfileInfo = profile('info');
        %profile viewer;
        try
            save(fullfile(outputPath,'unmix_debug.mat'),'clientCollectOperation','-v7.3');
        catch
            save(fullfile(outputPath,'unmix_debug.mat'),'unmixProfileInfo','-v7.3');
        end
    end
    
    save(fullfile(outputPath,'unmix_pipeline_settings.mat'),'pipelineMeta','-v7.3');
    
    completedTask = fullfile(outputPath)
    
    disp('Finished file writing');
    
    
    pause(2); % Give any file writing a fair chance of finishing. 
    
    
end

