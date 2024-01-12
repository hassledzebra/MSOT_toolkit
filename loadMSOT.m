
function loadMSOT(msotpath, path, modeltype)

%path=uigetdir(pwd);
%scanlist = [1,2];

% Input parameters.
% outputfolder = fullfile('C:\Users\McNally Group\Documents\breast_msot_data_clinical\MSOT_Data\Study_23\Scan_2','output');
outputfolder = fullfile(path,'output')
mkdir(outputfolder)
outputDirectory = outputfolder;
%msotFile = 'C:\Users\McNally Group\Documents\breast_msot_data_clinical\MSOT_Data\Study_23\Scan_2\Scan_2.msot';
sprintf('working on');
msotFile = msotpath

%% Preprocessing
util.addJavaLibraries;
warning('off','MATLAB:structOnObject')

% Scan interfacing.
loader = util.MSOTSignalLoader(msotFile);
meta = loader.Meta;
%%
% Populating settings from defaults.
pipelineMeta = util.settingsFromMeta(util.loadDefault('PipelineSettings'),meta);
%pipelineMeta.reconSettings.ModelType = 'backproject2D';
%pipelineMeta.reconSettings.ModelType = 'CDMMI';
pipelineMeta.reconSettings.ModelType = modeltype;

% Tuning speed of sound
pipelineMeta.reconSettings.SpeedOfSound = recon.backProjectTune(loader,meta,...
                                            max(pipelineMeta.reconSettings.N_x,300) ,...
                                            pipelineMeta.reconSettings.FieldOfView_X);
                                        
pipelineMeta.reconSettings.N_y = 256;
pipelineMeta.reconSettings.N_x = 256;


%% Constructors and initialization of processing objects                                    
% Creating prefilter which corrects for the impulse response function
preFilter = recon.MSOTPreFilter(pipelineMeta.filterSettings);

% Creating reconstructor.
reconstructor = recon.ReconSystem(pipelineMeta.reconSettings);

% Creating unmixer.
pipelineMeta.unmixSettings = util.importSpectra(pipelineMeta.unmixSettings);
unmixingArgs = pipelineMeta.unmixSettings;
unmixer = unmix.UnmixSystem(unmixingArgs);
Ncomps = numel(unmixingArgs.EndmemberNames);


% Defining reconstruction set. In this case it is every available scan frame.
% NFrames = numel(loader.Meta.ScanFrames);
% reconstructionSet = 1:NFrames;

reconstructionSet = 1:size(loader.LaserEnergy,2);
NFrames = numel(reconstructionSet);
NRep = meta.RepNum;
reconImageStructure = zeros( pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x , NFrames);
unmixImageStructure = zeros( pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x * Ncomps , NFrames);
reconLogStructure = cell(NFrames,1);


% % Creating writers
% reconWriter = filter.BinWriter(...
%                 outputDirectory,...
%                 'recon_images',...
%                 [pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x],...
%                 'double',...
%                 NFrames);
%             
% unmixWriter = filter.BinWriter(outputDirectory,...
%                 'unmix_images',...
%                 [pipelineMeta.reconSettings.N_y,pipelineMeta.reconSettings.N_x,Ncomps],...
%                 'double',...
%                 NFrames);
        
wavelength = zeros(NFrames,1);

% Serial reconstruction and unmixing
tic
h=waitbar(0, '');
for k = reconstructionSet
    dataFrame = loader(k,{});
    wavelength(k) = dataFrame.Meta.Wavelength;
    filteredSignal = preFilter(dataFrame);

    [reconVec,reconMessage] = reconstructor(filteredSignal);
    reconLog.message = reconMessage;
    reconLog.header = struct(dataFrame.Meta);
    
    recFrame.Data = reconVec;
    recFrame.Meta = reconLog;
    
    reconImageStructure(:,:,k) = reconVec;
    
    % reconWriter(recFrame,k);
    
    inputImage.Data = recFrame.Data;
    inputImage.Meta = recFrame.Meta.header;
    
    [unmixFrame,unmixLog] = unmixer(inputImage);
    
    % unmixWriter(unmixFrame,k);
    
    unmixImageStructure(:,:,k) = reshape(unmixFrame.Data,[pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x * Ncomps]);
    
%     figure(1);
%     imagesc(reconVec); title(num2str(k)); drawnow;
%     figure(2);
%     imagesc(reshape(unmixFrame.Data,[pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x * Ncomps]));
    
    waitbar(k/NFrames, h, sprintf('recon progress =  %d%%',round(k/NFrames*100)))
end
toc
delete(h)
%%
cd(outputDirectory)
% recon shape: 128x128xWavelengthxRepetition
recon = reshape(reconImageStructure,[size(reconImageStructure,1),size(reconImageStructure,2),size(reconImageStructure,3)/NRep, NRep]);
unmixed = reshape(unmixImageStructure,[size(reconImageStructure,1),size(unmixImageStructure,2),size(unmixImageStructure,3)/NRep, NRep]);
% recon average
recon_average = mean(recon, 4);
unmixed_average = mean(unmixed, 4);
w = wavelength(1:size(wavelength,1)/NRep);
save('msot_processed', 'recon','w')
save('unmix_processed', 'unmixed')
save('dataframe','dataFrame')

% slice = squeeze(recon(:,:,7,15));
% slice = medfilt2(imresize(slice,5,'nearest'),[10 10]);
% figure
% imagesc(slice,[0 0.02])
% colormap('jet')
% colorbar
% axis off
% 
% slice_unmix = unmixed(:,:,7, 15);
% slice_unmix = medfilt2(imresize(slice_unmix,[512, 512]),[3 3]);
% figure
% imagesc(slice_unmix, [0 0.0001] )
% colormap('jet')
% colorbar
% axis off


% figure
% plot(w, squeeze(recon_average(100, 100, :)))
end


