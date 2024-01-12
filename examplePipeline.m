

clearvars; close all; clc;

% Input parameters.
outputfolder = fullfile('E:\ZH 2-DG 2-LMP\Scan_7','output');
mkdir(outputfolder)
outputDirectory = outputfolder;
msotFile = 'E:\ZH 2-DG 2-LMP\Scan_7\Scan_7.msot';

%% Preprocessing
util.addJavaLibraries;
warning('off','MATLAB:structOnObject')

% Scan interfacing.
loader = util.MSOTSignalLoader(msotFile);
meta = loader.Meta;
%%
% Populating settings from defaults.
pipelineMeta = util.settingsFromMeta(util.loadDefault('PipelineSettings'),meta);
pipelineMeta.reconSettings.ModelType = 'backproject2D';
% pipelineMeta.reconSettings.ModelType = 'CDMMI';

% pipelineMeta.parallelSettings.parallelType = 'multiple';


% Tuning speed of sound
pipelineMeta.reconSettings.SpeedOfSound = recon.backProjectTune(loader,meta,...
                                            max(pipelineMeta.reconSettings.N_x,300) ,...
                                            pipelineMeta.reconSettings.FieldOfView_X);
                                        
pipelineMeta.reconSettings.N_y = 128;
pipelineMeta.reconSettings.N_x = 128;


pipelineMeta.unmixSettings.EndmemberNames = {'Deoxyhemoglobin','Oxyhemoglobin','IR800CW'};
% pipelineMeta.unmixSettings.EndmemberNames = {'Deoxyhemoglobin','Oxyhemoglobin'};
% pipelineMeta.unmixSettings.UnmixSolverType = 'nonNegAPCG';

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
unmixImageStructure = zeros( pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x * Ncomps , NRep);
reconLogStructure = cell(NFrames,1);
     
wavelength = zeros(NFrames,1);

% Serial reconstruction and unmixing
tic
h=waitbar(0, '');
numIterations = max(reconstructionSet);

% Then construct a ParforProgressbar object:
% parfor_progress(numIterations);

%      parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress;
%      end
%      parfor_progress(0);


for k = reconstructionSet
    dataFrame = loader(k,{});
    wavelength(k) = dataFrame.Meta.Wavelength;
    filteredSignal = preFilter(dataFrame);

    [reconVec,reconMessage] = reconstructor(filteredSignal);
%     reconLog.message = reconMessage;
%     reconLog.header = struct(dataFrame.Meta);
    
%     recFrame.Data = reconVec;
%     recFrame.Meta = reconLog;
%     
    reconImageStructure(:,:,k) = reconVec;
    
    % reconWriter(recFrame,k);
    
    inputImage.Data = reconVec;
    inputImage.Meta = struct(dataFrame.Meta);
    
    [unmixFrame,unmixLog] = unmixer(inputImage);
    
    % unmixWriter(unmixFrame,k);
    
    unmixImageStructure(:,:,k) = reshape(unmixFrame.Data,[pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x * Ncomps]);
    
    figure(1);
    imagesc(reconVec); title(num2str(k)); drawnow;
%     parfor_progress
    figure(2);
    imagesc(reshape(unmixFrame.Data,[pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x * Ncomps]));
    
    waitbar(k/NFrames, h, sprintf('progress =  %d%%',round(k/NFrames*100)))
end
toc
delete(h)
%%
% Num_w = NFrames/NRep;
% for rep = 1:NRep
%     inputImage.Data = reconImageStructure(:,:,rep+1 : rep*Num_w);
%     [unmixFrame,unmixLog] = unmixer(inputImage);
%     unmixImageStructure(:,:,Rep) = reshape(unmixFrame.Data,[pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x * Ncomps]);
%     figure(2);
%     imagesc(reshape(unmixFrame.Data,[pipelineMeta.reconSettings.N_y , pipelineMeta.reconSettings.N_x * Ncomps]));
% end

% parfor_progress(0);
%%

% recon shape: 128x128xWavelengthxRepetition
recon_img = reshape(reconImageStructure,[size(reconImageStructure,1),size(reconImageStructure,2),size(reconImageStructure,3)/NRep, NRep]);
unmix_img = reshape(unmixImageStructure,[size(unmixImageStructure,1),size(unmixImageStructure,2),size(unmixImageStructure,3)/NRep, NRep]);
% recon average
% recon_average = mean(recon_img, 4);
w = wavelength(1:size(wavelength,1)/NRep);
save(fullfile(outputfolder,'msot_processed.mat'), 'recon_img','w')
save(fullfile(outputfolder,'unmix_processed.mat'), 'unmix_img','w')

slice = reconImageStructure(:,:,1);
% slice = medfilt2(imresize(slice,[512, 512]),[2 2]);
figure
imagesc(slice)
colormap('jet')
colorbar


% figure
% plot(w, squeeze(recon_average(100, 100, :)))


