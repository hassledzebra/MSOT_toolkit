clear all; close all; clc;


% myFilt = filter.state.AlphaFilter; 
% testFrame.Data = 1;
% myFilt(testFrame)
% 
% t = 1:100;
% y = sin(t/10);
% 
% for k = 1:100
%     testFrame.Data = y(k);
%    output(k) = myFilt(testFrame);
%     
% end
% 
% 
% figure;
% plot(y);
% hold on;
% plot([output.Data]);
% legend('Orig','Alpha');

util.addJavaLibraries

msotFile = '/project/apps_database/MSOT/data/test/phantom/IPASC phantom scan/Scan_3/Scan_3.msot';

% Scan interfacing.
loader = util.MSOTSignalLoader(msotFile);
meta = loader.Meta;

pipelineMeta = util.settingsFromMeta(util.loadDefault('PipelineSettings'),loader.Meta);

pipelineMeta.unmixSettings.MultispectralStateFilterType = 'kalata';
pipelineMeta.unmixSettings = importSpectra(pipelineMeta.unmixSettings);
pipelineMeta.reconSettings.SpeedOfSound = 1520;

reconstructor = recon.ReconSystem(pipelineMeta.reconSettings);
prefilt = recon.MSOTPreFilter(pipelineMeta.filterSettings);


unmixer = unmix.UnmixSystem(pipelineMeta.unmixSettings);
for k = 1:100
frame = loader(k);

filtFrame = prefilt(frame);

recFrame = reconstructor(filtFrame);

testFrame.Data = recFrame;
testFrame.Meta = frame.Meta;

unmixFrame = unmixer(testFrame);

imagesc(reshape(unmixFrame.Data,[128 256])); pause(0.1);
end



function unmixSettings = importSpectra(unmixSettings)
    
    locComponents = unmixSettings.EndmemberNames;
    nComp = numel(unmixSettings.EndmemberNames);
    locWL = unmixSettings.UnmixingWavelengths;
    
    USER_SPECTRA_FOLDER = util.loadDefault('USER_SPECTRA_FOLDER');
    SYSTEM_SPECTRA_FOLDER = util.loadDefault('SYSTEM_SPECTRA_FOLDER');
    for k = 1 : nComp
                try
                    specOptions = cellfun(@(x) fullfile(USER_SPECTRA_FOLDER,[locComponents{k},'Spectrum.',x]),{'mat','csv','xls','xlsx'},'uni',false);
                    specFileExists = isfile(specOptions);
                    
                    if any(specFileExists)
                        specIndex = find(specFileExists,1,'first');
                        [locSpecWL,locSpecSpectrum,~] = unmix.loadAndCheckSpectrum(specOptions{specIndex});
                    else
                        
                        specOptions = cellfun(@(x) fullfile(SYSTEM_SPECTRA_FOLDER,[locComponents{k},'Spectrum.',x]),{'mat','csv','xls','xlsx'},'uni',false);
                        
                        specFileExists = isfile(specOptions);
                        if any(specFileExists)
                            
                        specIndex = find(specFileExists,1,'first');
                        [locSpecWL,locSpecSpectrum,~] = unmix.loadAndCheckSpectrum(specOptions{specIndex});
                        else
                            ME = MException('PreprocessingUnmix:NoSpectrumFile','No spectrum file found');
                            throw(ME);
                        end
                    end
                catch 
                    warning(['Unable to load ' locComponents{k} ' spectrum from user or system spectra folders. Skipping.']);
                end
                
                % Interpolate into the form we need.
                unmixSettings.UnmixingSpectrum{k} = interp1(locSpecWL,locSpecSpectrum,locWL);
                unmixSettings.EndmemberWavelengths{k} = locSpecWL;
                unmixSettings.EndmemberSpectrum{k}    = locSpecSpectrum;
            end
    
end

