function unmixSettings = importSpectra(unmixSettings)
    
    locComponents = unmixSettings.EndmemberNames;
    nComp = numel(unmixSettings.EndmemberNames);
    locWL = unmixSettings.UnmixingWavelengths;
    
    USER_SPECTRA_FOLDER = util.loadDefault('USER_SPECTRA_FOLDER')
    SYSTEM_SPECTRA_FOLDER = util.loadDefault('SYSTEM_SPECTRA_FOLDER')
    for k = 1 : nComp
                try
                    specOptions = cellfun(@(x) fullfile(USER_SPECTRA_FOLDER,[locComponents{k},'Spectrum.',x]),{'mat','csv','xls','xlsx'},'uni',false);
                    %disp(specOptions)
                    specFileExists = isfile(specOptions);
                    
                    if any(specFileExists)
                        specIndex = find(specFileExists,1,'first');
                        [locSpecWL,locSpecSpectrum,~] = unmix.loadAndCheckSpectrum(specOptions{specIndex});
                    else
                        
                        specOptions = cellfun(@(x) fullfile(SYSTEM_SPECTRA_FOLDER,[locComponents{k},'Spectrum.',x]),{'mat','csv','xls','xlsx'},'uni',false);
                        %disp(specOptions)
                        specFileExists = isfile(specOptions);
                        if any(specFileExists)
                            %disp('found')
                        specIndex = find(specFileExists,1,'first');
                        %disp(specOptions{specIndex})
                        [locSpecWL,locSpecSpectrum,~] = unmix.loadAndCheckSpectrum(specOptions{specIndex});
                        else
                            ME = MException('PreprocessingUnmix:NoSpectrumFile','No spectrum file found');
                            throw(ME);
                        end
                    end
                catch 
                    warning(['Unable to load ' locComponents{k} ' spectrum from user folder ', USER_SPECTRA_FOLDER,' or system spectra folder ', SYSTEM_SPECTRA_FOLDER,'. Skipping.']);
%                     warning(fprintf('User spectra folder is %s',USER_SPECTRA_FOLDER));
%                     warning(fprintf('System spectra folder is %s',SYSTEM_SPECTRA_FOLDER));
                end
                
                % Interpolate into the form we need.
                unmixSettings.UnmixingSpectrum{k} = interp1(locSpecWL,locSpecSpectrum,locWL);
                unmixSettings.EndmemberWavelengths{k} = locSpecWL;
                unmixSettings.EndmemberSpectrum{k}    = locSpecSpectrum;
            end
    
end