function [mixingModel] = linearMixingModel(inputCoords,outputCoords,varargin)
    %linearMixingModel builds a matrix representation between the
    %endmembers and the wavelength spaces.
    doKron = false;
    if numel(varargin)>0
    if strncmpi(varargin{1},'kron',4)
        doKron = true;
    end
    end
    
    mixSet.EndmemberNames = util.fixEndmembers(inputCoords.getCoordinateValues('Endmembers'));
    nComp = numel(mixSet.EndmemberNames);
    mixSet.UnmixingWavelengths = outputCoords.getCoordinateValues('Wavelengths');
    
    mixSet = util.importSpectra(mixSet);
    
    mixingModel = zeros(numel(mixSet.UnmixingWavelengths),numel(mixSet.EndmemberNames));
    
    for k = 1 : nComp
        mixingModel(:,k) = mixSet.UnmixingSpectrum{k};
    end
   
    if doKron
        mixingModel = kron(mixingModel,varargin{2});
    end
    
end



