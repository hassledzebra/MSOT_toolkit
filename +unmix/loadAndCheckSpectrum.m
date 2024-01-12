function [wavelengths,molarExtinctionCoeff,spectrumName] = loadAndCheckSpectrum(filename)
    %loadAndCheckSpectrum loads a file containing an absorption spectrum
    %for an endmember of interest. 
    %
    % Input: A path to a file, any of .mat, .xls, or .csv, which contains
    % a spectrum of an endmember of interest. The filename must end in
    % 'spectrum' or 'Spectrum', with all characters before taken as the
    % name of the endmember. The spectrum is required to have two columns
    % in the variable: The first contains the sampled wavelengths in
    % nanometers (nm, 10^-9 m) and the second contains the absorption
    % coefficients in cm^-1/mol.
    
    
    [fPath,fName,fExt] = fileparts(filename);
    
    
    switch fExt
        case '.mat'
            [wavelengths,molarExtinctionCoeff] = parseMatFile(filename);
        case '.xls'
            [wavelengths,molarExtinctionCoeff] = parseXLSFile(filename);
        case '.csv'
            [wavelengths,molarExtinctionCoeff] = parseCSVFile(filename);
        otherwise
            error(['Invalid file extension for spectrum file' filename]);
    end
    
    % strip off 'spectrum' from the name.
    try
        if strcmpi(fName((end-7):end),'spectrum')
           spectrumName = fName(1:(end-8)); 
        else
           spectrumName = fName;
        end
    catch
       error(['Spectrum file ' filename ' does not conform to naming standard.' ]);
    end
    
    
end


% Reads out the contents of a MAT file.
%   For now we dictate that the mat file should only have one variable inside,
%   with one column signifying wavelength in nanometers and the other the
%   molar extinction coefficient in cm-1/mol
function [wavelengths,molarExtinctionCoeff] = parseMatFile(filename)
    M = load(filename);
    
    Mfields = fieldnames(M);
    
    % Assertions to enforce compliance.
    assert(numel(Mfields)==1);
    
    spectrum = M.(Mfields{1});
    
    assert(size(spectrum,2) == 2);
    assert(isa(spectrum,'double'));
    
    wavelengths = spectrum(:,1);
    molarExtinctionCoeff = spectrum(:,2);
end


function [wavelengths,molarExtinctionCoeff] = parseXLSFile(filename)
    spectrum = xlsread(filename);
    
    % Assertions to enforce compliance.
    assert(isa(spectrum,'double'));
    assert(size(spectrum,2)==2);
    
    % Strip out any NaNs
    if any(isnan(spectrum(:)))
        spectrum = spectrum(sum(isnan(spectrum),2)==0,:);
        warning('Removed all rows that imported as NaN');
    end
    
    
    wavelengths = spectrum(:,1);
    molarExtinctionCoeff = spectrum(:,2);
end


function [wavelengths,molarExtinctionCoeff] = parseCSVFile(filename)
    try
    spectrum = csvread(filename);
    catch
        try
            spectrum = csvread(filename,1);
        catch
            error(['Error importing ' filename '. There may be too many header lines.']);
        end
    end
    
    % Assertions to enforce compliance.
    assert(isa(spectrum,'double'));
    assert(size(spectrum,2)==2);
    
    
    
    wavelengths = spectrum(:,1);
    molarExtinctionCoeff = spectrum(:,2);
end






