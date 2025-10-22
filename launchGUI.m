function launchGUI()
% launchGUI - Launch the MSOT Pipeline GUI
%
% Usage:
%   launchGUI()
%
% This function starts the MSOT Pipeline GUI for processing
% multispectral optoacoustic tomography data.
%
% The GUI provides an intuitive interface for:
%   - Loading MSOT data files
%   - Configuring reconstruction settings
%   - Setting up spectral unmixing parameters
%   - Processing data with live visualization
%   - Saving results
%
% See also: MSOTPipelineGUI, examplePipeline

    % Check for required toolboxes
    v = ver;
    hasImageProcessing = any(strcmp({v.Name}, 'Image Processing Toolbox'));
    hasSignalProcessing = any(strcmp({v.Name}, 'Signal Processing Toolbox'));

    if ~hasImageProcessing
        warning('Image Processing Toolbox not detected. Some features may not work.');
    end

    if ~hasSignalProcessing
        warning('Signal Processing Toolbox not detected. Some features may not work.');
    end

    % Launch the GUI
    fprintf('\n========================================\n');
    fprintf('  MSOT Pipeline GUI\n');
    fprintf('========================================\n');
    fprintf('Starting graphical user interface...\n\n');

    try
        MSOTPipelineGUI;
    catch ME
        fprintf('Error launching GUI: %s\n', ME.message);
        fprintf('Stack trace:\n');
        for i = 1:length(ME.stack)
            fprintf('  File: %s\n', ME.stack(i).file);
            fprintf('  Line: %d, Function: %s\n\n', ME.stack(i).line, ME.stack(i).name);
        end
    end
end
