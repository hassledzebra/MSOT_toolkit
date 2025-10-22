classdef MSOTPipelineGUI < matlab.apps.AppBase
    % MSOTPipelineGUI - GUI for MSOT Image Processing Pipeline
    % Based on the workflow in examplePipeline.m

    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup

        % Tab 1: Setup
        SetupTab                        matlab.ui.container.Tab
        InputFileLabel                  matlab.ui.control.Label
        InputFileEditField              matlab.ui.control.EditField
        BrowseInputButton               matlab.ui.control.Button
        OutputFolderLabel               matlab.ui.control.Label
        OutputFolderEditField           matlab.ui.control.EditField
        BrowseOutputButton              matlab.ui.control.Button
        LoadDataButton                  matlab.ui.control.Button
        MetadataTextArea                matlab.ui.control.TextArea
        MetadataLabel                   matlab.ui.control.Label

        % Tab 2: Reconstruction Settings
        ReconTab                        matlab.ui.container.Tab
        ModelTypeLabel                  matlab.ui.control.Label
        ModelTypeDropDown               matlab.ui.control.DropDown
        ResolutionLabel                 matlab.ui.control.Label
        NxEditField                     matlab.ui.control.NumericEditField
        NxLabel                         matlab.ui.control.Label
        NyEditField                     matlab.ui.control.NumericEditField
        NyLabel                         matlab.ui.control.Label
        SpeedOfSoundLabel               matlab.ui.control.Label
        AutoTuneCheckBox                matlab.ui.control.CheckBox
        SpeedOfSoundEditField           matlab.ui.control.NumericEditField
        TuneSOSButton                   matlab.ui.control.Button
        ParallelCheckBox                matlab.ui.control.CheckBox

        % Tab 3: Unmixing Settings
        UnmixTab                        matlab.ui.container.Tab
        EndmemberLabel                  matlab.ui.control.Label
        HbCheckBox                      matlab.ui.control.CheckBox
        HbO2CheckBox                    matlab.ui.control.CheckBox
        WaterCheckBox                   matlab.ui.control.CheckBox
        FatCheckBox                     matlab.ui.control.CheckBox
        IR800CWCheckBox                 matlab.ui.control.CheckBox
        ICGCheckBox                     matlab.ui.control.CheckBox
        CustomEndmemberEditField        matlab.ui.control.EditField
        CustomEndmemberLabel            matlab.ui.control.Label
        UnmixSolverLabel                matlab.ui.control.Label
        UnmixSolverDropDown             matlab.ui.control.DropDown
        StateFilterLabel                matlab.ui.control.Label
        StateFilterDropDown             matlab.ui.control.DropDown

        % Tab 4: Process & View
        ProcessTab                      matlab.ui.container.Tab
        RunButton                       matlab.ui.control.Button
        StopButton                      matlab.ui.control.Button
        ProgressGauge                   matlab.ui.control.LinearGauge
        ProgressLabel                   matlab.ui.control.Label
        StatusTextArea                  matlab.ui.control.TextArea
        ReconAxes                       matlab.ui.control.UIAxes
        UnmixAxes                       matlab.ui.control.UIAxes
        SaveResultsButton               matlab.ui.control.Button
        FrameSlider                     matlab.ui.control.Slider
        FrameSliderLabel                matlab.ui.control.Label
        CurrentFrameLabel               matlab.ui.control.Label
    end

    properties (Access = private)
        % Data properties
        loader                          % MSOTSignalLoader object
        meta                            % Metadata structure
        pipelineMeta                    % Pipeline settings
        preFilter                       % MSOTPreFilter object
        reconstructor                   % ReconSystem object
        unmixer                         % UnmixSystem object

        % Results
        reconImageStructure             % Reconstructed images
        unmixImageStructure             % Unmixed images
        wavelength                      % Wavelength array
        reconstructionSet               % Frames to reconstruct

        % Processing control
        stopProcessing = false          % Flag to stop processing
    end

    methods (Access = private)

        function createComponents(app)
            % Create UIFigure and components
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1200 800];
            app.UIFigure.Name = 'MSOT Pipeline GUI';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.UIFigure);
            app.TabGroup.Position = [10 10 1180 780];

            % ===== TAB 1: SETUP =====
            app.SetupTab = uitab(app.TabGroup);
            app.SetupTab.Title = '1. Setup';

            % Input file selection
            app.InputFileLabel = uilabel(app.SetupTab);
            app.InputFileLabel.Position = [20 700 100 22];
            app.InputFileLabel.Text = 'MSOT File:';

            app.InputFileEditField = uieditfield(app.SetupTab, 'text');
            app.InputFileEditField.Position = [130 700 850 22];

            app.BrowseInputButton = uibutton(app.SetupTab, 'push');
            app.BrowseInputButton.Position = [1000 700 100 22];
            app.BrowseInputButton.Text = 'Browse...';
            app.BrowseInputButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseInputButtonPushed, true);

            % Output folder selection
            app.OutputFolderLabel = uilabel(app.SetupTab);
            app.OutputFolderLabel.Position = [20 660 100 22];
            app.OutputFolderLabel.Text = 'Output Folder:';

            app.OutputFolderEditField = uieditfield(app.SetupTab, 'text');
            app.OutputFolderEditField.Position = [130 660 850 22];

            app.BrowseOutputButton = uibutton(app.SetupTab, 'push');
            app.BrowseOutputButton.Position = [1000 660 100 22];
            app.BrowseOutputButton.Text = 'Browse...';
            app.BrowseOutputButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseOutputButtonPushed, true);

            % Load Data button
            app.LoadDataButton = uibutton(app.SetupTab, 'push');
            app.LoadDataButton.Position = [20 610 150 30];
            app.LoadDataButton.Text = 'Load Data';
            app.LoadDataButton.FontSize = 14;
            app.LoadDataButton.FontWeight = 'bold';
            app.LoadDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadDataButtonPushed, true);

            % Metadata display
            app.MetadataLabel = uilabel(app.SetupTab);
            app.MetadataLabel.Position = [20 570 100 22];
            app.MetadataLabel.Text = 'Scan Metadata:';
            app.MetadataLabel.FontWeight = 'bold';

            app.MetadataTextArea = uitextarea(app.SetupTab);
            app.MetadataTextArea.Position = [20 50 1120 510];
            app.MetadataTextArea.Editable = 'off';

            % ===== TAB 2: RECONSTRUCTION SETTINGS =====
            app.ReconTab = uitab(app.TabGroup);
            app.ReconTab.Title = '2. Reconstruction';

            % Model Type
            app.ModelTypeLabel = uilabel(app.ReconTab);
            app.ModelTypeLabel.Position = [20 680 120 22];
            app.ModelTypeLabel.Text = 'Model Type:';
            app.ModelTypeLabel.FontWeight = 'bold';

            app.ModelTypeDropDown = uidropdown(app.ReconTab);
            app.ModelTypeDropDown.Items = {'backproject2D', 'CDMMI', 'CDMMI2', 'dIMMI'};
            app.ModelTypeDropDown.Value = 'backproject2D';
            app.ModelTypeDropDown.Position = [150 680 200 22];

            % Resolution Settings
            app.ResolutionLabel = uilabel(app.ReconTab);
            app.ResolutionLabel.Position = [20 630 120 22];
            app.ResolutionLabel.Text = 'Image Resolution:';
            app.ResolutionLabel.FontWeight = 'bold';

            app.NxLabel = uilabel(app.ReconTab);
            app.NxLabel.Position = [150 630 30 22];
            app.NxLabel.Text = 'N_x:';

            app.NxEditField = uieditfield(app.ReconTab, 'numeric');
            app.NxEditField.Position = [190 630 80 22];
            app.NxEditField.Value = 128;
            app.NxEditField.Limits = [64 512];

            app.NyLabel = uilabel(app.ReconTab);
            app.NyLabel.Position = [290 630 30 22];
            app.NyLabel.Text = 'N_y:';

            app.NyEditField = uieditfield(app.ReconTab, 'numeric');
            app.NyEditField.Position = [330 630 80 22];
            app.NyEditField.Value = 128;
            app.NyEditField.Limits = [64 512];

            % Speed of Sound
            app.SpeedOfSoundLabel = uilabel(app.ReconTab);
            app.SpeedOfSoundLabel.Position = [20 580 150 22];
            app.SpeedOfSoundLabel.Text = 'Speed of Sound (m/s):';
            app.SpeedOfSoundLabel.FontWeight = 'bold';

            app.AutoTuneCheckBox = uicheckbox(app.ReconTab);
            app.AutoTuneCheckBox.Position = [180 580 100 22];
            app.AutoTuneCheckBox.Text = 'Auto-Tune';
            app.AutoTuneCheckBox.Value = true;
            app.AutoTuneCheckBox.ValueChangedFcn = createCallbackFcn(app, @AutoTuneCheckBoxValueChanged, true);

            app.SpeedOfSoundEditField = uieditfield(app.ReconTab, 'numeric');
            app.SpeedOfSoundEditField.Position = [290 580 100 22];
            app.SpeedOfSoundEditField.Value = 1520;
            app.SpeedOfSoundEditField.Enable = 'off';

            app.TuneSOSButton = uibutton(app.ReconTab, 'push');
            app.TuneSOSButton.Position = [410 580 120 22];
            app.TuneSOSButton.Text = 'Tune Now';
            app.TuneSOSButton.ButtonPushedFcn = createCallbackFcn(app, @TuneSOSButtonPushed, true);

            % Parallel Processing
            app.ParallelCheckBox = uicheckbox(app.ReconTab);
            app.ParallelCheckBox.Position = [20 530 200 22];
            app.ParallelCheckBox.Text = 'Enable Parallel Processing';
            app.ParallelCheckBox.Value = false;

            % Info Panel
            infoText = uilabel(app.ReconTab);
            infoText.Position = [20 50 1000 450];
            infoText.Text = sprintf(['Reconstruction Settings:\n\n' ...
                'Model Type:\n' ...
                '  - backproject2D: Fast backprojection (default)\n' ...
                '  - CDMMI: Coordinate Descent MMI\n' ...
                '  - CDMMI2: CDMMI variant 2\n' ...
                '  - dIMMI: Damped Iterative MMI\n\n' ...
                'Image Resolution:\n' ...
                '  - N_x, N_y: Number of pixels (typically 128-250)\n' ...
                '  - Higher resolution = better quality but slower\n\n' ...
                'Speed of Sound:\n' ...
                '  - Auto-Tune: Automatically optimize (recommended)\n' ...
                '  - Manual: Set custom value (default 1520 m/s)\n' ...
                '  - Click "Tune Now" to run optimization with current data\n\n' ...
                'Parallel Processing:\n' ...
                '  - Enable to use multiple CPU cores (requires Parallel Computing Toolbox)']);
            infoText.VerticalAlignment = 'top';

            % ===== TAB 3: UNMIXING SETTINGS =====
            app.UnmixTab = uitab(app.TabGroup);
            app.UnmixTab.Title = '3. Unmixing';

            % Endmember selection
            app.EndmemberLabel = uilabel(app.UnmixTab);
            app.EndmemberLabel.Position = [20 680 200 22];
            app.EndmemberLabel.Text = 'Select Endmembers:';
            app.EndmemberLabel.FontWeight = 'bold';

            yPos = 640;
            app.HbCheckBox = uicheckbox(app.UnmixTab);
            app.HbCheckBox.Position = [40 yPos 200 22];
            app.HbCheckBox.Text = 'Deoxyhemoglobin (Hb)';
            app.HbCheckBox.Value = true;

            yPos = yPos - 30;
            app.HbO2CheckBox = uicheckbox(app.UnmixTab);
            app.HbO2CheckBox.Position = [40 yPos 200 22];
            app.HbO2CheckBox.Text = 'Oxyhemoglobin (HbO2)';
            app.HbO2CheckBox.Value = true;

            yPos = yPos - 30;
            app.WaterCheckBox = uicheckbox(app.UnmixTab);
            app.WaterCheckBox.Position = [40 yPos 200 22];
            app.WaterCheckBox.Text = 'Water';
            app.WaterCheckBox.Value = false;

            yPos = yPos - 30;
            app.FatCheckBox = uicheckbox(app.UnmixTab);
            app.FatCheckBox.Position = [40 yPos 200 22];
            app.FatCheckBox.Text = 'Fat';
            app.FatCheckBox.Value = false;

            yPos = yPos - 30;
            app.IR800CWCheckBox = uicheckbox(app.UnmixTab);
            app.IR800CWCheckBox.Position = [40 yPos 200 22];
            app.IR800CWCheckBox.Text = 'IR800CW';
            app.IR800CWCheckBox.Value = true;

            yPos = yPos - 30;
            app.ICGCheckBox = uicheckbox(app.UnmixTab);
            app.ICGCheckBox.Position = [40 yPos 200 22];
            app.ICGCheckBox.Text = 'ICG (Indocyanine Green)';
            app.ICGCheckBox.Value = false;

            % Custom endmember
            yPos = yPos - 40;
            app.CustomEndmemberLabel = uilabel(app.UnmixTab);
            app.CustomEndmemberLabel.Position = [40 yPos 150 22];
            app.CustomEndmemberLabel.Text = 'Custom Endmembers:';

            yPos = yPos - 30;
            app.CustomEndmemberEditField = uieditfield(app.UnmixTab, 'text');
            app.CustomEndmemberEditField.Position = [40 yPos 400 22];
            app.CustomEndmemberEditField.Placeholder = 'Enter comma-separated names (e.g., MyDye1, MyDye2)';

            % Unmixing solver
            yPos = yPos - 50;
            app.UnmixSolverLabel = uilabel(app.UnmixTab);
            app.UnmixSolverLabel.Position = [40 yPos 120 22];
            app.UnmixSolverLabel.Text = 'Unmixing Solver:';
            app.UnmixSolverLabel.FontWeight = 'bold';

            app.UnmixSolverDropDown = uidropdown(app.UnmixTab);
            app.UnmixSolverDropDown.Items = {'linear', 'nnls', 'nonNegAPCG'};
            app.UnmixSolverDropDown.Value = 'linear';
            app.UnmixSolverDropDown.Position = [170 yPos 200 22];

            % State filter
            yPos = yPos - 40;
            app.StateFilterLabel = uilabel(app.UnmixTab);
            app.StateFilterLabel.Position = [40 yPos 120 22];
            app.StateFilterLabel.Text = 'State Filter Type:';
            app.StateFilterLabel.FontWeight = 'bold';

            app.StateFilterDropDown = uidropdown(app.UnmixTab);
            app.StateFilterDropDown.Items = {'alphabeta', 'alpha', 'kalata', 'slidingFunction'};
            app.StateFilterDropDown.Value = 'alphabeta';
            app.StateFilterDropDown.Position = [170 yPos 200 22];

            % Info text
            infoText2 = uilabel(app.UnmixTab);
            infoText2.Position = [500 300 600 400];
            infoText2.Text = sprintf(['Spectral Unmixing Information:\n\n' ...
                'Endmembers are the pure spectral signatures of biological\n' ...
                'chromophores that will be separated from the multispectral data.\n\n' ...
                'Common Endmembers:\n' ...
                '  - Deoxyhemoglobin (Hb): Deoxygenated blood\n' ...
                '  - Oxyhemoglobin (HbO2): Oxygenated blood\n' ...
                '  - Water: Tissue water content\n' ...
                '  - Fat: Lipid content\n' ...
                '  - IR800CW: Near-infrared fluorescent dye\n' ...
                '  - ICG: Indocyanine green contrast agent\n\n' ...
                'Solver Types:\n' ...
                '  - linear: Standard linear least squares\n' ...
                '  - nnls: Non-negative least squares\n' ...
                '  - nonNegAPCG: Non-negative accelerated projected conjugate gradient\n\n' ...
                'State Filters are used for temporal smoothing of multispectral data.']);
            infoText2.VerticalAlignment = 'top';

            % ===== TAB 4: PROCESS & VIEW =====
            app.ProcessTab = uitab(app.TabGroup);
            app.ProcessTab.Title = '4. Process & View';

            % Run/Stop buttons
            app.RunButton = uibutton(app.ProcessTab, 'push');
            app.RunButton.Position = [20 700 150 40];
            app.RunButton.Text = 'RUN PROCESSING';
            app.RunButton.FontSize = 14;
            app.RunButton.FontWeight = 'bold';
            app.RunButton.BackgroundColor = [0.2 0.8 0.2];
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);

            app.StopButton = uibutton(app.ProcessTab, 'push');
            app.StopButton.Position = [190 700 100 40];
            app.StopButton.Text = 'STOP';
            app.StopButton.FontSize = 14;
            app.StopButton.FontWeight = 'bold';
            app.StopButton.BackgroundColor = [0.9 0.3 0.3];
            app.StopButton.Enable = 'off';
            app.StopButton.ButtonPushedFcn = createCallbackFcn(app, @StopButtonPushed, true);

            % Progress gauge
            app.ProgressLabel = uilabel(app.ProcessTab);
            app.ProgressLabel.Position = [310 720 100 22];
            app.ProgressLabel.Text = 'Progress:';

            app.ProgressGauge = uigauge(app.ProcessTab, 'linear');
            app.ProgressGauge.Position = [310 700 500 20];
            app.ProgressGauge.Limits = [0 100];

            app.CurrentFrameLabel = uilabel(app.ProcessTab);
            app.CurrentFrameLabel.Position = [830 700 250 22];
            app.CurrentFrameLabel.Text = 'Frame: 0/0';

            % Status text area
            app.StatusTextArea = uitextarea(app.ProcessTab);
            app.StatusTextArea.Position = [20 540 500 140];
            app.StatusTextArea.Editable = 'off';
            app.StatusTextArea.Value = {'Ready to process...'};

            % Image axes
            app.ReconAxes = uiaxes(app.ProcessTab);
            app.ReconAxes.Position = [20 50 500 470];
            title(app.ReconAxes, 'Reconstructed Image');

            app.UnmixAxes = uiaxes(app.ProcessTab);
            app.UnmixAxes.Position = [550 50 600 470];
            title(app.UnmixAxes, 'Unmixed Components');

            % Frame slider
            app.FrameSliderLabel = uilabel(app.ProcessTab);
            app.FrameSliderLabel.Position = [550 580 100 22];
            app.FrameSliderLabel.Text = 'View Frame:';

            app.FrameSlider = uislider(app.ProcessTab);
            app.FrameSlider.Position = [650 590 480 3];
            app.FrameSlider.Limits = [1 1];
            app.FrameSlider.Value = 1;
            app.FrameSlider.Enable = 'off';
            app.FrameSlider.ValueChangedFcn = createCallbackFcn(app, @FrameSliderValueChanged, true);

            % Save results button
            app.SaveResultsButton = uibutton(app.ProcessTab, 'push');
            app.SaveResultsButton.Position = [550 630 200 40];
            app.SaveResultsButton.Text = 'Save Results';
            app.SaveResultsButton.FontSize = 14;
            app.SaveResultsButton.Enable = 'off';
            app.SaveResultsButton.ButtonPushedFcn = createCallbackFcn(app, @SaveResultsButtonPushed, true);

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end

        % ===== CALLBACK FUNCTIONS =====

        function BrowseInputButtonPushed(app, ~, ~)
            [file, path] = uigetfile('*.msot', 'Select MSOT File');
            if file ~= 0
                app.InputFileEditField.Value = fullfile(path, file);
            end
        end

        function BrowseOutputButtonPushed(app, ~, ~)
            folder = uigetdir('', 'Select Output Folder');
            if folder ~= 0
                app.OutputFolderEditField.Value = folder;
            end
        end

        function LoadDataButtonPushed(app, ~, ~)
            try
                % Get file path
                msotFile = app.InputFileEditField.Value;
                if isempty(msotFile) || ~isfile(msotFile)
                    uialert(app.UIFigure, 'Please select a valid MSOT file.', 'Error');
                    return;
                end

                % Add Java libraries
                util.addJavaLibraries;
                warning('off','MATLAB:structOnObject');

                % Load data
                app.loader = util.MSOTSignalLoader(msotFile);
                app.meta = app.loader.Meta;

                % Display metadata
                metaStr = sprintf('MSOT Scan Metadata\n');
                metaStr = [metaStr sprintf('==================\n\n')];
                metaStr = [metaStr sprintf('File: %s\n\n', msotFile)];

                if isfield(app.meta, 'ScanName')
                    metaStr = [metaStr sprintf('Scan Name: %s\n', app.meta.ScanName)];
                end
                if isfield(app.meta, 'RepNum')
                    metaStr = [metaStr sprintf('Number of Repetitions: %d\n', app.meta.RepNum)];
                end
                if isfield(app.meta, 'ScanFrames')
                    metaStr = [metaStr sprintf('Total Frames: %d\n', numel(app.meta.ScanFrames))];
                end
                if isfield(app.meta, 'SamplingFrequency')
                    metaStr = [metaStr sprintf('Sampling Frequency: %.2f MHz\n', app.meta.SamplingFrequency/1e6)];
                end
                if isfield(app.meta, 'ZPositions')
                    metaStr = [metaStr sprintf('Z Positions: %d\n', numel(app.meta.ZPositions))];
                end

                % Get wavelengths
                if ~isempty(app.loader.LaserEnergy)
                    nFrames = size(app.loader.LaserEnergy, 2);
                    metaStr = [metaStr sprintf('\nTotal Acquisition Frames: %d\n', nFrames)];

                    % Try to get wavelength info from first frame
                    try
                        dataFrame = app.loader(1, {});
                        if isfield(dataFrame.Meta, 'Wavelength')
                            metaStr = [metaStr sprintf('Wavelength Range: multispectral\n')];
                        end
                    catch
                    end
                end

                app.MetadataTextArea.Value = metaStr;

                % Initialize pipeline settings
                app.pipelineMeta = util.settingsFromMeta(util.loadDefault('PipelineSettings'), app.meta);

                uialert(app.UIFigure, 'Data loaded successfully!', 'Success', 'Icon', 'success');

            catch ME
                uialert(app.UIFigure, sprintf('Error loading data: %s', ME.message), 'Error');
            end
        end

        function AutoTuneCheckBoxValueChanged(app, ~, ~)
            if app.AutoTuneCheckBox.Value
                app.SpeedOfSoundEditField.Enable = 'off';
            else
                app.SpeedOfSoundEditField.Enable = 'on';
            end
        end

        function TuneSOSButtonPushed(app, ~, ~)
            if isempty(app.loader) || isempty(app.meta)
                uialert(app.UIFigure, 'Please load data first.', 'Error');
                return;
            end

            try
                % Get resolution
                Nx = max(app.NxEditField.Value, 300);
                FOV = 0.030; % 30mm default

                % Run tuning
                app.UIFigure.Pointer = 'watch';
                drawnow;

                tunedSOS = recon.backProjectTune(app.loader, app.meta, Nx, FOV);

                app.SpeedOfSoundEditField.Value = tunedSOS;
                app.AutoTuneCheckBox.Value = false;
                app.SpeedOfSoundEditField.Enable = 'on';

                app.UIFigure.Pointer = 'arrow';

                uialert(app.UIFigure, sprintf('Speed of sound tuned to: %.2f m/s', tunedSOS), ...
                    'Tuning Complete', 'Icon', 'success');

            catch ME
                app.UIFigure.Pointer = 'arrow';
                uialert(app.UIFigure, sprintf('Error tuning speed of sound: %s', ME.message), 'Error');
            end
        end

        function RunButtonPushed(app, ~, ~)
            if isempty(app.loader) || isempty(app.meta)
                uialert(app.UIFigure, 'Please load data first.', 'Error');
                return;
            end

            try
                % Disable run button, enable stop button
                app.RunButton.Enable = 'off';
                app.StopButton.Enable = 'on';
                app.stopProcessing = false;

                % Update status
                app.StatusTextArea.Value = {'Starting processing...', ''};
                drawnow;

                % Configure settings from GUI
                app.pipelineMeta.reconSettings.ModelType = app.ModelTypeDropDown.Value;
                app.pipelineMeta.reconSettings.N_x = app.NxEditField.Value;
                app.pipelineMeta.reconSettings.N_y = app.NyEditField.Value;

                % Speed of sound
                if app.AutoTuneCheckBox.Value
                    updateStatus(app, 'Auto-tuning speed of sound...');
                    app.pipelineMeta.reconSettings.SpeedOfSound = recon.backProjectTune(...
                        app.loader, app.meta, ...
                        max(app.pipelineMeta.reconSettings.N_x, 300), ...
                        app.pipelineMeta.reconSettings.FieldOfView_X);
                    updateStatus(app, sprintf('Speed of sound: %.2f m/s', ...
                        app.pipelineMeta.reconSettings.SpeedOfSound));
                else
                    app.pipelineMeta.reconSettings.SpeedOfSound = app.SpeedOfSoundEditField.Value;
                end

                % Parallel processing
                if app.ParallelCheckBox.Value
                    app.pipelineMeta.parallelSettings.parallelType = 'multiple';
                else
                    app.pipelineMeta.parallelSettings.parallelType = 'single';
                end

                % Get selected endmembers
                endmembers = {};
                if app.HbCheckBox.Value
                    endmembers{end+1} = 'Deoxyhemoglobin';
                end
                if app.HbO2CheckBox.Value
                    endmembers{end+1} = 'Oxyhemoglobin';
                end
                if app.WaterCheckBox.Value
                    endmembers{end+1} = 'Water';
                end
                if app.FatCheckBox.Value
                    endmembers{end+1} = 'Fat';
                end
                if app.IR800CWCheckBox.Value
                    endmembers{end+1} = 'IR800CW';
                end
                if app.ICGCheckBox.Value
                    endmembers{end+1} = 'ICG';
                end

                % Add custom endmembers
                customStr = strtrim(app.CustomEndmemberEditField.Value);
                if ~isempty(customStr)
                    customList = strsplit(customStr, ',');
                    for i = 1:numel(customList)
                        endmembers{end+1} = strtrim(customList{i});
                    end
                end

                if isempty(endmembers)
                    uialert(app.UIFigure, 'Please select at least one endmember.', 'Error');
                    app.RunButton.Enable = 'on';
                    app.StopButton.Enable = 'off';
                    return;
                end

                app.pipelineMeta.unmixSettings.EndmemberNames = endmembers;
                app.pipelineMeta.unmixSettings.UnmixSolverType = app.UnmixSolverDropDown.Value;
                app.pipelineMeta.unmixSettings.MultispectralStateFilterType = app.StateFilterDropDown.Value;

                % Create processing objects
                updateStatus(app, 'Creating processing objects...');
                app.preFilter = recon.MSOTPreFilter(app.pipelineMeta.filterSettings);
                app.reconstructor = recon.ReconSystem(app.pipelineMeta.reconSettings);

                app.pipelineMeta.unmixSettings = util.importSpectra(app.pipelineMeta.unmixSettings);
                unmixingArgs = app.pipelineMeta.unmixSettings;
                app.unmixer = unmix.UnmixSystem(unmixingArgs);
                Ncomps = numel(unmixingArgs.EndmemberNames);

                updateStatus(app, sprintf('Unmixing components: %s', strjoin(endmembers, ', ')));

                % Define reconstruction set
                app.reconstructionSet = 1:size(app.loader.LaserEnergy, 2);
                NFrames = numel(app.reconstructionSet);
                NRep = app.meta.RepNum;

                % Pre-allocate arrays
                app.reconImageStructure = zeros(app.pipelineMeta.reconSettings.N_y, ...
                    app.pipelineMeta.reconSettings.N_x, NFrames);
                app.unmixImageStructure = zeros(app.pipelineMeta.reconSettings.N_y, ...
                    app.pipelineMeta.reconSettings.N_x * Ncomps, NRep);
                app.wavelength = zeros(NFrames, 1);

                % Update frame slider
                app.FrameSlider.Limits = [1 NFrames];
                app.FrameSlider.Value = 1;

                updateStatus(app, sprintf('Processing %d frames...', NFrames));

                % Main processing loop
                for k = app.reconstructionSet
                    if app.stopProcessing
                        updateStatus(app, 'Processing stopped by user.');
                        break;
                    end

                    % Update progress
                    app.ProgressGauge.Value = (k / NFrames) * 100;
                    app.CurrentFrameLabel.Text = sprintf('Frame: %d/%d', k, NFrames);

                    % Load frame data
                    dataFrame = app.loader(k, {});
                    app.wavelength(k) = dataFrame.Meta.Wavelength;

                    % Apply prefilter
                    filteredSignal = app.preFilter(dataFrame);

                    % Reconstruct
                    [reconVec, ~] = app.reconstructor(filteredSignal);
                    app.reconImageStructure(:,:,k) = reconVec;

                    % Unmix
                    inputImage.Data = reconVec;
                    inputImage.Meta = struct(dataFrame.Meta);
                    [unmixFrame, ~] = app.unmixer(inputImage);
                    app.unmixImageStructure(:,:,k) = reshape(unmixFrame.Data, ...
                        [app.pipelineMeta.reconSettings.N_y, ...
                        app.pipelineMeta.reconSettings.N_x * Ncomps]);

                    % Display images
                    imagesc(app.ReconAxes, reconVec);
                    axis(app.ReconAxes, 'image');
                    colormap(app.ReconAxes, 'jet');
                    colorbar(app.ReconAxes);
                    title(app.ReconAxes, sprintf('Frame %d (%.0f nm)', k, app.wavelength(k)));

                    unmixDisplay = reshape(unmixFrame.Data, ...
                        [app.pipelineMeta.reconSettings.N_y, ...
                        app.pipelineMeta.reconSettings.N_x * Ncomps]);
                    imagesc(app.UnmixAxes, unmixDisplay);
                    axis(app.UnmixAxes, 'image');
                    colormap(app.UnmixAxes, 'jet');
                    colorbar(app.UnmixAxes);
                    title(app.UnmixAxes, 'Unmixed Components');

                    drawnow;

                    if mod(k, 10) == 0
                        updateStatus(app, sprintf('Processed frame %d/%d', k, NFrames));
                    end
                end

                if ~app.stopProcessing
                    app.ProgressGauge.Value = 100;
                    updateStatus(app, 'Processing complete!');
                    app.FrameSlider.Enable = 'on';
                    app.SaveResultsButton.Enable = 'on';
                    uialert(app.UIFigure, 'Processing completed successfully!', ...
                        'Success', 'Icon', 'success');
                end

            catch ME
                uialert(app.UIFigure, sprintf('Error during processing: %s', ME.message), 'Error');
                updateStatus(app, sprintf('ERROR: %s', ME.message));
            end

            % Re-enable buttons
            app.RunButton.Enable = 'on';
            app.StopButton.Enable = 'off';
        end

        function StopButtonPushed(app, ~, ~)
            app.stopProcessing = true;
            updateStatus(app, 'Stopping processing...');
        end

        function FrameSliderValueChanged(app, ~, ~)
            if isempty(app.reconImageStructure)
                return;
            end

            frameNum = round(app.FrameSlider.Value);

            % Display reconstructed image
            reconVec = app.reconImageStructure(:,:,frameNum);
            imagesc(app.ReconAxes, reconVec);
            axis(app.ReconAxes, 'image');
            colormap(app.ReconAxes, 'jet');
            colorbar(app.ReconAxes);
            title(app.ReconAxes, sprintf('Frame %d (%.0f nm)', frameNum, app.wavelength(frameNum)));

            % Display unmixed image
            unmixVec = app.unmixImageStructure(:,:,frameNum);
            imagesc(app.UnmixAxes, unmixVec);
            axis(app.UnmixAxes, 'image');
            colormap(app.UnmixAxes, 'jet');
            colorbar(app.UnmixAxes);
            title(app.UnmixAxes, sprintf('Unmixed Components - Frame %d', frameNum));
        end

        function SaveResultsButtonPushed(app, ~, ~)
            if isempty(app.OutputFolderEditField.Value)
                uialert(app.UIFigure, 'Please specify an output folder.', 'Error');
                return;
            end

            try
                outputFolder = app.OutputFolderEditField.Value;
                if ~exist(outputFolder, 'dir')
                    mkdir(outputFolder);
                end

                % Reshape results
                NRep = app.meta.RepNum;
                NFrames = size(app.reconImageStructure, 3);

                recon_img = reshape(app.reconImageStructure, ...
                    [size(app.reconImageStructure,1), size(app.reconImageStructure,2), ...
                    NFrames/NRep, NRep]);
                unmix_img = reshape(app.unmixImageStructure, ...
                    [size(app.unmixImageStructure,1), size(app.unmixImageStructure,2), ...
                    NFrames/NRep, NRep]);

                w = app.wavelength(1:NFrames/NRep);

                % Save to MAT files
                save(fullfile(outputFolder, 'msot_processed.mat'), 'recon_img', 'w');
                save(fullfile(outputFolder, 'unmix_processed.mat'), 'unmix_img', 'w');

                % Save settings
                pipelineMeta = app.pipelineMeta;
                save(fullfile(outputFolder, 'pipeline_settings.mat'), 'pipelineMeta');

                updateStatus(app, sprintf('Results saved to: %s', outputFolder));
                uialert(app.UIFigure, sprintf('Results saved successfully to:\n%s', outputFolder), ...
                    'Success', 'Icon', 'success');

            catch ME
                uialert(app.UIFigure, sprintf('Error saving results: %s', ME.message), 'Error');
            end
        end

        function updateStatus(app, message)
            currentStatus = app.StatusTextArea.Value;
            if numel(currentStatus) > 20
                currentStatus = currentStatus(end-19:end);
            end
            app.StatusTextArea.Value = [currentStatus; {sprintf('[%s] %s', ...
                datestr(now, 'HH:MM:SS'), message)}];
            % Scroll to bottom
            scroll(app.StatusTextArea, 'bottom');
            drawnow;
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = MSOTPipelineGUI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
