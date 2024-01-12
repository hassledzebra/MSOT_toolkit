classdef MSOTSignalLoader < matlab.System
    % MSOTLoader.m Flexible interface for MSOT metadata
    % The main purpose of the MSOTSignalLoader is to provide a means by which
    % the pipeline may load MSOT datasets one frame at a time. In a sense, it is
    % a Dataset object, whose invocation yields a calibrated ShotFrame.
 
  %%%%%%%%%%%%%%%%
  %% PROPERTIES %%
  %%%%%%%%%%%%%%%%
  
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% Public Properties%%%
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    % Public, non-tunable properties
    properties(Nontunable)
        
        % Source filenames.
        HeaderFilename
        DataFilename
        
        % Active data stores.
        Meta@util.msotData
        MemoryMap
        
        % Size of the data panels to load. This should eventually be replaced
        % with the Coordinate attached to the SignalLoader. Remember, even the
        % loader is supposed to be a Filter of a special sort. 
        DataPanelSizeY
        DataPanelSizeX
        
        % Data for corrections. Note: Even if we refactor the corrections
        % downward, it might be helpful to turn these into Dependent properties.
        CorrectionFactors
        LaserEnergy
        
        FactoryLaserEnergy
        
        WaterPathLength
        WaterAbsCoefficient
        
        % Impulse response function. Note: This could theoretically attach to
        % the dataset, but the IRF could be different for each frame, so leave
        % this. 
        IRF
        
        % Correction variables. 
        UndoScalingCorrection@logical = true;
        UndoLaserEnergyCorrection@logical = false;
        CorrectWithAverageLaserEnergy@logical = false;
    end

    
    
    
 
  %%%%%%%%%%%%%%%%
  %%  METHODS   %%
  %%%%%%%%%%%%%%%%
    
    methods
        
        % Constructor
        function obj = MSOTSignalLoader(varargin)
            % TODO: Check filename, determine if .msot or .mat or .json
            % TODO: Handle different loading schemes. 
            if nargin>0 && ~isa(varargin{1},'struct')
                setProperties(obj,nargin,varargin{:},'HeaderFilename')
                
                % Gets the basic file-level information for the metadata.
                numel(varargin)
                obj.Meta=util.msotData(obj.HeaderFilename); 
                % disp(obj)
                % disp(obj.Meta)
                
                % Build the metadata file structure. 
                calculateStructure(obj.Meta);

                setupImpl(obj);
                
            elseif nargin==1 && isa(varargin{1},'struct')
                structFieldNames = fieldnames(varargin{1});
                structFieldValues = struct2cell(varargin{1});
                
                [~,B] = ismember('MsotFile',structFieldNames);
                
                if B~=0
                   structFieldNames{B} = 'HeaderFilename';
                    
                end
                
                objProps = properties(obj);
                [~,keepInds] = ismember(objProps,structFieldNames);
                
                concatStruct = [structFieldNames(keepInds(keepInds~=0)),structFieldValues(keepInds(keepInds~=0))].';
                
                setProperties(obj,numel(concatStruct),concatStruct{:})
                % Gets the basic file-level information for the metadata.
                numel(varargin)
                obj.Meta=util.msotData(obj.HeaderFilename); 
                
                % Build the metadata file structure. 
                calculateStructure(obj.Meta);

                setupImpl(obj);
                
            else
                % Nothing
            end
            
                
        end
        
        
        % Loads the impulse response function for the data. 
        function irf = get.IRF(obj)
               filename=strrep(obj.HeaderFilename,'.msot','.irf');
               fID=fopen(filename,'r','l');
               irf=fread(fID,Inf,'double');
               fclose(fID);
            switch obj.Meta.DataModelVersion
                case {'2.2','2.3'} % v3.8
                    irf = fftshift(irf);
                otherwise
            end
               
        end
        
        
    end
    
    
    
    
    
    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Set up the memory map.
            meta = obj.Meta;
            hwdesc = obj.Meta.HWDesc
            
            obj.DataPanelSizeY=obj.Meta.MeasurementDesc.RecordLength;
            obj.DataPanelSizeX = obj.Meta.HWDesc.NumDetectors;
            %obj.DataPanelSizeX = 512;
           
            obj.DataFilename=obj.Meta.FileName;
            obj.MemoryMap=memmapfile(obj.DataFilename,'Format',...
                {'uint16',double([obj.DataPanelSizeY,obj.DataPanelSizeX]),'Data'});
            
            
            % Populate all the correction data. 
            obj.CorrectionFactors=[obj.Meta.ScanFrames(:).CorrectionFactor];
            obj.LaserEnergy=[obj.Meta.ScanFrames(:).LaserEnergy];
            
            measDesc = obj.Meta.MeasurementDesc;
            
            obj.FactoryLaserEnergy=[measDesc.FactoryLaserEnergyTable];
            obj.WaterPathLength=[measDesc.PathLengthInWater];
            obj.WaterAbsCoefficient=[measDesc.WaterAbsorptionCoeff];
            
            % Create the coordinate system to affix to the individual frames.
            % 
            %    sourceTimes =
            %    "eval:(transpose(double(1:DataPanelSizeY).*1./double(meta:HWDesc.SamplingFrequency)))"
            %    sourceXdcr = eval:(transpose(double(1:DataPanelSizeX)));        
        end
        
        function [y,loadLog] = stepImpl(obj,u,varargin)
            loadLog = {['Loaded ' num2str(u)]};
            
            if numel(u)==1 % Load a single data frame. 
                y.Meta=obj.Meta.ScanFrames(u); % Metadata.
                y.Data=obj.MemoryMap.Data(y.Meta.IDOffset+1).Data; % Raw signal.
                % y.Coordinates = obj.frameCoordinateSystem();
                y = obj.applyCorrections(y);
            else % Load the array (reverse allocation).
                for k=numel(u):-1:1
                    y(k).Meta=obj.Meta.ScanFrames(u(k));
                    y(k).Data=obj.MemoryMap.Data(y(k).Meta.IDOffset+1).Data;
                    % y(k).Coordinates = obj.frameCoordinateSystem();
                    y(k) = obj.applyCorrections(y(k));
                end
            end
            
        end
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Operator functions %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Applies whatever scaling and laser energy corrections need to be
        % applied.
        % TODO: Act on arrays. 
        % TODO: Refactor downward so that ScanFrames can calibrate themselves.
        function correctedFrame = applyCorrections(obj,rawFrame)
            frameData = rawFrame.Data;
            frameMeta = rawFrame.Meta;
            
            assert(isa(frameData,'uint16')); % Raw frame data should ALWAYS be uint16.
            
            frameData = double(frameData); 
            
            %% Corrects for energy variations. 
            if obj.UndoScalingCorrection
                frameData = frameData./frameMeta.CorrectionFactor;
            end
            
            if obj.UndoLaserEnergyCorrection
                frameData = frameData.*frameMeta.LaserEnergy;
            end
            
            if obj.CorrectWithAverageLaserEnergy
                frameData = frameData./frameMeta.LaserEnergy;
            end
            
            correctedFrame.Data = frameData;
            correctedFrame.Meta = frameMeta;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Administrative functions %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        function resetImpl(obj)
            % Reset loader to raw state. Right now this serves no function.
        end

        function releaseImpl(obj)
            % Release resources, such as file handles. 
        end

        function num = getNumInputsImpl(obj)
            % Define total number of inputs for system with optional inputs
            num = 2;
            % if obj.UseOptionalInput
            %     num = 2;
            % end
        end

        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s

            loadObjectImpl@matlab.System(obj,s,wasLocked);
            % Set private and protected properties

            % Set public properties and states
            obj.Meta = matlab.System.loadObject(s.props.Meta);
            obj.MemoryMap = matlab.System.loadObject(s.props.MemoryMap);
        end
        
        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj
            
            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);
            s.props.Meta=matlab.System.saveObject(obj.Meta);
            s.props.MemoryMap=matlab.System.saveObject(obj.MemoryMap);
            
        end
        
    end
    
    
%% Helper methods for performing operations. 
    methods (Static)
    end
end



















