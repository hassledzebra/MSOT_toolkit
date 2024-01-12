classdef MSOTReconLoader < matlab.System
    % MSOTLoader.m Flexible interface for MSOT metadata
    
   
    % Public, tunable properties
    properties
        LoadType = 'Recons'
    end

    % Public, non-tunable properties
    properties(Nontunable)
        HeaderFilename
        DataFilename
        Meta@util.msotData
        MemoryMap
        ReconNodeID
        TimeOfDay
        
        IRF
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)
        DataPanelSizeY
        DataPanelSizeX
        CorrectionFactors
    end
    
    %
    properties (Hidden,Transient)
        LoadTypeSet = ...
            matlab.system.StringSet({'MetaOnly','Signals','Recons','MSPs'})
    end

    methods
        % Constructor
        function obj = MSOTReconLoader(varargin)
            % Support name-value pair arguments when constructing object
            if nargin<2
                obj.ReconNodeID=1;
            end
            setProperties(obj,nargin,varargin{:},'HeaderFilename','ReconNodeID')
            obj.Meta=util.msotData(obj.HeaderFilename); 

            % Build the metadata file structure. 
            calculateStructure(obj.Meta);

            setupImpl(obj);
            
        end
        
        
        
        % Loads the impulse response function for the data. 
        function irf = get.IRF(obj)
               filename=strrep(obj.HeaderFilename,'.msot','.irf');
               fID=fopen(filename,'r','l');
               irf=fread(fID,Inf,'double');
               fclose(fID);
        end
    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            % Set up the memory map.
            obj.Meta=util.msotData(obj.HeaderFilename);
            RNs=obj.Meta.ReconNode;
            obj.DataPanelSizeY=RNs(obj.ReconNodeID).Resolution;
            obj.DataPanelSizeX=obj.Meta.ReconNode(obj.ReconNodeID).Resolution;
            
            % Navigate to the recon data file. 
            reconDir=strrep(obj.HeaderFilename,[obj.Meta.FriendlyName '.msot'],'RECONs');
            reconBin=fullfile(reconDir,[obj.Meta.ReconNode(obj.ReconNodeID).GUID '.bin']);
            
            
            obj.DataFilename=reconBin;
            
            %% If it's signals...
            obj.MemoryMap=memmapfile(obj.DataFilename,'Format',...
                {'double',double([obj.DataPanelSizeY,obj.DataPanelSizeX]),'recons'},'Repeat',numel(obj.Meta.ReconNode.IDLookup));

%             obj.CorrectionFactors=[obj.Meta.ScanFrames(:).CorrectionFactor];

            
        end

        function y = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            
            y.img=obj.MemoryMap.Data(obj.Meta.ReconNode(obj.ReconNodeID).IDLookup(u)).recons;
            y.meta=obj.Meta.ReconNode(obj.ReconNodeID).Frames(u);
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            
        end
%         function restartMap(obj)
%             
%         end
        function releaseImpl(obj)
            % Release resources, such as file handles
%             obj.MemoryMap.recons=[];
%             obj.MemoryMap=memmapfile(obj.DataFilename,'Format',...
%                 {'double',double([obj.DataPanelSizeY,obj.DataPanelSizeX]),'recons'},'Repeat',numel(obj.Meta.ReconNode.IDLookup));

        end

        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);

            % Set private and protected properties
            %s.myproperty = obj.myproperty;
        end

        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s

            % Set private and protected properties
            % obj.myproperty = s.myproperty; 

            % Set public properties and states
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end

        %% Advanced functions
%         function validateInputsImpl(obj,u)
%             % Validate inputs to the step method at initialization
%         end

%         function validatePropertiesImpl(obj)
%             % Validate related or interdependent property values
%         end

        function ds = getDiscreteStateImpl(obj)
            % Return structure of properties with DiscreteState attribute
            ds = struct([]);
        end

        function processTunedPropertiesImpl(obj)
            % Perform calculations if tunable properties change while
            % system is running
        end

        function flag = isInputSizeLockedImpl(obj,index)
            % Return true if input size is not allowed to change while
            % system is running
            
            % We want to be able to access multiple panels at once, so we
            % should be able to just pass in a vector of IDs that we want
            % to average.
            if index==1
                flag = false;
            end
        end

        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
end
