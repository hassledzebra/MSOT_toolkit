classdef SolverSystem < filter.Filter
    %InterpFilter Interpolates the given dataFrame's Data object into a new coordinate
    %system. 
    %
    %
    
    properties
        ModelType
        ModelHandle = []
        ModelOperator
        
        SolverType
        SolverHandle
        MODEL_SCALING_FACTOR = 1;
    end
    
    methods
        function obj = SolverSystem(varargin)
            if (nargin > 0) && ( ~isa( varargin{1} , 'struct' ) )
                setProperties(obj,nargin,varargin{:})
                % Otherwise, if there's only one input and it's a struct,
                % assume that it contains the fields we need to create the
                % filter.
            elseif (nargin==1) && (isa(varargin{1},'struct'))
                structFieldNames = fieldnames(varargin{1});
                structFieldValues = struct2cell(varargin{1});
                objProps = properties(obj);
                [~,keepInds] = ismember(objProps,structFieldNames);
                
                concatStruct = [structFieldNames(keepInds(keepInds~=0)),structFieldValues(keepInds(keepInds~=0))].';
                setProperties(obj,numel(concatStruct),concatStruct{:})
            else
                % Nothing should happen in other cases.
            end
        end
        
        %% Coordinate system management.
        function outputCoords = getOutputCoordinateSystemImpl(obj)
            
            outputCoords = 0;
        end
        
        function inputCoords = getInputCoordinateSystemImpl(obj)
            
            inputCoords  = 0;
            
        end
    end
    methods(Access = protected)
        
        function doSetup(obj,frame,varargin)
            if isempty(obj.ModelHandle)
               obj.ModelHandle = str2func(obj.ModelType);
               % TODO: Add a capture cell group based on nargout which can
               % give us the stateLog.
            end
            if isempty(obj.SolverHandle)
               obj.SolverHandle = str2func(obj.SolverType); 
            end
            
            
            
            obj.ModelOperator = feval(obj.ModelHandle,obj.OutputCoordinateSystem,obj.InputCoordinateSystem,obj.InitArgs{:})./obj.MODEL_SCALING_FACTOR;
        end
        
        function [outFrame,stateLog] = doStep(obj,frame,varargin)
            outFrame = frame;
            
            % Output of solver will in general have some information in it
            % so we will need to generalize this to handle varargout and/or
            % variable capture (e.g. for LSQR if we want the lsvec only
            % sometimes)
            [outFrame.Data,solverLog] = feval(obj.SolverHandle,obj.ModelOperator,frame.Data(:),obj.RunningArgs{:});
            outFrame.Data = reshape(outFrame.Data,obj.OutputCoordinateSystem.getTensorSize')./obj.MODEL_SCALING_FACTOR;
            outFrame.CoordinateSystem = obj.OutputCoordinateSystem;
            
            stateLog.solverLog = solverLog;
            
        end
    end
end

