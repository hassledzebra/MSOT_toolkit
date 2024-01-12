classdef AlphaBetaFilter < filter.state.StateFilter
    % AlphaFilter Performs alpha filtering on the input data. 
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        AlphaCoefficient = 0.5;
        BetaCoefficient = 0.05;
        OutputDataSize = 1;
    end
    
    properties(SetAccess = private,GetAccess = public)
    end

    % Pre-computed constants
    properties(Access = private)
    end

    
    methods
        function obj = AlphaBetaFilter(varargin)
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
    end
    methods(Access = protected)
        function doSetup(obj,u,varargin)
            disp("AlphaBeta Setup");
            defaultSize = (ndims(obj.OutputDataSize)==1)&&(obj.OutputDataSize == 1);
            
                
            if isnumeric(u) && defaultSize % Numeric, no size instruction
                rawData = u;
                obj.OutputDataSize = size(u);
            elseif isnumeric(u) && ~defaultSize % Numeric, size instruction
                rawData = u;
            elseif ~isnumeric(u) && defaultSize % Structured, no size instr.
                rawData = u.Data;
                obj.OutputDataSize = ~size(u.Data); % Structured, size instr.
            else
                rawData = u.Data;
                if ~isinf(obj.OutputCoordinateSystem)
                    obj.OutputDataSize = obj.OutputCoordinateSystem.getTensorSize;
                else
                end
            end
            
                obj.InternalState.state = zeros(obj.OutputDataSize(:)');
                obj.InternalState.d_state = zeros(obj.OutputDataSize(:)');
            
            try
                inputChannelIndex = u.Meta.WLNum;
            catch
                inputChannelIndex = 1;
            end
            
            obj.InternalState.state(:,:,inputChannelIndex) = rawData;
        end

        function [y,statelog] = doStep(obj,u,varargin)
%             disp("AlphaBeta Step");
            statelog = '';
            % Filtering.
            try y = copy(u);
            catch y = u; 
            end
            
            if isnumeric(u)
                rawData = u;
            else
                rawData = u.Data;
            end
            
            
            try
                inputChannelIndex = u.Meta.WLNum;
            catch
                inputChannelIndex = 1;
            end
            
            try
                dt = u.Meta.RelTime - obj.LastInputFrame.Meta.RelTime;
            catch
                dt = 1;
            end
            
            
            
            % Apply kinematics.
            obj.InternalState.state = obj.InternalState.state + dt.*obj.InternalState.d_state;
            obj.InternalState.d_state = obj.InternalState.d_state;
            
            
            % Calculate residual.
            resid = rawData - obj.InternalState.state(:,:,inputChannelIndex);
            
            
            % Apply correction step.
            obj.InternalState.state(:,:,inputChannelIndex) = obj.InternalState.state(:,:,inputChannelIndex) + obj.AlphaCoefficient.*resid;
            obj.InternalState.d_state(:,:,inputChannelIndex) = obj.InternalState.d_state(:,:,inputChannelIndex) + obj.BetaCoefficient.*resid;

%             obj.InternalState(:,:,inputChannelIndex) = (1-obj.AlphaCoefficient).*obj.InternalState.state(:,:,inputChannelIndex) + ...
%                 obj.AlphaCoefficient .* rawData;
            
            if isnumeric(u)
                y = obj.InternalState.state;
            else
                y.Data = obj.InternalState.state;
            end
            obj.LastInputFrame = u;
            obj.LastOutputFrame = y;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end