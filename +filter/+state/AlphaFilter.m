classdef AlphaFilter < filter.state.StateFilter
    % AlphaFilter Performs alpha filtering on the input data. 
    %


    % Public, tunable properties
    properties
        AlphaCoefficient = 0.2;
        OutputDataSize = 1;
    end
    
    properties(SetAccess = private,GetAccess = public)
    end

    % Pre-computed constants
    properties(Access = private)
    end
methods
        function obj = AlphaFilter(varargin)
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
            end
            
                obj.InternalState = zeros(obj.OutputDataSize);
            
            try
                inputChannelIndex = u.Meta.WLNum;
            catch
                inputChannelIndex = 1;
            end
            
            obj.InternalState(:,:,inputChannelIndex) = rawData;
        end

        function [y,statelog] = doStep(obj,u,varargin)
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
            
            obj.LastInputFrame = u;
            
            try
                inputChannelIndex = u.Meta.WLNum;
            catch
                inputChannelIndex = 1;
            end
            
            obj.InternalState(:,:,inputChannelIndex) = (1-obj.AlphaCoefficient).*obj.InternalState(:,:,inputChannelIndex) + ...
                obj.AlphaCoefficient .* rawData;
            
            if isnumeric(u)
                y = obj.InternalState;
            else
                y.Data = obj.InternalState;
            end
            obj.LastOutputFrame = y;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end