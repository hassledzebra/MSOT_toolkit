classdef Iterator < filter.Filter
    % Iterator Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        IterationList
        IteratingFilter;
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)

    end

    
    methods 
    function obj = Iterator(varargin)
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
        function doSetup(obj,u,varargin)
            % Perform one-time calculations, such as computing constants
            % Set up the iterating filter. 
        end

        function [outFrame,stateLog] = doStep(obj,u,varargin)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            locIterList = obj.IterationList;
            locFilter = obj.IteratingFilter;
            DEBUG = obj.DoDebug;
            
            
            for k = locIterList
               if DEBUG
                   u.IterationIndex = k;
               end
               
               [outFrame,stateLog] = feval(locFilter{1},k,k);
                
            end
                
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        function releaseImpl(obj)
           release(obj.IteratingFilter{1}); 
        end
    end
end
