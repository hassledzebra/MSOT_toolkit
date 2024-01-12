classdef EvalFilter < filter.Filter
    %EvalFilter Performs an 'eval' operation on input data according to an eval
    %string, provided that string follows basic format guidelines.
    %
    % This can obviously be used for built-ins to be used arbitrarily.
    % Can this be used to allow users to run their own filters?
    %
    
    properties
        Property1
    end
    
    methods
        function obj = EvalFilter(varargin)
           obj@filter.Filter(varargin{:}); 
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
            
        end
        
        function [outFrame,stateLog] = doStep(obj,frame,varargin)
            outFrame = frame;
            outFrame.Data = feval(obj.RunningArgs{:});
            
            stateLog = "Evaluation successfully run";
            
        end
    end
end

