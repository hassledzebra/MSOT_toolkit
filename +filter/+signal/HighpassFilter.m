classdef HighpassFilter < filter.Filter
    
    
    properties
        FilterArgs
        DesignMethod
        ApplyMethod
    end
    
    methods 
    function obj = HighpassFilter(varargin)
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
            obj.FilterArgs = {feval(obj.DesignMethod{:})};
        end
        function [outFrame,stateLog] = doStep(obj,frame,varargin)
            stateLog = {};
            outFrame = frame;
            
            outFrame.Data = feval(obj.ApplyMethod{:},obj.FilterArgs{:},frame.Data);
        end
        
    end




end