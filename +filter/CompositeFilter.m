classdef CompositeFilter < filter.Filter
    % Iterator Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        FilterList
        ChildFilters
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)

    end
    methods
        function obj = CompositeFilter(varargin)
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
            % Perform one-time calculations, such as computing constants
            % Set up the iterating filter. 
        end

        function [itFrame,stateLog] = doStep(obj,frame,varargin)
            % When doStep is called, it should always have a frame as its
            % first entry. If there are variable arguments remaining, there
            % should be N_childfilters arguments, and they should be cell
            % arrays that are dumped as part of the calling. 
            
%             nVargs = numel(varargin);
            % If we have more than one variable argument, the first should
            % still be the frame, but we will want to past the rest of the vargs
            % in to the child filters
%             if nVargs>1 
%                 frame = varargin{1};
%                 nVargs = nVargs-1;
%                 varargin = varargin(2:end);
%             elseif nVargs==1
%                 frame = varargin{1};
%                 nVargs = 0;
%                 varargin = {};
%             else
%                 frame = 0;
%                 varargin = {};
%             end
             
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            locFilters = obj.ChildFilters;
            DEBUG = obj.DoDebug;
            % For parallel version, set up parallel environment.
            itFrame = frame;
            stateLog = {obj.FilterTag};
%             itFrame.tag = varargin;
            
            for k = 1:numel(locFilters) 
%                 if isa(locFilters{k},'filter.Filter')
%                     [itFrame,itLog] = feval(locFilters{k},itFrame);
%                 else
                [itFrame,itLog] = feval(locFilters{k},itFrame,varargin{:});
%                 end
%                [itFrame,itLog] =
%                feval(locFilters{k},itFrame,varargin{k}{:}); % New call. 
               statelogCop = stateLog;
               stateLog = [stateLog {itLog}];
                
            end
            
                
            stateLog = [{''} stateLog];
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        function releaseImpl(obj)
            for k = numel(obj.ChildFilters):-1:1 
                
               release(obj.ChildFilters{k});
            end
        end
    end
end
