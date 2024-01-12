classdef Filter < matlab.System
    % StateFilter Transforms DataFrames to other DataFrames, modulated by the
    % filter's internal state. 'Filtering' in this sense implies an estimation
    % of some output from the input. Transforms, for example, are bijective and
    % fully determined the output given the input. Solvers are also filters, but
    % proceed via an iterative process.
    
    % Public, tunable properties
    properties
        
        RunningArgs_Recipe = {};
        RunningArgs = {}
        
        InputCoordinateSystem = Inf;
        OutputCoordinateSystem = Inf;
        
    end
    
    % Public, non-tunable properties
    properties(Nontunable)
        ParentHandle = {} % Points to the encapsulating object, which should be some form of Filter, or empty.
        ChildHandles = {} % Array of Filters that this Filter is encapsulating and/or using, or empty.
        PrevHandle = {}  % Array of Filters that this Filter receive input from, or empty (if a source filter like a loader)
        NextHandle = {} % Array of Filters that this Filter outputs to, or empty (if a sink filter like a terminal writer)
        
        % Composite filters have a relationship where the first subfilter has
        % both its ParentHandle and its PrevHandle as the same object, while the
        % last subfilter has both its ParentHandle and its NextHandle as the
        % same object.
        InitArgs = {}
        
    end
    
    % Pre-computed constants
    properties(SetAccess = protected,GetAccess = public)
        
        FilterID = '';
        FilterAlias = '';
        FilterTag = '';
        
        
        
        LastInputFrame = {} % Last input DataFrame to the StateFilter.
        LastInputArgs = {} % Last set of input arguments.
        
        LastOutputFrame = {} % Last output of the StateFilter.
        LastOutputArgs = {} % Last set of output arguments
        
        
        DoDebug = false; 
    end
    properties(Dependent)
        
        
    end
    
    
    methods
        function obj = Filter(varargin)
            %UNMIXSYSTEM Construct an instance of this class
            %
            
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
        
        function outputCoords = get.OutputCoordinateSystem(obj)
            
            % If we have an output system already, return that.
            if ~isempty(obj.OutputCoordinateSystem)
                outputCoords = obj.OutputCoordinateSystem;
            else % otherwise we should start asking around.
                outputCoords = obj.getOutputCoordinateSystemImpl;
            end
        end
        
        function inputCoords = get.InputCoordinateSystem(obj)
            if ~isempty(obj.InputCoordinateSystem)
                inputCoords = obj.InputCoordinateSystem;
            else % otherwise we should start asking around.
                inputCoords = obj.getInputCoordinateSystemImpl;
            end
        end
        
    end
    
    methods(Abstract)
        % getOutputCoordinateSystem % Coordinate systems can only be
        % defined, not necessarily redefined. The input coordinate system
        % for a Solver is the same as the output coordinate system of its
        % Model, which itself is defined by ModelArgs/InitArgs and which
        % gets ITS input args from the output args of the Solver.
        % getInputCoordinateSystem
        
        getOutputCoordinateSystemImpl(obj)
        getInputCoordinateSystemImpl(obj)
    end
    
    methods(Abstract, Access = protected)
        
        doSetup(obj,u,varargin)
        % Initialize all sub-filters by calling THEIR setups and passing the
        % InitArgs through. This could probably be moved to the Protected block
        % and the details of initialization left up to each implementation.
        
        doStep(obj,u,varargin)
        
        
        
    end
    
    methods(Access = protected)
        
        
        % Setting up the Filter object. The object will gain an alias based
        % on its class, and a tag to refer to it which is derived from the
        % UUID generated at instantiation time. 
        function setupImpl(obj,varargin)
            
             obj.writeLog('Entering setup');
            
            obj.FilterAlias = genvarname(regexprep(class(obj),'\.','_DOT_'));
            obj.FilterID = genvarname(['id' regexprep(char(java.util.UUID.randomUUID),'-','_')]);
                shortTag = regexp(obj.FilterID,'id([\w]+?)_','tokens');
                shortName = strsplit(obj.FilterAlias,'_DOT_');
            obj.FilterTag = strcat(shortName{end}, '_',shortTag{1});
            if iscell(obj.FilterTag)
                obj.FilterTag = obj.FilterTag{1};
            end
            
             obj.writeLog(sprintf('Filter Alias: %s',obj.FilterAlias));
             obj.writeLog(sprintf('Filter ID: %s',obj.FilterID));
             obj.writeLog(sprintf('Filter Tag: %s',obj.FilterTag));
            
            % Replace the initialization arguments 
            
            doSetup(obj,varargin{:});
            
        end
        
        
        
        
        % Perform the actual behavior of the filter. 
        function [y,modArgs] = stepImpl(obj,frame,varargin)
            if obj.DoDebug
             obj.writeLogNoNewLine('Running step...');
            end
            % Replace any frame: calls in the RunningArgs. 
            for k = 1:numel(obj.RunningArgs_Recipe)
                repSig = obj.RunningArgs_Recipe{k};
                if isstring(repSig) || (ischar(repSig))
                    if contains(repSig,'frame:')
                        evalSig = regexprep(repSig,"frame:\((.*)\)|frame:(.*)",'frame.$1');
                        obj.RunningArgs{k} = eval(evalSig);
                    else
                        obj.RunningArgs{k} = obj.RunningArgs_Recipe{k};
                    end
                else
                    obj.RunningArgs{k} = obj.RunningArgs_Recipe{k};
                end
            end
            
            
            
            tic;
            try
                [y,modArgs] = doStep(obj,frame,varargin);
            catch ME
               rethrow(ME); 
            end
            filterRecord.opTime = toc;
            if obj.DoDebug
            fprintf('%s\n',sprintf(' Ran doStep in %d seconds',filterRecord.opTime));
            end
%             obj.writeLog(sprintf('\b Ran doStep in %d seconds',filterRecord.opTime));
            
            
            
            
            
            
            % TODO: Extend this and streamline a little (Class folder?) to
            % allow for more thorough diagnostics. 
            if isempty(modArgs)
                modArgs = filterRecord;
            else
                if iscell(modArgs)
                    modArgs = {modArgs {'FILTER_DIAGNOSTICS'},{filterRecord}};
                elseif isstruct(modArgs)
                    if ismember('FILTER_DIAGNOSTICS',fieldnames(modArgs))
                        error('Illegal use of FILTER_DIAGNOSTICS field in filter implementation!');
                    else
                        modArgs.FILTER_DIAGNOSTICS = filterRecord;
                    end
                end
            end
            obj.LastInputFrame = frame;
            obj.LastInputArgs = varargin;
            
            obj.LastOutputFrame = y;
            obj.LastOutputArgs = modArgs;
            
            
            
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
        %% Other functions not yet considered.
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
        
        function validateInputsImpl(obj,u,varargin)
            % Validate inputs to the step method at initialization
        end
        
        function validatePropertiesImpl(obj)
            % Validate related or interdependent property values
        end
        
        function processTunedPropertiesImpl(obj)
            % Perform calculations if tunable properties change while
            % system is running. Since this should pretty much only consist of
            % changes to the RunningArgs, this should be handling how the
            % RunningArgs should be modified based on changed input.
            % For example, if we have many RunningArgs set already, but we want
            % to change just a few, this should handle the overwriting.
            
            
            
        end
        
        function flag = isInputSizeLockedImpl(obj,index)
            % Return true if the input cannot change size during execution.
            % The second argument is the RunningArgs and should be allowed to change
            % (either as cell arrays or as a structure)
            
            flag = false;
        end
        
        function numIn = getNumInputsImpl(obj)
            numIn = 2;
        end
        
        
        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object
            % configuration, for the command line and System block dialog
            flag = false;
        end
        
       function writeLog(obj,printstr)
           fprintf('[%s]:(%s): %s\n',datestr(now),obj.FilterTag,printstr);
        end 
       function writeLogNoNewLine(obj,printstr)
           fprintf('[%s]:(%s): %s',datestr(now),obj.FilterTag,printstr);
       end 
        
        
    end
    
end
