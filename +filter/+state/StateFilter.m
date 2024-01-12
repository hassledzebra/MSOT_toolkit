classdef StateFilter < filter.Filter
    % StateFilter Transforms DataFrames to other DataFrames, modulated by the
    % filter's internal state. 'Filtering' in this sense implies an estimation
    % of some output from the input. Transforms, for example, are bijective and
    % fully determined the output given the input. Solvers are also filters, but
    % proceed via an iterative process. 

    % Public, tunable properties
    properties

    end

    % Public, non-tunable properties
    properties(Nontunable)
        % Filters whose execution is chained together to yield an overall filter
        % output. 
        SubFilters 
    end

    % Pre-computed constants
    properties(SetAccess = protected,GetAccess = public)
        
        InternalState % State of this current filter. Somehow relates to the
                      % InternalState s of the SubFilters?
                      
    end
    
    properties(Dependent)
        SubFilterInternalStates % Should dynamically return the internal states,
                                % if available, of the SubFilters. For things
                                % like called fevals, this would correspond to
                                % whatever RunningArgs are currently set up.
        
    end
    methods(Access = public)
        function outputCoords = getOutputCoordinateSystemImpl(obj)
            
             outputCoords = 0;
        end
        
        function inputCoords = getInputCoordinateSystemImpl(obj)
            
            inputCoords  = 0;
            
        end
                
    end  
    methods(Abstract, Access = protected)
        
        doSetup(obj,u,varargin)
        % Initialize all sub-filters by calling THEIR setups and passing the
        % InitArgs through. This could probably be moved to the Protected block
        % and the details of initialization left up to each implementation. 
        
    end
    
    methods(Access = protected)

        % Performs the actual implementation of the filtering algorithm. This is
        % generalized to take a DataFrame object and pass it through each of the
        % filtering steps in the SubFilter chain. 
        % 
        % Input should be a DataFrame, along with any RunningArgs which might
        % have changed.
        
        function [y,modArgs] = doStep(obj,u,varargin)
            % Implement algorithm. 
            y = u;
            
            
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

        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
end
