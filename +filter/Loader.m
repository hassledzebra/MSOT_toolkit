classdef Loader < filter.Filter
    % Loader takes as an argument an index and uses that to load a DataFrame. 

    % InitArgs should lead to a configuration state where some means of mapping
    % into a dataset occurs. 
    %
    % Support should be prioritized for .bin files (memmapfile), .mat files, and
    % .h5 files. 
    %
    % Abstractly, Loaders are filters that act on Datasets to yield DataFrames.
    % The Dataset is therefore necessarily communicated in InitArgs in order to
    % build the association. 
    % The Loader should be an abstract class still that has an enforced calling
    % structure, but we should allow the data link to be whatever it needs to
    % be. 

    
    methods
        function obj = Loader(varargin)
           obj@filter.Filter;  
           obj.InputCoordinateSystem = [1];
        end
        function outputCoords = getOutputCoordinateSystemImpl(obj)
            outputCoords = 'Under Construction';
        end
        function inputCoords = getInputCoordinateSystemImpl(obj)
            inputCoords = [1]; % Denotes scalar.
        end
        
    end

    
    methods(Access = protected)


                
        function [y,modArgs] = doSetup(obj,u,varargin)
            y = u;
            modArgs = varargin;
        end
        
        function [y,modArgs] = doStep(obj,u,varargin)
             
            y = u;
            modArgs = varargin;
            
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

        function validateInputsImpl(obj,u)
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

            flag = true;
        end

        function flag = isInactivePropertyImpl(obj,prop)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
end
