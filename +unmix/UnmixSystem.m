classdef UnmixSystem < matlab.System
    %UNMIXSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Filesystem location properties
        % TODO: These should be moved out of the UnmixSystem code. Instead, we
        % should add the actual spectra to the InitArgs for the UnmixSystem
        SYSTEM_SPECTRA_FOLDER = util.loadDefault('SYSTEM_SPECTRA_FOLDER');
        USER_SPECTRA_FOLDER = ''; % on Astrocyte this should be the 'work' folder for a Nextflow sesh.
        
        % Image properties
        % Input space
        Ny_Input
        Nx_Input 
        Nc_Input = 1; 
        
        % Output space
        Ny_Output
        Nx_Output
        Nc_Output
        
        
        % Endmembers. Used to create/describe the model, but will actually be
        % derived from the Output coordinate system (see Nc_Output above)
        EndmemberFilenames
        EndmemberNames
        EndmemberWavelengths
        EndmemberWavelength_Unit = 'nm';
        EndmemberSpectrum
        EndmemberSpectrum_Unit = 'cm^{-1}/{mol}'
        
        % Stuff for operation
        % Any running filter which is used to provide a complete
        % multispectral image at each point in time, given the data up
        % to this point. 
        MultispectralStateFilterType = 'slidingFunction';
        MultispectralStateFilter
        StateFilterInitArgs = {{}};
        StateFilterRunningArgs = {{}};
        
        % Actual unmixing procedures.
        UnmixingWavelengths
        UnmixingSpectrum
        UnmixingInitArgs = {{}};
        UnmixingRunningArgs = {{}};
        
        
        MixingModelType = 'linear';
        UnmixSolverType = 'linear';
        
        LastFrame;
    end
    
    properties(Access = private)
        MixingModel
        SolverHandle      
    end
    
    methods
        function obj = UnmixSystem(varargin)
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
        
        function S = struct(obj)
            S = builtin('struct',obj);
            S = rmfield(S,'MixingModel');
            S.Initialized = obj.isLocked;
            S.SolverHandle = func2str(obj.SolverHandle);
        end
    end
    
    methods(Access = protected)
        function setupImpl(obj,firstInput)
            
            % Number of pixels and number of input channels. 
            % TODO: Since this should be a DataFrame, we should be expecting the
            % correct form. 
            if isnumeric(firstInput)
                [nY,nX,niC] = size(firstInput);
            else
                try
                    [nY,nX,niC] = size(firstInput.Data);
                catch
                    [nY,nX,niC] = size(firstInput.data);
                    warning('Deprecated use of "data" instead of "Data"');
                end
            end
            obj.Ny_Input = nY;
            obj.Nx_Input = nX;
            obj.Nc_Input = niC;
            
            obj.Ny_Output = nY;
            obj.Nx_Output = nX;
            
            locWL = obj.UnmixingWavelengths;
            locComponents = obj.EndmemberNames;
            nComp = numel(locComponents);
            
            obj.Nc_Output = nComp;
           
            
            %%% By this point we should know the dimensions and configuration of
            %%% the DataFrame that the multispectral filter should be providing
            %%% us with. We should also know the operator that links together
            %%% each data frame with its multispectral 'parent'.
            
            %% Setup the multispectral state filter.
            filterType = obj.MultispectralStateFilterType;
            
            % setting up an object would be passing in the InitArgs  to the
            % system initialization. The setupImpl would then be called at the
            % end of this setupImpl once we have the rest of the chain
            % configured. Or do we even need to do that, if setupImpl then
            % passes into stepImpl?
            % Moreover, since 
            switch filterType
                case 'sliding'
                case 'slidingFunction'
                    obj.MultispectralStateFilter = str2func('unmix.slidingFunction');
                    feval(obj.MultispectralStateFilter,firstInput,[nY,nX,numel(locWL)]);
                case 'alpha'
                    obj.MultispectralStateFilter = filter.state.AlphaFilter;
                    obj.MultispectralStateFilter.OutputDataSize = [nY,nX,numel(locWL)];
                    feval(obj.MultispectralStateFilter,firstInput,{});
                case 'alphabeta'
                    obj.MultispectralStateFilter = filter.state.AlphaBetaFilter;
                    obj.MultispectralStateFilter.OutputDataSize = [nY, nX, numel(locWL)];
                    feval(obj.MultispectralStateFilter,firstInput,{{}}); 
                case 'kalata'
                    obj.MultispectralStateFilter = filter.state.KalataFilter;
                    obj.MultispectralStateFilter.OutputDataSize = [nY, nX, numel(locWL)];
                    feval(obj.MultispectralStateFilter,firstInput,{{}}); 
                    
                case 'kalman' 
                otherwise
                    error('Unrecognized state filter type');
            end
            
            %% Configuring the actual unmixing process.
            %Setup the unmixing model.
            mixModelType = obj.MixingModelType;
            
            switch mixModelType
                case {'linear'}
                    
                    % Calculate the mixing matrix (N_wl x N_components)
                    locMixMat = zeros(numel(locWL),numel(locComponents));
                    
                    for k = 1 : nComp
                        locMixMat(:,k) = obj.UnmixingSpectrum{k};
                    end
                    
                    obj.MixingModel = locMixMat;
                otherwise
                    error('Unrecognized mixing model');
            end
            
            % Setup the unmixing procedure.
            unmixSolverType = obj.UnmixSolverType;
            
            switch unmixSolverType
                case {'linear','direct'}
                    obj.SolverHandle = str2func('unmix.linearUnmix');
                    obj.MixingModel = obj.MixingModel;
                case {'nnls'}
                    obj.SolverHandle = str2func('unmix.nnlsUnmix');
                    obj.MixingModel = kron(obj.MixingModel,speye(nY*nX));
                case {'nonNegAPCG'}
                    obj.SolverHandle = str2func('unmix.nonNegAPCG_unmix');
                    obj.MixingModel = kron(obj.MixingModel,speye(nY*nX));
                otherwise
                    error('Unrecognized unmixing solver.');
            end
            
            
        end
        
        
        
        function validatePropertiesImpl(obj)
            %% Validate related or interdependent property values
        end
        
        function validateInputsImpl(obj,u)
            % Validate inputs to the step method at initialization
            
        end
        
        
        %%% Loading a saved UnmixSystem.
        function loadObjectImpl(obj,s,wasLocked)
            
            % Set public properties and states, instantiate model.
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end
        
        %%% Saving an UnmixSystem
        function s = saveObjectImpl(obj)
            
            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);
            
        end
        
        function [unmixedImage, unmixedLog] = stepImpl(obj,dataFrame)
            % The input dataframe can either be a single-wavelength image
            % or a multispectral image the size we need for unmixing. In
            % either case, filter it.
            
            multispectralImage = obj.MultispectralStateFilter( dataFrame , obj.StateFilterRunningArgs{:});
            unmixedImage = dataFrame;
            % Multispectral image will be of the form Ny x Nx x Nc.
            [unmixedImage.Data,unmixedLog] = feval(obj.SolverHandle,obj.MixingModel,multispectralImage.Data(:),obj.UnmixingRunningArgs{:});
            
            obj.LastFrame = unmixedImage;
            obj.UnmixingRunningArgs{1} = unmixedImage.Data;
            
        end
    end
end
