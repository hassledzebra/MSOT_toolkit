classdef ReconSystem < matlab.System
    %RECONSYSTEM Reconstruction system for converting signal panels to
    %reconstructed images. 

    
    %% PROPERTIES
    properties
        Name
        
        
        % Properties describing the raw data space.
        SamplingFrequency % Sampling rate of the ADC.
        N_Samples                % Number of time points collected by the transducer
        N_Transducers      % Number of transducers in the raw data.
        
        TransducerCoordinatesCartesian_X       % Transducer coordinates in the X direction
        TransducerCoordinatesCartesian_Y       % Transducer coordinates in the Y direction
        
        TransducerCoordinatesCylindrical_Theta      % Transducer coordinates in polar coordinates (theta)
        TransducerCoordinatesCylindrical_R       % Transducer coordinates in polar coordinates (rho or radius)
        
        
        % Properties describing the image space.
        
        N_x                % Number of pixels in the reconstructed Y direction
        N_y                % Number of pixels in the reconstructed X direction
        
        FieldOfView_X              % Lateral field of view (X-dir)
        FieldOfView_Y              % Lateral field of view (Y-Dir)
        
        
        
        % Parameters describing the operator. Some of these will be used by the
        % model, others by the interpolation process.
        SpeedOfSound      % Needs to be tunable at some point
        
        PixelOversampling = 3; % Number of model samples per pixel. 
        N_DiscreteRays = 512;   % For raycast-type models, how many rays to use per projection. 
        N_Projections = 300;    % Number of projections in the model. Effectively 'interpolated transducers'.
        
        
        
        % TODO: Add PreFilter processing.
        PreFilterInitArgs = {};
        PreFilterRunningArgs = {};
        
        % TODO: Add Interpolation settings. 
        InterpolationInitArgs = {};
        InterpolationRunningArgs = {};
        
        
        % RECONSTRUCTION PROCESSING
        ModelType = 'dIMMI';  % The model operator to use for the forward modelling step from image --> data.
                              %   In general, this should either be a
                              %   matrix or a function handle which follows
                              %   the usage rules of LSQR.m, i.e. that it
                              %   passes 'transpose' or 'notranspose' to
                              %   evaluate A.' * x or A * x, respectively. 
        ModelInitArgs = {};
        ModelRunningArgs = {};
        
                              
        ModelSolver = 'lsqr'; % The solver used for reconstruction, which should conform to the following format style:
                              % SOLVER(A,b,xo,params) where:
                              %    - A is a model operator (matrix or
                              %    function)
                              %    - b is the data used for reconstruction.
                              %    - xo is an initial guess used for
                              %    initialization of the reconstruction
                              %    process.
                              %    - params is a structure of name-value
                              %    pairs which dictate how the solver runs.
                              %    This may include such things as
                              %    testCase.regularization constants or maximum
                              %    iterations. This corresponds to the options
                              %    structure that many builtin solvers use.
        SolverInitArgs = {};
        SolverRunningArgs = {};
        SolverHandle          % The handle to the solver, used for invoking it as an feval call. 
        
        LastImage             % The last image reconstructed by the ReconSystem. 
       
                
        dump = true;
        
    end
    
    
    properties (Access=private)
        ModelMatrix                % The substance of the model operator. TODO: Rewrite for handles/general operators.
                                   % Because the general model is a (2,2)
                                   % tensor, we can write it as something like
                                   % D(t,r) = M(t,r;y,x)*I(y,x) which would
                                   % implicitly define not only the operation
                                   % but also the dimensions of the inputs, and
                                   % correspondingly what the dimensions of the
                                   % model matrix itself should be. 
        
        ModelProjectionCoordinates % The coordinates of the projections.
        ModelTimes                 % The timepoints the model uses for sampling. 
        
        ModelDump                  % The dump of model parameters that were inputs to the operator initialization.
        
        
        RawCoordinateGrids
        InterpCoordinateGrids
        InterpolationMethod = 'linear';
        
        MODEL_SCALE_FACTOR=1e6;
    end
        
    
    
    methods
        function obj = ReconSystem(varargin)
            %RECONSYSTEM Construct an instance of this class
            disp('Creating ReconSystem');
            if (nargin > 0) && ( ~isa( varargin{1} , 'struct' ) )
                setProperties(obj,nargin,varargin{:})
            % Otherwise, if there's only one input and it's a struct,
            % assume that it contains the fields we need to create the
            % filter.
            % TODO: Put this in a utility. Something with a call like
            % fun(?ReconSystem,varargin).
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
            
            % If the ModelType is backproject2D and the solver doesn't
            % match, fix that. 
            if strcmp(obj.ModelType , 'backproject2D') && ~strcmp(obj.ModelSolver,'backproject2DSolver')
                warning('Overriding selected solver for backprojection')
                obj.ModelSolver = 'backproject2DSolver';
            end
            
            if ~strcmp(obj.ModelType , 'backproject2D') && strcmp(obj.ModelSolver,'backproject2DSolver')
                warning('Overriding selected model for backprojection')
                obj.ModelType = 'backproject2D';
            end
            
            %% Select which solver we're using. 
            % Should do some more checking around here to ensure that the
            % function exists, along with its corresponding solver. 
            % TODO: UnmixSystem uses RunningArgs{1} for the model itself.
            switch obj.ModelSolver
                case 'lsqr'
                    obj.SolverHandle=str2func('recon.lsqrWrapper');
                case 'nonNegAPCG'
                    obj.SolverHandle=str2func('recon.nonNegAPCGWrapper');
                case 'backproject2DSolver'
                    obj.SolverHandle = str2func(['recon.' obj.ModelSolver]);
                    obj.MODEL_SCALE_FACTOR = 1;
                otherwise
                    obj.SolverHandle = str2func(['recon.' obj.ModelSolver]);
            end
            
            obj.LastImage = zeros(obj.N_y,obj.N_x); % This should be in setup
            obj.N_DiscreteRays = ceil(sqrt(obj.N_y*obj.N_x)); % This needs to move, but where?
            disp('ReconSystem created');
        end
        
        % Rather than these, make a properties block that's private set access
        % and public get access. 
        
        function dump = getModelDump(obj)
            dump = obj.ModelDump;
        end
        function dump = getModelTimes(obj)
            dump = obj.ModelTimes;
        end
        function dump = getModelMatrix(obj)
            dump = obj.ModelMatrix;
        end
        
        
        function S = struct(obj)
            S = builtin('struct',obj);
            S = rmfield(S,'ModelMatrix');
            S.ModelDump = obj.ModelDump;
            S.Initialized = obj.isLocked;
            S.SolverHandle = func2str(obj.SolverHandle);
        end
        
        function constructModel(obj)
                        
            % TODO: Find the code used for writing the model, and save that to a
            % string. Possibly using 'which'?
            % Save all of the parameters that went in, if dump is enabled. 
            
            testTime = deriveModelTimes(obj); % How to generalize this step?
            xdcrCoords = interpXdcrCoords(obj); % This should be the sort of thing handled by Coordinate objects,
                                                %   interpolating in the
                                                %   'natural' space and yielding
                                                %   Cartesian coordinates. 
            [pixX,pixY] = deriveImageCoords(obj,'face'); % Again, something that seems harder to generalize.
                                                         
            modelArgs = populateModelArgs(obj); % ModelArgs should be a structure so it can be referenced order-free.
             inputCoords = {pixY,pixX}; % TODO: Inherit from OutputCoordinateSystem once refactored. 
             outputCoords = {testTime,xdcrCoords};
            
            ModelCode = ['recon.',obj.ModelType];
            modelOutArgs = nargout(ModelCode);
            
%             if nargout(ModelCode)<0
                outCap = cell(1,abs(modelOutArgs)-1);
%             end
            disp(sprintf('Creating model %s....',obj.ModelType));
            [obj.ModelMatrix,outCap{:}] = feval(ModelCode,inputCoords,outputCoords,modelArgs);
            obj.ModelMatrix = obj.ModelMatrix/obj.MODEL_SCALE_FACTOR;
            obj.ModelTimes = testTime;
            disp('Model created.');
            
        end
    end
    
    
    
    
    %% METHODS 
    methods(Access = protected)
        
        
        function setupImpl(obj,dataFrame)
            % If this is the first call, and the model isn't already
            % constructed, construct the model.
            % TODO: Add checks for model existing. 
            disp('Entered setup');
            obj.constructModel;
            
            % Since we can have an initial input which tells us the dimensions
            % of the frame we're using, create the interpolation components. 
            
            if isnumeric(dataFrame)
                rawData = dataFrame;
            else
                rawData = dataFrame.Data;
            end
            
            % TODO: Replace this section with an InterpFilter.
                % TODO: Get this information from dataFrame.
                sourceTimes = (double(1:size(rawData,1)).*1./double(obj.SamplingFrequency)).';
                sourceXdcr = double(1:obj.N_Transducers).';

                targetTimes=double(obj.ModelTimes);
                targetNxdcrs=double(linspace(1,obj.N_Transducers,obj.N_Projections));

                % TODO: Change this so that we're capturing a number of coordinates
                % which is consistent with how we're trying to interpolate.
                rawCoords = cell(2,1);
                interpCoords = cell(2,1);

                [rawCoords{:}] = ndgrid(sourceTimes,sourceXdcr);
                [interpCoords{:}] = ndgrid(targetTimes,targetNxdcrs);

                obj.RawCoordinateGrids = rawCoords;
                obj.InterpCoordinateGrids = interpCoords;
            
            % TODO: Replace this with a call for replaceVals into RunningArgs.
            solverArgs = populateSolverArgs(obj);
            disp('Leaving setup');
        end
        
        
        
        %%% Loading a saved ReconSystem.
        function loadObjectImpl(obj,s,wasLocked)
            
            % Set public properties and states, instantiate model if it was
            % already instantiated. 
            loadObjectImpl@matlab.System(obj,s,wasLocked);
            % TODO: If this is in setupImpl, do we need to use this here? 
            obj.constructModel;
        end
        
        %%% Saving a ReconSystem. Effectively compresses the model to its
        %%% parametric form so that it can be transferred over network
        %%% efficiently. 
        function s = saveObjectImpl(obj)
            
            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);
            
        end
        
        
        
        
        function [reconstructedImage, reconstructionLog] = stepImpl(obj,dataFrame)
            
            %% TODO: Split out this interpolation process. 
            %       Realistically the dataFrame itself should carry info about
            %       the data interior to it in some metadata component. Since
            %       that would require a bit of a refactor, we should rely on
            %       the fact that the mapping is consistent between
            %       reconstructions, i.e. this could be mostly done in the
            %       setupImpl. 
            
            if isnumeric(dataFrame)
                rawData = dataFrame;
            else
                rawData = dataFrame.Data;
            end
            
            % TODO: Add a call to do replaceVals on obj.RunningArgs with the
            % dataFrame as an argument. 
            
            % TODO: Run the PreFilter, then the Interpolator, then the
            % reconstruction.
            % preOperation =
            % 'feval(@interpn,{obj:RawCoordinateGrids},frame:Data,{obj:InterpCoordinateGrids},obj:InterpolationMethod)';
            % or
            % 'eval(interpn(obj.RawCoordinateGrids{:},dataFrame.Data,obj.InterpCoordinateGrids{:},obj.InterpolationMethod))';
            dataPanel = interpn(obj.RawCoordinateGrids{:},rawData,obj.InterpCoordinateGrids{:},obj.InterpolationMethod);
            
            % stepOperation =
            % 'feval:2:(obj:SolverHandle,obj:ModelMatrix,dataPanel(:),{obj:RunningArgs})'
            [y,reconstructionLog] = feval(obj.SolverHandle,obj.ModelMatrix,dataPanel(:),obj.SolverRunningArgs{:});
            
            reconstructedImage = reshape(y, obj.N_y, obj.N_x)/obj.MODEL_SCALE_FACTOR; 
            
            % Collect the image if, for example, we're running the
            % reconstruction as a distributed job. 
            try
                obj.LastImage = reconstructedImage;
            catch
                obj.LastImage = gather(reconstructedImage);
            end
            
        end
        
    end
    
    methods
        % Return the time vector for the model's sampling.
        %   It's unclear how this could be easily generalized. We may want to
        %   let each model figure out its own 'natural' variables, so e.g.
        %   SpeedOfSound and TimeRes would be in InitArgs as variables to be
        %   substituted. The model would therefore need to be bound to an
        %   InitArgs list when the ReconSystem is created. 
        function [timeVector] = deriveModelTimes(obj)
            speedOfSound = obj.SpeedOfSound;
            timeRes      = obj.PixelOversampling;
            
            fieldOfView = obj.FieldOfView_X;
            n_x = obj.N_x;
            r = obj.TransducerCoordinatesCylindrical_R;
            
            rFov = fieldOfView*sqrt(2)/2;
            rClose = min(r-rFov);
            rFar   = max(r+rFov);
            
            timeVector = linspace(rClose/speedOfSound,rFar/speedOfSound,n_x*timeRes);
            timeVector = timeVector(:);
        end
        
        
        % Interpolate the transducer coordinates to match the number of
        % projections we're looking to use.
        %   As before, this is a challenge to generalize, but should probably be
        %   moved down to the Coordinate object for the ReconSystem, once we add
        %   the InterpSystem to the front-end of it. 
        function [xdcrCoords] = interpXdcrCoords(obj)
            
            xdcrXStart = obj.TransducerCoordinatesCartesian_X;
            xdcrYStart = obj.TransducerCoordinatesCartesian_Y;
            N_proj = obj.N_Projections;

            [theta,rho] = cart2pol(xdcrXStart,xdcrYStart);
            [thetaNew] = interp1(1:numel(theta),unwrap(theta),linspace(1,numel(theta),N_proj));
            [rhoNew] = interp1(1:numel(rho),rho,linspace(1,numel(theta),N_proj));
            [sensorx,sensory] = pol2cart(thetaNew,rhoNew);
            
            xdcrCoords = [sensorx(:),sensory(:)];
            
            
        end
        
        function [pixX,pixY] = deriveImageCoords(obj,cornerOrFace)
            
            fovX = obj.FieldOfView_X;
            fovY = obj.FieldOfView_Y;
            Nx = obj.N_x;
            Ny = obj.N_y;
            
            switch cornerOrFace
                case 'corner'
                    pixX = linspace(-fovX/2,fovX/2,Nx);
                    pixY = linspace(-fovY/2,fovY/2,Ny);
                case 'face'
                    dx = fovX/Nx;
                    dy = fovY/Ny;
                    pixX = linspace((-fovX/2+dx/2),(fovX/2-dx/2),Nx);
                    pixY = linspace((-fovY/2+dy/2),(fovY/2-dy/2),Ny);
                otherwise
            end
            pixX = pixX(:);
            pixY = pixY(:);
        end
        
        function [args] = populateModelArgs(obj)
            args.SpeedOfSound = obj.SpeedOfSound;
            switch obj.ModelType
                case 'dIMMI'
            args.N_DiscreteRays = obj.N_DiscreteRays;
                case 'CDMMI'
                case 'backproject2D'
                otherwise
                    error('Unrecognized model');
            end
            
            
            
%             for k = 1:numel(lpfa)
%                 if strncmp(lpfa{k},'obj:',4)
%                    lpfa{k} = obj.(lpfa{k}(5:end));
%                 end
%             end
            obj.ModelInitArgs{1} = args;
        end
        
        function [args] = populateSolverArgs(obj)
            args = struct();
            switch obj.ModelSolver
                case 'lsqr'
%                     args.LastImage = 
                case 'nonNegAPCG'
                     %args.LastImage = zeros((obj.N_x.*obj.N_y),1);
                case 'backproject2DSolver'
                    args.Nproj = obj.N_Projections;
                    args.timeArray = obj.InterpCoordinateGrids{1}(:,1);
                    args.Nt = numel(args.timeArray);
                otherwise
                    error('Unrecognized solver');
            end
            
            
            
%             for k = 1:numel(lpfa)
%                 if strncmp(lpfa{k},'obj:',4)
%                    lpfa{k} = obj.(lpfa{k}(5:end));
%                 end
%             end
            obj.SolverRunningArgs{1} = args;
        end
        
    end
end

