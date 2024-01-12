classdef DeconvolutionFilter < filter.Filter
    % DeconvolutionFilter Removes a given IRF from the input.
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        Kernel
        IgnoreRegion
        SignalRegion
        RegularizationConstant
        Method
        SignalIRF
        SignalMTF
    end

    % Pre-computed constants
    properties(Access = private)
    end

    methods
        function obj = DeconvolutionFilter(varargin)
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
            % Assign the IRF and calculate the MTF.
            % TODO: Perform a switch based on the method name, then setup
            % the function handle to the target.
            try
                obj.SignalIRF = obj.Kernel./sum(obj.Kernel.^2);
                obj.SignalMTF = fft(obj.SignalIRF);
            catch
                
            end
        end

        function [y,stateLog] = doStep(obj,frame,varargin)
            % Deconvolution process. Takes a DataFrame and deconvolves the
            % impulse response function from it. 
            
            % For a pure signal P(t,r) , temporal convolution with the IRF I(t)
            % yields an output signal Q(t,r) of the same size and shape.
            % Deconvolution is therefore stating that the above process is
            % taking place, but we want to retrieve P from Q. 
            y=frame;
            frameData = frame.Data;
                sigvar=mean(var(frameData(obj.SignalRegion,:),[],1)); 
                noisevar=mean(var(frameData(~obj.SignalRegion&~obj.IgnoreRegion,:),[],1));
                testCorrection = 1./obj.SignalMTF.*((abs(obj.SignalMTF).^2)./(abs(obj.SignalMTF).^2+obj.RegularizationConstant*(noisevar+0.001)./sigvar));

                 y.Data = ifftshift(ifft(testCorrection.*fft(frameData,[],1),[],1),1);
            
           
            
            
            stateLog = 0;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end