classdef MSOTPreFilter < matlab.System
    % MSOTPREFILTER preprocesses MSOT data for passing to a reconstructor
    %
    
    % Public, tunable properties
    %#function cheby1 cheby2 designfilt butter sosfilt
    properties
        
        % Signal parameters
        N_Samples
            SampleRate
        N_Transducers
        
        SignalRegion % denotes where the expected signal is generated. Used for Wiener estimation. Should be same as for reconstructor
        IgnoreRegion % denotes where transients are expected and so should be ignored for e.g. calculations.
        
        
        % Wiener filter parameters (only used if RemoveIRF=true)
        WienerLambda
        SignalIRF
        SignalMTF
        
        
        % High-pass filter parameters
        HpFilterFrequency
        HpFilterOrder
        HpFilterRipple
        HpFilterAttenuation
        
        
        % Low-pass filter parameters
        LpFilterFrequency
        LpFilterOrder
        LpFilterRipple
        LpFilterAttenuation
        
        
        % Water atten. correction 
        WaterAbsorptionCoeffs
        PathLengthInWater
        
        HighpassFilterArgs = {"designfilt",'highpassfir', ...        % Response type
                                'FilterOrder','obj:HpFilterOrder',...
                                'PassbandFrequency',5000, ...
                                'StopbandFrequency',4500, ...
                                'SampleRate','obj:SampleRate'};
                            
%         HighpassFilterArgs = {'designfilt','highpassfir', ...        % Response type
%                                 'FilterOrder','obj:HpFilterOrder',...
%                                 'PassbandFrequency','obj:HpFilterFrequency', ...
%                                 'StopbandFrequency','eval:(obj.HpFilterFrequency*0.9)', ...
%                                 'SampleRate','obj:SampleRate'};
%                             
% %                 'PassbandFrequency','obj:LpFilterFrequency', ...
% %                 'StopbandFrequency','obj:LpFilterFrequency', ...
        LowpassFilterArgs = {'designfilt', 'lowpassfir', ...        % Response type
                'FilterOrder','obj:LpFilterOrder',...
                'PassbandFrequency',7e6, ...
                'StopbandFrequency',7.7e6, ...
                'SampleRate','obj:SampleRate'};
            
%         HighpassFilterArgs = {"cheby1", 4, 0.01, 2*5000/40E6*1.46, "high"};
%         LowpassFilterArgs = {"cheby1", 8,0.3150,"low"};
%         filter1_LPF = {"cheby1", 8, 0.01, 2*f_LPF/Fs * 0.9, "low"};
%         filter1_HPF = {"cheby1", 4, 0.01, 2*f_HPF/Fs * 1.46, "high"};

        doDeMean = true;
        doRemoveIRF                = true;
        doHighpass = true;
        doLowpass = true;
        doCorrectWaterAttenuation  = true;

    end
    
    % Pre-computed constants
    properties(Access = private)
        HiFilter 
            % Cell array containing all of the parameters necessary to apply
            % filtfilt. 
        LoFilter
            % Cell array containing all of the parameters necessary to
            % apply filtfilt.
        
    end
    


    methods
        function obj = MSOTPreFilter(varargin)
            % Support name-value pair arguments when constructing object
            % if we have any inputs and the first input isn't a struct,
            % assume that we're loading up in Name/Value mode. 
            
            if nargin>0 && ~isa(varargin{1},'struct') 
                setProperties(obj,nargin,varargin{:})
                
            % Otherwise, if there's only one input and it's a struct,
            % assume that it contains the fields we need to create the
            % filter.
            elseif nargin==1 && isa(varargin{1},'struct')
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
        
    end
    
    methods(Access = protected)
        function setupImpl(obj)
            % TODO: Give more options here, similar to iThera's usage. 
            
            
            % Replace all instances of 'obj:Prop' with the corresponding
            % property.
            lpfa = util.replaceVals(obj,obj.LowpassFilterArgs);
            hpfa = util.replaceVals(obj,obj.HighpassFilterArgs);
            
%             lpfa2 = obj.LowpassFilterArgs;
%             
%             
%             C = 
%             
%             
%             for k = 1:numel(lpfa)
%                 if strncmp(lpfa{k},'obj:',4)
%                    lpfa{k} = obj.(lpfa{k}(5:end));
%                 end
%             end
%             
%             
%             for k = 1:numel(lpfa)
%                 if strncmp(lpfa{k},'obj:',4)
%                    lpfa{k} = obj.(lpfa{k}(5:end));
%                 end
%             end
%             
%             for k = 1:numel(lpfa2)
%                 if strncmp(lpfa2{k},'obj:',4)
%                    lpfa2{k} = eval([strrep(lpfa2{k},':','.'),'*1.1']);
%                 end
%             end
%             
%             
%             hpfa_2 = util.replaceValues(obj,obj.HighpassFilterArgs{:});
%             lpfa_2 = util.replaceValues(obj,obj.LowpassFilterArgs{:});
%             
%             hpfa = obj.HighpassFilterArgs;
%             for k = 1:numel(hpfa)
%                 if strncmp(hpfa{k},'obj:',4)
%                    hpfa{k} = obj.(hpfa{k}(5:end));
%                 end
%             end
            
            % Allocate the right number of outputs to capture.
            switch lpfa{1}
                case {'butter','cheby1','cheby2','sosfilt'}
                    obj.LoFilter = cell(2,1);
                case 'designfilt'
                    obj.LoFilter = cell(1);
                otherwise
            end
            
            switch hpfa{1}
                case {'butter','cheby1','cheby2','sosfilt'}
                    obj.HiFilter = cell(2,1);
                case 'designfilt'
                    obj.HiFilter = cell(1);
                otherwise
            end
            
            %% Creating filters
            [obj.LoFilter{:}]= feval(lpfa{:});
            [obj.HiFilter{:}]= feval(hpfa{:});
            
            

            
            % TODO: Make this more robust; do we even want to do the
            % normalization? 
            try
                obj.SignalIRF = obj.SignalIRF./sum(obj.SignalIRF.^2);
                obj.SignalMTF = fft(obj.SignalIRF);
            catch
                
            end
            
        end
        
        function y = stepImpl(obj,u)
            
            [structFlag,dataFrameFlag] = deal(0);
            
            % Handle a couple of different cases.
            %   If the system is passed a numeric datatype, treat that as the
            %   data and output a corresponding datatype. 
            %
            %   If the system is passed either a structure or a DataFrame
            %   object, it should have a Data field which contains the
            %   corresponding numeric. The output should also be either a
            %   structure or a DataFrame, with the metadata preserved.
            
            if isstruct(u)
                structFlag = 1;
                try uTemp = u.Data; datField = 'Data'; 
                catch uTemp = u.signal; datField = 'signal'; warning('Deprecated use of "signal" instead of "Data"'); end
            elseif isa(u,'dataObject.DataFrame')
                dataFrameFlag = 1;
                error('Unimplemented');
            else
                uTemp = u;
            end
            
            % Remove the means of each column to correct for offset.  
            % TODO: Make this another step/filter, something like #
            % eval('detrend',uTemp)
            % or
            % eval(@(x) x - mean(x,1); % Probably better for generality. 
            if obj.doDeMean
                uTemp = uTemp - mean(uTemp);
            end
            
            % Perform IRF removal. 
            % TODO: Replace with deconvolution filter.
            if obj.doRemoveIRF
                sigvar=mean(var(uTemp(obj.SignalRegion,:),[],1)); 
                noisevar=mean(var(uTemp(~obj.SignalRegion&~obj.IgnoreRegion,:),[],1));
                testCorrection = 1./obj.SignalMTF.*((abs(obj.SignalMTF).^2)./(abs(obj.SignalMTF).^2+obj.WienerLambda*(noisevar+0.001)./sigvar));

                uTemp = ifftshift(ifft(testCorrection.*fft(uTemp,[],1),[],1),1);
            end
            
            % Perform filtering
            % TODO: Replace with LowpassFilter and HighpassFilter systems.
            
            if obj.doLowpass
            uTemp = filtfilt(obj.LoFilter{:},uTemp);
            end
            if obj.doHighpass
            uTemp = filtfilt(obj.HiFilter{:},uTemp);
            end
            % Apply water corrections. 
            % TODO: Replace with WaterAttenuationCorrectionFilter (or whatever
            % it's called)
            if obj.doCorrectWaterAttenuation
                try
                 water_ua = double(obj.WaterAbsorptionCoeffs(u.Meta.WLNum));
                 water_z  = obj.PathLengthInWater;
                 atten = exp(-water_ua.*water_z);
                 uTemp = uTemp./atten;
                catch
                   warning("Tried and failed to correct for water attenuation correction"); 
                end
            end
            
            
            % Rewrite the filtered data so that the output format matches the
            % input format. 
            if structFlag
                y = u;
                y.(datField) = uTemp;
            elseif dataFrameFlag
                y = copy(u);
                y.Data = uTemp;
                error('Unimplemented');
            else
                y = uTemp;
            end
                
            
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end

        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s

            % Set private and protected properties
            % obj.myproperty = s.myproperty; 

            % Set public properties and states
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end

        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);

            % Set private and protected properties
            %s.myproperty = obj.myproperty;
        end

    end
end
