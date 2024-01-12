classdef KalataFilter < filter.Filter
    % KalataFilter Superclass to handle alpha, alpha-beta, and alpha-beta-gamma filtering on the input data.
    %
    
    % Public, tunable properties
    properties
        
    end
    
    properties(SetAccess = public,GetAccess = public)
        
        TrackingIndex
        ProcessNoiseVariance = 2^2;
        MeasurementNoiseVariance = 70^2;
        FilterOrder = 1;
        DeltaT = 1;
        BurnIn = 3;
        
        OutputDataSize % TODO: Refactor to OutputCoordinateSystem using Filter interface.
        % InputDataSize % This should not change for each frame; the
        % measurement operator, however, will. 
        
        ChildFilter % TODO: Refactor as CompositeFilter
        
        LastFrame % TODO: Refactor using Filter interface.
        
        lastTimeSample = 0;
        NTimesSeenChannel = 0;
    end
    
    % Pre-computed constants
    properties(Access = private)
    end
    
    
    methods(Access = public)
        function outputCoords = getOutputCoordinateSystemImpl(obj)
            
             outputCoords = 0;
        end
        
        function inputCoords = getInputCoordinateSystemImpl(obj)
            
            inputCoords  = 0;
            
        end
                
    end    
    methods(Access = protected)
        function doSetup(obj,u)
            disp("DoSetup");
            % Calculate the Kalata coefficients.
            switch obj.FilterOrder
                case 0                    
                    obj.ChildFilter = filter.state.AlphaFilter;
                    obj.ChildFilter.OutputDataSize =  obj.OutputDataSize;
                    feval(obj.ChildFilter,u);
                case 1
                    obj.ChildFilter = filter.state.AlphaBetaFilter;
                    obj.ChildFilter.OutputDataSize = obj.OutputDataSize;
                    feval(obj.ChildFilter,u); 
                case 2
                otherwise
                    error('Filter order for Kalata filter must be 0, 1, or 2.');                    
            end
            
            thisWL = u.Meta.WLNum;
            
            obj.lastTimeSample = -Inf.*ones(obj.OutputDataSize(3),1);
            obj.NTimesSeenChannel = zeros(obj.OutputDataSize(3),1);
            obj.LastFrame = u;
        end
        
        function [y,statelog] = doStep(obj,u,varargin)
            disp("KalataStep");
            channelInd = u.Meta.WLNum;
            channelTimeDiff = u.Meta.RelTime - obj.lastTimeSample(channelInd);
            trackingIndex = filter.state.KalataFilter.calcLambda(channelTimeDiff,obj.ProcessNoiseVariance,obj.MeasurementNoiseVariance);
            
            if obj.NTimesSeenChannel(channelInd)<obj.BurnIn
                switch obj.FilterOrder
                    case 0
                        obj.ChildFilter.AlphaCoefficient = 1;
                    case 1
                        obj.ChildFilter.AlphaCoefficient = 1;
                        obj.ChildFilter.BetaCoefficient = 0.0;
                    case 2
                    otherwise
                end
            else
            
            
                switch obj.FilterOrder
                    case 0
                        [italpha_opt] = filter.state.KalataFilter.calcABG_from_lambda(trackingIndex);
                        [alpha] = filter.state.KalataFilter.calc_scheduling(obj.NTimesSeenChannel(channelInd),italpha_opt);
                        obj.ChildFilter.AlphaCoefficient = alpha;
                    case 1
                        [italpha,itbeta] = filter.state.KalataFilter.calcABG_from_lambda(trackingIndex);
                        [alpha,beta] = filter.state.KalataFilter.calc_scheduling(obj.NTimesSeenChannel(channelInd),italpha,itbeta)
                        obj.ChildFilter.AlphaCoefficient = alpha;
                        obj.ChildFilter.BetaCoefficient = beta;
                    case 2
                    otherwise
                        error('Filter order for Kalata filter must be 0, 1, or 2.');
                end
            end
            obj.lastTimeSample(channelInd) = u.Meta.RelTime;
            obj.NTimesSeenChannel(channelInd) = obj.NTimesSeenChannel(channelInd) + 1;
            obj.TrackingIndex = trackingIndex;
            [y,statelog] = feval(obj.ChildFilter,u);
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end

        function processTunedPropertiesImpl(obj)
            % Perform calculations if tunable properties change while
            % system is running
        end

        function flag = isInputSizeLockedImpl(obj,index)
            % Return true if input size is not allowed to change while
            % system is running
            flag = false;
        end
        
        
    end
    methods(Static)
        
        
        
        function lambda = calcLambda(dt,process_var,noise_var)
            
            sw = sqrt(process_var);
            sn = sqrt(noise_var);
            
            lambda = sw.*(dt.^2)./sn;
            
        end
        
        % calculating a-b-g
        function [alpha, beta, gamma] = calcABG_from_lambda(lambda)
            
            switch nargout
                case 1
                    alpha = (-(lambda.^2) + sqrt(lambda.^4 + 16.*(lambda.^2)))./8;
                case 2
                    r = (4 + lambda - sqrt(8.*lambda + lambda.^2))./(4);
                    alpha = 1 - r.^2;
                    beta = 2.*(2-alpha) - 4.*sqrt(1-alpha);
                    
                case 3
                    b = lambda/2 - 3;
                    c = lambda/2 +3;
                    d = -1;
                    p = c - (b.^2)/3;
                    q = (2.*b.^3)./27 - (b.*c)/3 + d;
                    v = sqrt(q.^2 + (4.*p.^3)/27);
                    z = -(q + v/2).^(1/3);
                    s = z - p./(3.*z)-b./3;
                    alpha = 1 - s.^2;
                    beta = 2.*(1-s).^2;
                    gamma = (beta.^2)./(2.*alpha);
                otherwise
                    disp('Number of output arguments specifies filter order.');
            end
            
        end
        
        % constructing time-advance matrices F from given noise model and dt.
        
        % performing predictive time step.
        
        % Constructing gain matrix
        
        % Applying gain matrix
        
        % Calculating the scheduling based on the sample index.
        function [varargout] = calc_scheduling(k,varargin)
            
            nvarg = numel(varargin);
            switch nvarg
                case 1
                    [alpha] = varargin{:};
                    if isnan(alpha)
                        alpha = -Inf;
                    end
                    
                    itValue_alpha = 1./(k+1);
                    alpha_k = max(itValue_alpha,alpha);
                    varargout = {alpha_k};
                case 2
                    [alpha, beta] = varargin{:};
                    if isnan(alpha)
                        alpha = -Inf;
                    end
                    if isnan(beta)
                        beta = -Inf;
                    end
                    
                    itValue_alpha = (2.*(2.*k + 1)./((k+1).*(k+2)));
                    itValue_beta = (6./((k+1).*(k+2)));
                    
                    alpha_k = max(itValue_alpha,alpha);
                    beta_k = max(itValue_beta,beta);
                    varargout = {alpha_k beta_k};
                case 3
                    [alpha, beta, gamma] = varargin{:};
                    if isnan(alpha)
                        alpha = -Inf;
                    end
                    if isnan(beta)
                        beta = -Inf;
                    end
                    if isnan(gamma)
                        gamma = -Inf;
                    end
                    
                    itValue_alpha = (3.*(3.*k.^2 + 3.*k + 2)./((k+1).*(k+2).*(k+3)));
                    itValue_beta = (18.*(2.*k + 1)./((k+1).*(k+2).*(k+3)));
                    itValue_gamma = (60)./((k+1).*(k+2).*(k+3));
                    
                    alpha_k = max(itValue_alpha,alpha);
                    beta_k = max(itValue_beta,beta);
                    gamma_k = max(itValue_gamma,gamma);
                    varargout = {alpha_k beta_k gamma_k};
                otherwise
            end
        end
        
    end
end