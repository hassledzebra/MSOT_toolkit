classdef AlphaBetaGammaFilter < filter.KalataFilter
    % AlphaFilter Performs alpha filtering on the input data. 
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        AlphaCoefficient
        BetaCoefficient
        GammaCoefficient
    end
    
    properties(SetAccess = private,GetAccess = public)
        State_Velocity;
        State_Acceleration;
    end

    % Pre-computed constants
    properties(Access = private)
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Calculate the Kalata coefficients.
        end

        function y = stepImpl(obj,u)
            % Filtering.
            try y = copy(u);
            catch y = u; 
            end
            
            y.Data = (1-obj.AlphaCoefficient).*InternalState + obj.AlphaCoefficient .* u.Data;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
end