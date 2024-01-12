classdef ModelOperator < matlab.System & matlab.mixin.Copyable
    % ModelOperator Class for describing mappings between spaces.
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.

    % Public, tunable properties
    properties
        Model = 5;
    end

    properties(DiscreteState)

    end

    % Pre-computed constants
    properties(Access = private)

    end
    methods
        function y = plus(a,b)
            y = a.Model + b;
        end
        
        function y = mtimes(a,b)
            if isa(a,'recon.model.ModelOperator')&&isa(b,'recon.model.ModelOperator')
                y = copy(a);
                y.Model = a.Model * b.Model;
            elseif isa(a,'recon.model.ModelOperator')
               y = a.Model * b; 
            elseif isa(b,'recon.model.ModelOperator')
               y = a * b.Model; 
               
            
            end
        end
    end
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
        end

        function y = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            y = obj * u;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
    end
end
