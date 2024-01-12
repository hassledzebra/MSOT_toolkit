function [xSolution,convInfo,iterInfo,resVec,iterationHistory] = nonNegativeAPCG(A,b,varargin)
    % Rewrite testbed of non-negative accelerated projected conjugate gradient.
    
    
    
    %% Input parsing 
        % Check type of model input.
        
    %% Initializing variables and output parameters
    
    
    %% Outer loop
    
    for outer = 1:1
        % Calculate residual
        
        % Store previous gradient
        
        % Calculate new gradient
        
        % Find values which can have a valid improvement.
        
        % Restrict operators and variables to the reduced set.
        
        for inner = 1:1
            % Calculate the beta correction coefficient.
            
            % Update the conjugate direction.
            
            % Project the conjugate direction to the data space.
            
            % Calculate scaling coefficient.
            
            % Update step
            
            % Update residual
            
            % Update gradient
            
            % Record convergence of inner loop. 
            
        end
        
        % Find all indices which are in I1 and where x is greater than zero or
        % where the step direction moves into the positive orthant.
        
        % Along that step direction, and given the whole of x, determine step
        % size.
        
        
    end
        
end


%%
% Helper Functions

function [parsed] = parseInputs(A,b,varargin)
   % Check number of inputs
   
   % Check if options structure or if we have a bunch of name-value pairs.
   
   
   % Check consistency
end

function [parsedOpts] = parseOptions(varargin)
    
    
end

function convergeInfo = calcConvergence(v)
    
end

function beta = calculateBeta(p)
    
end

function step = calculateStepSize(p)
    
end

