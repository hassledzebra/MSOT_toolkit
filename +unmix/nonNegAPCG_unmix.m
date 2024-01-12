function [x,infoLog] = nonNegAPCG_unmix(A,b,varargin)
    %
    nVargs = numel(varargin);
    
    %% Input parsing
%     [sy,sx,sw] = size(b);
%     b = permute(b,[3 1 2]);
%     b = b(:);
    %% Iterating variables
    [m,n] = size(A);
%     sc = n./(sy*sx);
    
    maxIter_outer = 10;
    maxIter_inner = 10;
    
    g = zeros(n,1);
    d = zeros(n,1);
    p = zeros(n,1);
    q = zeros(m,1);
    
    if nVargs == 1
        x = varargin{1}(:);
        if numel(x) ~= size(A,2)
            warning('Initial estimate of X is of wrong dimensions');
            x = A.'*b;
        end
    else
            x = A.'*b;
    end
    x = max(x,0);
    
    %% Convergence variables
    gradientEpsilon = 1E-4;
    xEpsilon = 2.5E-5;
    
    %% Monitoring variables.
    blankPanel = zeros(maxIter_inner,maxIter_outer);
    Linfnorm_G = blankPanel;
    
    %% Setting up log.
    infoLog.convergence = 'Max Iterations Reached';
    
    
    
    %% Outer loop.
    
    for outerIndex = 1:maxIter_outer
        
        r = b - A*x;
%         r = A*x - b;
        g_old = g;
        g = -(A.'*r);
        
        if outerIndex == 1
            g0norm = LinfNorm(g);
        end
        
        % break if the Inf norm for G is reduced.
        if (LinfNorm(g)/g0norm) < xEpsilon
            disp('outer gradient Linf converged');
            break
        end
        
        % Find all indices where x is >0 or where following the gradient
        % downwards will reduce the objective.
        I1 = (x > 0) | ((-g) > 0) ;
        
        d(:) = 0; % reset the step direction.
        
        gouter = g;
        gpreNorm = LinfNorm(g);
        
        %% Inner loop on the reduced set of variables.
        for innerIndex = 1:maxIter_inner
            
            % calculate the correction coefficient.
            beta_cg=max(0, ...
                -( (g(I1).')*(g(I1) - g_old(I1)))./((g_old(I1).')*g_old(I1)));
            
            % Update the conjugate direction.
            p(I1) = -g(I1) + beta_cg*p(I1);
            p(~I1) = 0;
            
            % Project to data space.
            q = A*p; %%%%%%%%%%%%%%%%
            
            alpha_cg = dot(-g,p)./dot(q,q);
            
            d(I1) = d(I1) + alpha_cg * p(I1);
            r = r-alpha_cg*q;
            g_old = g;
            g = -(A.'*r);
            
            Linfnorm_G(innerIndex,outerIndex) = LinfNorm(g);
            
            %% Break if we violate invariants or satisfy convergence.
            if norm(g-g_old) < gradientEpsilon*norm(g)
%                 disp('gradient L2 converged');
                break
            end
            
            if (Linfnorm_G(innerIndex,outerIndex)/gpreNorm) < gradientEpsilon
                disp('gradient Linf converged');
                break
            end
            
            %
            if norm(g)>norm(g_old)
                break
            end
            
            
        end
        
        % Find all indices in I1 where x is >0 or where the step leads us
        % towards positive values.
        I2 = I1 & (x > 0 | d > 0);
        
        searchVector = d.*I2;
        stepSize = 1;
        xprev = x;
        
        alp = calculateStepsize(A,b,x,gouter,searchVector,stepSize);
        
        x = max(x + alp.*searchVector,0);
        
        
        
        %% Convergence
        if norm(xprev-x) < xEpsilon*norm(x)
            infoLog.convergence = 'Solution converged';
            break;
        end
        
        
    end
%     x = permute(reshape(x,[sc,sy,sx]),[2 3 1]);
end




%% Helper functions

function armijoAlpha = calculateStepsize(A,b,x,g,pdir,alpha)
    tau = 0.5;
    c = 0.1;
    
    Ax = A*x;
    Ap = A*pdir;
    
    linObj  = dot(Ax,Ap) + dot(Ap,Ax) - dot(Ap,b) - dot(b,Ap);
    quadObj = dot(Ap,Ap);
    
    mslope = dot(pdir,g);
    
    expectedLoss = -mslope*c;
    
    calcLossApp = -alpha*(linObj+alpha*quadObj);
    
    while(calcLossApp<expectedLoss)
        alpha = alpha*tau;
        calcLossApp = -alpha*(linObj+alpha*quadObj);
        if alpha<1E-3
            break
        end
    end
    armijoAlpha = alpha;
    
end


function [parsed] = parseInputs(A,b,varargin)
    % Check number of inputs
    
    % Check if options structure or if we have a bunch of name-value pairs.
    
    
    % Check consistency of sizes.
end

function [parsedOpts] = parseOptions(varargin)
    
    
end

function convergeInfo = calcConvergence(v)
    
end

function beta = calculateBeta(p)
    
end

function step = calculateStepSize(p)
    
end

function n = LinfNorm(vec)
    n = max(vec(:));
end

function n = L2Norm(vec)
    n = norm(vec(:),2);
end

function n = L1Norm(vec)
    n = norm(vec(:),1);
end

