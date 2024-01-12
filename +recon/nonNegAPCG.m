function [x,infoLog] = nonNegAPCG(A,b,varargin)
    %
    
    
    %% Input parsing
    prs=inputParser;
    
    %Validation functions
    boundOnUnit=@(zzz) ((zzz>=0)&(zzz<=1));
    
    %Require that the inputs are numeric.
    addRequired(prs,'A',@(x) validateattributes(x,{'numeric' 'sparse'},{'2d'}));
    addRequired(prs,'b',@(x) validateattributes(x,{'numeric' 'sparse'},{'2d','column'}));
    
    
    
    
    %Set defaults that are passed to the function.
    
    addParameter(prs,'NMaxIterOuter',11);
    addParameter(prs,'NMaxIterInner',11);
    addParameter(prs,'GradientTol',1E-6);
    addParameter(prs,'ConvergeTol',1E-6);
    addParameter(prs,'c_ls',0.1,boundOnUnit);
    addParameter(prs,'tau_ls',0.5,boundOnUnit);
    addParameter(prs,'nMax_ls',100,@isinteger);
    addParameter(prs,'a_step',10);
    addParameter(prs,'verbose',false);
    addParameter(prs,'regularized',false);
    addParameter(prs,'lambda',0.01);
    addParameter(prs,'x0',zeros(size(A,2),1));
    addParameter(prs,'DEBUG',false);
    
    parse(prs,A,b,varargin{:});
    DEBUG = prs.Results.DEBUG;
    
    %% Iterating variables
    [m,n] = size(A);
    
    maxIter_outer = prs.Results.NMaxIterOuter;
    maxIter_inner = prs.Results.NMaxIterInner;
    
    g = zeros(n,1);
    d = zeros(n,1);
    p = zeros(n,1);
    q = zeros(m,1);
    
    x = A.'*b;
    x = max(x,0);
    
    %% Convergence variables
    gradientEpsilon = prs.Results.GradientTol;
    xEpsilon = prs.Results.ConvergeTol;
    
    %% Monitoring variables.
    blankPanel = zeros(maxIter_inner,maxIter_outer);
    Linfnorm_G = blankPanel;
    
    %% Setting up log.
    infoLog.convergence = 'Max Iterations Reached';
    
    if DEBUG
       z = nan(maxIter_inner,maxIter_outer);
       c = cell(maxIter_inner,maxIter_outer);
       
       normStruc.L0 = z;
       normStruc.L1 = z;
       normStruc.L2 = z;
       normStruc.Linf = z;
       
       g_debug.norms = normStruc;
       g_debug.images = c;
       
       r_debug.norms = normStruc;
       r_debug.images = c;
       
       p_debug.norms = normStruc;
       p_debug.images = c;
       
       d_debug.norms = normStruc;
       d_debug.images = c;
       
       q_debug.norms = normStruc;
       q_debug.images = c;
       
       alphacg_debug = z;
       betacg_debug = z;
       
       I1_debug = cell(1,maxIter_outer);
       I2_debug = cell(1,maxIter_outer);
       x_debug = cell(1,maxIter_outer);
       step_debug = zeros(1,maxIter_outer);
       relNorm = zeros(1,maxIter_outer);
       
        
    end
    
    
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
            
            if DEBUG
                
                [g_debug.norms.L0(innerIndex,outerIndex) , ...
                    g_debug.norms.L1(innerIndex,outerIndex),...
                    g_debug.norms.L2(innerIndex,outerIndex),...
                    g_debug.norms.Linf(innerIndex,outerIndex)] = calcNorms(g);
                
                [r_debug.norms.L0(innerIndex,outerIndex) , ...
                    r_debug.norms.L1(innerIndex,outerIndex),...
                    r_debug.norms.L2(innerIndex,outerIndex),...
                    r_debug.norms.Linf(innerIndex,outerIndex)] = calcNorms(r);
                
                [p_debug.norms.L0(innerIndex,outerIndex) , ...
                    p_debug.norms.L1(innerIndex,outerIndex),...
                    p_debug.norms.L2(innerIndex,outerIndex),...
                    p_debug.norms.Linf(innerIndex,outerIndex)] = calcNorms(p);
                
                [q_debug.norms.L0(innerIndex,outerIndex) , ...
                    q_debug.norms.L1(innerIndex,outerIndex),...
                    q_debug.norms.L2(innerIndex,outerIndex),...
                    q_debug.norms.Linf(innerIndex,outerIndex)] = calcNorms(q);
                
                [d_debug.norms.L0(innerIndex,outerIndex) , ...
                    d_debug.norms.L1(innerIndex,outerIndex),...
                    d_debug.norms.L2(innerIndex,outerIndex),...
                    d_debug.norms.Linf(innerIndex,outerIndex)] = calcNorms(d);
                
                g_debug.images{innerIndex,outerIndex} = g;
                r_debug.images{innerIndex,outerIndex} = r;
                p_debug.images{innerIndex,outerIndex} = p;
                q_debug.images{innerIndex,outerIndex} = q;
                d_debug.images{innerIndex,outerIndex} = d;
                
                alphacg_debug(innerIndex,outerIndex) = alpha_cg;
                betacg_debug(innerIndex,outerIndex) = beta_cg;
                
            end
            
            %% Break if we violate invariants or satisfy convergence.
            if norm(g-g_old) < gradientEpsilon*norm(g)
                disp('gradient L2 converged');
                break
            end
            
            if (Linfnorm_G(innerIndex,outerIndex)/gpreNorm) < gradientEpsilon
                disp('gradient Linf converged');
                break
            end
            
            %
            %if norm(g)>norm(g_old)
            %    break
            %end
            
            
        end
        
        % Find all indices in I1 where x is >0 or where the step leads us
        % towards positive values.
        I2 = I1 & (x > 0 | d > 0);
        
        searchVector = d.*I2;
        stepSize = 1;
        xprev = x;
        
        alp = calculateStepsize(A,b,x,gouter,searchVector,stepSize);
        
        x = max(x + alp.*searchVector,0);
        
        if DEBUG
           I1_debug{outerIndex} = I1;
           I2_debug{outerIndex} = I2; 
           
           x_debug{outerIndex} = x; 
           step_debug(outerIndex) = alp;
            relNorm(outerIndex) = norm(b-A*x)./norm(b);
        end
        
        
        
        %% Convergence
        if norm(xprev-x) < xEpsilon*norm(x)
            infoLog.convergence = 'Solution converged';
            break;
        end
        
        
    end
    
    if DEBUG
       infoLog.debug.x = x_debug; 
       infoLog.debug.d = d_debug; 
       infoLog.debug.g = g_debug; 
       infoLog.debug.r = r_debug;
       infoLog.debug.p = p_debug;
       infoLog.debug.q = q_debug;
       
       
       infoLog.debug.I1 = I1_debug; 
       infoLog.debug.I2 = I2_debug; 
       
       infoLog.debug.cg_alpha = alphacg_debug;
       infoLog.debug.cg_beta = betacg_debug;
       infoLog.debug.step = step_debug;
       
       infoLog.relNorm = relNorm;
    end
    
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


function convergeInfo = calcConvergence(v)
    
end

function beta = calculateBeta(p)
    
end

function step = calculateStepSize(p)
    
end

function [L0,L1,L2,Linf] = calcNorms(vec)
   L0 = L0Norm(vec(:)); 
   L1 = L1Norm(vec(:)); 
   L2 = L2Norm(vec(:)); 
   Linf = LinfNorm(vec(:)); 
end

function n = L0Norm(vec)
   n = sum(vec(:)~=0); 
end

function n = LinfNorm(vec)
    n = max(abs(vec(:)));
end

function n = L2Norm(vec)
    n = norm(vec(:),2);
end

function n = L1Norm(vec)
    n = norm(vec(:),1);
end

