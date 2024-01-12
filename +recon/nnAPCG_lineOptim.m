function [x] = nnAPCG_lineOptim(A,b,xo)
    %nnAPCG reference implementation for correctness.
    
    if nargin == 0
        load testOptimizerInput.mat;
    end
    A = A./1E6;
%     profile off;
%     profile on; 
%     A = gpuArray(A);
    %% Iterating variables
    [~,n] = size(A);
    
    maxIter_outer = 10;
    maxIter_inner = 10;
    
    g = zeros(n,1);
    d = zeros(n,1);
    p = zeros(n,1);
    
%     x = xo(:);%ones(n,1);
    x = A.'*b;
    x = max(x,0);
%     b = b;
%     xprev = x;
    
    %% Convergence variables
    gradientEpsilon = 1E-2;
    xEpsilon = 2.5E-3;
    
    %% Monitoring variables.
    blankPanel = zeros(maxIter_inner,maxIter_outer);
%     blankOuter = zeros(1,maxIter_outer);
    
%     L1norm_G = blankPanel;
%     L2norm_G = blankPanel;
    Linfnorm_G = blankPanel;
    
%     L1norm_R = blankPanel;
%     L2norm_R = blankPanel;
%     Linfnorm_R = blankPanel;
%     
%     L1norm_D = blankPanel;
%     L2norm_D = blankPanel;
%     Linfnorm_D = blankPanel;
%     L1norm_X = blankOuter;
%     L2norm_X = blankOuter;
%     Linfnorm_X = blankOuter;
%     
%     beta_vars = blankPanel;
%     alpha_vars = blankPanel;
    
    
%     N_negative = blankOuter;
    
    
  
    
    
    
    %% Outer loop.
    
    for outerIndex = 1:maxIter_outer
    
          r = b - A*x; 
        g_old = g;
%         g = -1.*((r.')*A).';
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
%             g = -(A.'*r); %%%%%%%%%%
%             g = -1.*((r.')*A).';
            g = -(A.'*r);
            
            % Record the data values for various parameters.
%             L1norm_G(innerIndex,outerIndex) = L1Norm(g);
%             L2norm_G(innerIndex,outerIndex) = L2Norm(g);
            Linfnorm_G(innerIndex,outerIndex) = LinfNorm(g);
%             
%             L1norm_R(innerIndex,outerIndex) = L1Norm(r);
%             L2norm_R(innerIndex,outerIndex) = L2Norm(r);
%             Linfnorm_R(innerIndex,outerIndex) = LinfNorm(r);
%             
%             L1norm_D(innerIndex,outerIndex) = L1Norm(d);
%             L2norm_D(innerIndex,outerIndex) = L2Norm(d);
%             Linfnorm_D(innerIndex,outerIndex) = LinfNorm(d);
%             
%             beta_vars(innerIndex,outerIndex) = beta_cg;
%             alpha_vars(innerIndex,outerIndex) = alpha_cg;
%             
%             
%             figure(1);
%             subplot(121);
%             imagesc(reshape(g,[1 1]*sqrt(numel(g)))); title('g');
%             subplot(122);
%             imagesc(reshape(d,[1 1]*sqrt(numel(d)))); title('d');
            
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
            if norm(g)>norm(g_old)
                disp('inner gradient norm increased')
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
        
        
%         originalObjective = x.'*A.'*A*x - x.'*A.'*b -b.'*A*x +b.'*b;
%         anaObjective = (A*x-b).'*(A*x-b);
%         assert((originalObjective-anaObjective)./((originalObjective+anaObjective)/2)<sqrt(eps))
%         
%         linearCorrection = (A*x).'*(A*psearch)  +  (A*psearch).'*(A*x)  -  (A*psearch).'*(b)  -  (b).'*(A*psearch);
%         quadraticCorrection = (A*psearch).'*(A*psearch);
%         
%         objectiveRaw = @(A,x,b,ps,alpha) (A*(x+alpha*ps)-b).'*(A*(x+alpha*ps)-b);
%         objectiveExpand = @(alpha) originalObjective + alpha*(linearCorrection) + (alpha.^2)*(quadraticCorrection);
        
        
%         xcopy = x;
%         x(I2) = max( (x(I2)+alp.*searchVector(I2)),0);
        x = max(x + alp.*searchVector,0);
        
        
        
%             figure(2);
%             subplot(131);
% %             imagesc(reshape(x,[1 1]*sqrt(numel(d)))); title('d');
%             subplot(132);
% %             imagesc(reshape(alp.*searchVector(I2),[1 1]*sqrt(numel(alp.*searchVector(I2))))); title('update');
%             subplot(133);
%             imagesc(reshape(x,[1 1]*sqrt(numel(x)))); title('x');
        
%         L1norm_X(outerIndex) = L1Norm(x);
%         L2norm_X(outerIndex) = L2Norm(x);
%         Linfnorm_X(outerIndex) = LinfNorm(x);
        
%         N_negative(outerIndex) = sum(x<0);
        
        %% Convergence
        if norm(xprev-x) < xEpsilon*norm(x)
           disp('x converged');
%            profile viewer;
           return; 
        end
        
        
    end
%     profile viewer;
%     s = 0;
    
        function armijoAlpha = calculateStepsize(A,b,x,g,pdir,alpha)
            tau = 0.5;
            c = 0.1;
            
            Ax = A*x;
            Ap = A*pdir;
            
%             anaObj = @(a) (A*(x+a*pdir)-b).'*(A*(x+a*pdir)-b);
            
%             baseObj = dot(Ax,Ax) - dot(Ax,b) - dot(b,Ax) + dot(b,b);
            linObj  = dot(Ax,Ap) + dot(Ap,Ax) - dot(Ap,b) - dot(b,Ap);
            quadObj = dot(Ap,Ap);
            
            mslope = dot(pdir,g);
            
            expectedLoss = -mslope*c;
            
%             calcLossAna = anaObj(0)-anaObj(alpha);
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
    

end

%% Helper functions

function n = LinfNorm(vec)
   n = max(vec(:));
end

function n = L2Norm(vec)
   n = norm(vec(:),2);
end

function n = L1Norm(vec)
   n = norm(vec(:),1);
end