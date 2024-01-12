function [x] = nnAPCG_reference(A,inputBundle)
    %nnAPCG reference implementation for correctness.
  

    b = inputBundle.b;
    try
    A = gpuArray(A);
    catch
    end
    %% Iterating variables
    [m,n] = size(A);
    
    maxIter_outer = 10;
    maxIter_inner = 10;
    
    g = zeros(n,1);
    d = zeros(n,1);
    p = zeros(n,1);
    
    q = zeros(m,1);
    
    try
        x = inputBundle.xo;
    catch
        x = zeros(n,1);
    end
    
    %% Monitoring variables.
    blankPanel = zeros(maxIter_inner,maxIter_outer);
    blankOuter = zeros(1,maxIter_outer);
    
    L1norm_G = blankPanel;
    L2norm_G = blankPanel;
    Linfnorm_G = blankPanel;
    
    L1norm_R = blankPanel;
    L2norm_R = blankPanel;
    Linfnorm_R = blankPanel;
    
    L1norm_D = blankPanel;
    L2norm_D = blankPanel;
    Linfnorm_D = blankPanel;
    L1norm_X = blankOuter;
    L2norm_X = blankOuter;
    Linfnorm_X = blankOuter;
    
    beta_vars = blankOuter;
    alpha_vars = blankOuter;
    
    
    N_negative = blankOuter;
    
    
    
    
    
    
    %% Outer loop.
    
    for outerIndex = 1:maxIter_outer
    
        r = b - A*x; % Residual in data space
        g_old = g; % Storage of previous gradient.
        g = - (A.'*(r)); % Projection of residual into image space.
        
        % Find all indices where x is >0 or where following the gradient
        % downwards will reduce the objective.
        %      -If x is greater than zero, it can either increase, decrease, or
        %       stay the same and will still be in the valid region, so all of
        %       those values should be optimized over.
        %      -If x is zero, it should not be able to decrease in value, but
        %       it should be allowed to either increase in value or stay at zero,
        %       so the gradient should point towards the positive half-space.
        %      - x should never be negative, because we set it to be
        %      non-negative at the end of every iteration. 
        I1 = (x > 0) | ((-g) > 0) ;
        
        d(:) = 0; % reset the step direction.
        
        %% Inner loop on the reduced set of variables. 
        for innerIndex = 1:maxIter_inner
            
            % calculate the correction coefficient.
            % Here we're looking to see the normalized projection of the
            % gradient vs. the old gradient, 
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
            g = -(A.'*r); %%%%%%%%%%
            
            % Record the data values for various parameters.
            L1norm_G(innerIndex,outerIndex) = L1Norm(g);
            L2norm_G(innerIndex,outerIndex) = L2Norm(g);
            Linfnorm_G(innerIndex,outerIndex) = LinfNorm(g);
            
            L1norm_R(innerIndex,outerIndex) = L1Norm(r);
            L2norm_R(innerIndex,outerIndex) = L2Norm(r);
            Linfnorm_R(innerIndex,outerIndex) = LinfNorm(r);
            
            L1norm_D(innerIndex,outerIndex) = L1Norm(d);
            L2norm_D(innerIndex,outerIndex) = L2Norm(d);
            Linfnorm_D(innerIndex,outerIndex) = LinfNorm(d);
            
            beta_vars(innerIndex,outerIndex) = beta_cg;
            alpha_vars(innerIndex,outerIndex) = alpha_cg;
        end
        
        % Find all indices in I1 where x is >0 or where the step leads us
        % towards positive values.
        I2 = I1 & (x > 0 | d > 0); 
        
        searchVector = d.*I2;
        stepSize = 1;
        
        x(I2) = max( (x(I2)+stepSize.*searchVector(I2)),0);
        
        
        L1norm_X(outerIndex) = L1Norm(x);
        L2norm_X(outerIndex) = L2Norm(x);
        Linfnorm_X(outerIndex) = LinfNorm(x);
        
        N_negative(outerIndex) = sum(x<0);
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