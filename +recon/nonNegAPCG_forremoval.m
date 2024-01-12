function [x,infostruc] = nonNegAPCG(A,b,varargin)
    %ACPG.m ACCELERATED PROJECTED CONJUGATE GRADIENT
    %Given a model matrix A, an initial guess x, and data b, finds the
    %non-negative solution of x to argmin((Ax-b)^2)
%     persistent xi;

    if nargin == 0
        load testOptimizerInput.mat;
    end
    
    A = A./1E6;
    varargin = {};
    %% Preconditioning checks
    prs=inputParser;
    
    %Validation functions
    boundOnUnit=@(zzz) ((zzz>=0)&(zzz<=1));
    
    %Require that the inputs are numeric.
    addRequired(prs,'A',@(x) validateattributes(x,{'numeric' 'sparse'},{'2d'}));
    
    
    addRequired(prs,'b',@(x) validateattributes(x,{'numeric' 'sparse'},{'2d','column'}));
    
    
    
    
    %Set defaults that are passed to the function.
    
    addParameter(prs,'NMaxIterOuter',10);
    addParameter(prs,'NMaxIterInner',10);
    addParameter(prs,'GradientTol',1E-4);
    addParameter(prs,'ConvergeTol',1E-6);
    addParameter(prs,'c_ls',0.1,boundOnUnit);
    addParameter(prs,'tau_ls',0.5,boundOnUnit);
    addParameter(prs,'nMax_ls',100,@isinteger);
    addParameter(prs,'a_step',10);
    addParameter(prs,'verbose',false);
    addParameter(prs,'regularized',false);
    addParameter(prs,'lambda',0.01);
    addParameter(prs,'x0',zeros(size(A,2),1));
    
    parse(prs,A,b,varargin{:});
    
    goodMatMultiply=(size(A,1)==size(b,1));%&&size(A,2)==size(x0,1));
    
    
    %Sizes of A,x, and b compatible.
    if ~goodMatMultiply
        error('Size of inputs do not match.');
    end
    
    
    %%%%%%%%%%%%%%%
    
    % If regularization is on, construct the kernel and append.
    if prs.Results.regularized
       kern=(1/9)*[-1 -1 -1;
                   -1  8 -1;
                   -1 -1 -1]; %high pass filter
       L=prs.Results.lambda.*convmtx2(kern,sqrt(size(A,2)),sqrt(size(A,2)));
       A=[A;L];
       b=[b;zeros(size(L,1),1)];
    else
%         A=A;
    end
    
    
    %% Pre-allocation of storage bits.
    
    % Matrix A is of size m x n
    [m,n]=size(A);
    % xi is of size n x 1
%     xi=abs(lsqr(A,b,[],5));
    xi=zeros(size(A,2),1);
    
    %% Get the initial guess, slightly regularized.

    d=zeros(n,1);
    g=zeros(n,1);
    q=zeros(m,1);
    r=zeros(m,1);
    p=d;
    
    %%Outer CG loop
    for k=1:prs.Results.NMaxIterOuter
        
        % r=b-Ax
        r(:)=b(:)-(A*xi);%%%%
        
        %g_old<-g
        g_old=g;%%%%
        
        %g<- -A.'*r
        g=-((r.')*A).';%%%%
            normG_L2_0=norm(g);
                normG_L1_0=sum(abs(g(:)));
        
        
        I1=logical(((xi>0)|((-g)>0)));%%%%
        d(:)=0;%%%%
        %% Inner CG loop
        
        for it=1:prs.Results.NMaxIterInner
            
            beta_cg=max(0, ...
                -( (g(I1).')*(g(I1) - g_old(I1)))./((g_old(I1).')*g_old(I1)));
            
            p(:)= (-g + beta_cg.*p).*I1;%%%%
%             p(~I1)=0;%%%%
            
            q=A*(p);%%%%
            alpha_cg=(-((g.')*(p)))./((q.')*q);%%%%
            
            d(I1)=d(I1)+alpha_cg.*p(I1);%%%%
            
            r=r-(alpha_cg*q);%%%%
            
            g_old=g; %%%%
            g=-((A.')*r); %%%%
            
            

                normG_L2=norm(g);
                normG_L1=sum(abs(g(:)));
                innerConverge_L2=normG_L2./normG_L2_0;
                innerConverge_L1=normG_L1./normG_L1_0;
                
            
            if innerConverge_L2<prs.Results.GradientTol
                disp('Inner L2 convergence!');
                infostruc.converge(it)='L2';
                break;
            elseif innerConverge_L1<prs.Results.GradientTol
                disp('Inner L1 convergence!');
                infostruc.converge(it)='L1';
                break;
            end
                
            
            
            % Convergence criterion
            
%             %% loop termination
%             if it>1
%                 dotvec=(g.'*d)./(norm(d));
%                 gg=g'*(g.*I1);
%                 
%                 
%                 if gg<=(prs.Results.GradientTol)
%                     disp(['Gradient Convergence' num2str(-dotvec./gg)]);
%                     badInner=0;%%%%%%%set this back to 1
%                     break;
%                 end
%                 
%                 if dotvec>0
%                     disp('Inner non-descent!')
%                     it;
%                     d=-g;
%                     continue;
%                 end
%                 
%                 
%                 %% Full parameters
%                 aaa_angle=(-g'*d)./(norm(g)*norm(d));
%                 
%                 
%                 if dotvec<=(-inner_descent_tol*gg)
%                     disp(['Descent break: ' num2str(-dotvec./gg)]);
%                     badInner=0;
%                     break;
%                 end
%                 
%                 if aaa_angle>=inner_angle_tol
%                     disp(['Angle break: ' num2str(aaa_angle)]);
%                     badInner=0;
%                     break;
%                 end
%             end
            
            
            
        end
        %%
        
        
        
        I2=( ((xi.*I1) > 0) | (d.*I1 > 0)); %Check if lagrange multipliers are positive.
        
        
        
        
        a_step=1;%prs.Results.a_step;

        %% Update and project
        
        x_last=xi;
        xi=xi.*I2+a_step.*(d.*I2);
        xi(xi<0)=0;
        
%         imagesc(reshape(xi,[sqrt(numel(xi)) sqrt(numel(xi))])); colormap gray; drawnow;
        
        if k==1
            projNorm0=max(r(r>0));
            norm_ratio=Inf;
        else
            projNormCur=max(r(r>0));
            norm_ratio=projNormCur./projNorm0;
        end
        
        
        
        
        if norm_ratio<eps
            disp('Convergent in Inf norm');
            break;
        end
    end
    infostruc.nIter=k;
    x=xi;
end

