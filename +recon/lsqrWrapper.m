function [varargout] = lsqrWrapper(modelMatrix,b,varargin)

    nout = nargout;
    
    prs=inputParser;
    
    
    addParameter(prs,'Tol',1e-6);
    addParameter(prs,'Iter',20);
    addParameter(prs,'M1',[]);
    addParameter(prs,'M2',[]);
    addParameter(prs,'x0',[]);
    
    parse(prs,varargin{:});
    
%     [varargout{1:nargout}] = lsqr(modelMatrix,b,[],[],[],[],reconObj.LastImage(:));
    [varargout{1:nargout}] = lsqr(modelMatrix,b,prs.Results.Tol,prs.Results.Iter,prs.Results.M1,prs.Results.M2,prs.Results.x0);
    
    % reconObj.SolverRunningArgs{:} at the end. 
end

