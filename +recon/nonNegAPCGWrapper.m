function [varargout] = nonNegAPCGWrapper(modelMatrix,b,varargin)

    [varargout{1:nargout}] = recon.nonNegAPCG(modelMatrix,b,varargin{:});
    
    
end

