


function [unmixedOutput,unmixLog] = nnlsUnmix(msImage,mixModel,varargin)
% nnlsUnmix applies a constrained linear solver to the unmixing process.
% The model is assumed to be a linear mixture model, expanded as the sparse
% Kronecker product of the original mixing matrix and an identity matrix of size
% (Ny*Nx)x(Ny*Nx). 
%
% This implementation should generally be taken as a reference point; it is not
% particularly efficient.
    
    [Ny,Nx,Nwl]= size(msImage);
    [NW_,NC] = size(mixModel);
    
    Nwl_ = NW_/(Ny*Nx);
    Nc = NC/(Ny*Nx);
    
    if ( (Nwl~=Nwl_) && (Nc==Nwl) ) 
        mixModel = mixModel.';
       [Nc,Nwl_] = size(mixModel);
    end
    
    tempIm = reshape(permute(msImage,[3 1 2]),[Nwl,Ny*Nx]);
    
     output = lsqlin(mixModel,tempIm(:),[],[],[],[],zeros(numel(tempIm)/Nwl*Nc,1),[]);
    outputStitch = reshape(output,[Nc Ny Nx]);

    unmixedOutput = permute(outputStitch,[2 3 1]);
    unmixLog = '';
end