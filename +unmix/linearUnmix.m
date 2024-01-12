function [unmixOutput,unmixLog] = linearUnmix(mixModel,msImage, varargin)
    %linearUnmix performs simple unmixing based on the given mixing model.
    % This is the 'most basic' unmixing solution and will not, in general,
    % provide physically-meaningful results for all pixels. 
    % 
    % The mixModel is assumed to be a single N_components x N_wavelengths matrix
    % describing the mixture of the multiple components into a single
    % multispectral signature for a given pixel. 
    
%     [Ny,Nx,Nwl]= size(msImage);
    [Nwl,Nc] = size(mixModel);
    Nel = numel(msImage);
    Npix = Nel/Nwl;
    
    tempIm = reshape(msImage,[Npix,Nwl]);
    
    unmixTemp = mrdivide(tempIm,mixModel.');
    
    unmixOutput = unmixTemp(:);
    unmixLog = '';
    
end

