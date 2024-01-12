function [x, infoStruc] = backproject2DSolver(backprojectionModel,inputData,varargin)
    
    % backprojectionModel should be a (NxNy x Nxdcr) array of times to
    % each pixel in the image. 
    
    if isstruct(varargin{1})&&(numel(varargin)==1)
       solverArgs = varargin{1}; % Assume this is a struct of arg pairs.
    else
        solverArgs = struct(varargin{:}); % Convert to name/value pairs
    end
    
    Nt = solverArgs.Nt;
    Nproj = solverArgs.Nproj;
    timeArray = solverArgs.timeArray;
    Fs = 1./mean(diff(timeArray));
%     
%     Nx = reconstructor.N_x;
%     Ny = reconstructor.N_y;
    
    inputData = reshape(inputData,[Nt,Nproj]);
    
%     Nsamp = Nt*Nproj;
    % Make sure that the data is de-trended.
    centeredData = inputData - mean(inputData); % apply centering matrix which acts on the columns.
    p1 = centeredData;
    
%     centMat = eye(Nt)-(1/Nt)*ones(Nt);
%     bigCent = kron(eye(Nproj),centMat);
%     p1_2 = bigCent*inputData(:);
    
    psmooth = centeredData;
%     dMat = convmtx([1;-1],Nt);
%     dMat = sparse(dMat(1:(end-1),:));
%     dMat = sparse([dMat(2:(end),:)])*Fs;
%     dMats = repmat({dMat},[Nproj,1]);
%     bigdMat = blkdiag(dMats{:});

    pderivative = [ diff(psmooth) ; zeros(1,Nproj) ] .* Fs; % apply derivative to columns
%     pderiv_2 = dMat*psmooth;
    p2 = timeArray(:) .* pderivative; 
    
%     bigOp = spdiags(repmat(timeArray(:),[Nproj,1]),0,[Nproj*Nt],[Nproj*Nt]);
%     p2_2 = bigOp * pderivative(:);
    
    
    b = p1 - p2;
%     b_2 = [speye(Nsamp) speye(Nsamp)]*[bigCent;-bigOp*bigdMat]*inputData(:);
%     superOp = [speye(Nsamp) speye(Nsamp)]*[bigCent;-bigOp*bigdMat];
%     b_2 = 
     interpt = backprojectionModel;
%    x = backprojectionModel*b(:);
     x = zeros( size(backprojectionModel,1) , 1 );
%     x_2 = zeros( size(backprojectionModel,1) , 1 );

    % For each value in interp_it, there are values in timeArray which are
    % just above and just below it.
%     [tarr,darr] = ndgrid(timeArray,1:Nproj);
%     ndarr = repmat(1:Nproj,[size(backprojectionModel,1) 1]);
% x = griddata(timeArray,1:Nproj,b,interpt
        F = griddedInterpolant( timeArray , b(:,1) , 'linear' );
    for iter = 1 : Nproj
        F.Values = b(:,iter);
        interp_it = interpt(:,iter);
        xiter = F( interp_it );
        x = x + xiter;
        
%         x = interpn(tarr,darr,b,backprojectionModel,ndarr);

        
%         lowEdge = arrayfun(@(x) find(x>=timeArray,1,'last'),interp_it);
%         highEdge = arrayfun(@(x) find(x<=timeArray,1,'first'),interp_it);
% 
%         lowDiff = abs(timeArray(lowEdge)-interp_it);
%         highDiff = abs(timeArray(highEdge)-interp_it);
%         
%         interpMat = zeros(10000,1033);
%         
%         highSub = sub2ind(size(interpMat),(1:10000)',highEdge);
%         lowSub = sub2ind(size(interpMat),(1:10000)',lowEdge);
%         interpMat(highSub) = lowDiff./(lowDiff+highDiff);
%         interpMat(lowSub) = highDiff./(lowDiff+highDiff);
%         
%         xiter_2  = interpMat*b_it;
        
    end
    
%         lowEdge = arrayfun(@(x) find(x>=timeArray,1,'last'),interpt(:));
%         highEdge = arrayfun(@(x) find(x<=timeArray,1,'first'),interpt(:));
% 
%         lowDiff = abs(timeArray(lowEdge)-interpt(:));
%         highDiff = abs(timeArray(highEdge)-interpt(:));
%         
% %         interpMat = sparse(10000,1033*50);
%         
%         pixCo = repmat((1:10000)',[1 50]);
%         transCo = repmat((1:50),[10000,1]);
%         
%         highJ = sub2ind([1033 50],highEdge,transCo(:));
%         lowJ = sub2ind([1033 50],lowEdge,transCo(:));
%         
% %         highSub = sub2ind([10000,1033,50],pixCo(:),highEdge,transCo(:));
% %         lowSub = sub2ind([10000,1033,50],pixCo(:),lowEdge,transCo(:));
% %         interpMat(highSub) = lowDiff./(lowDiff+highDiff);
% %         interpMat(lowSub) = highDiff./(lowDiff+highDiff);
%         highVal = lowDiff./(lowDiff+highDiff);
%         lowVal = highDiff./(lowDiff+highDiff);
%         
%         interpMat = sparse([pixCo(:);pixCo(:)],[highJ(:);lowJ(:)],[highVal(:);lowVal(:)],10000,1033*50);
%         tic;
% %         completeDat = interpMat*(superOp*inputData(:));
%         toc;
%         tic;
%         completeOp = interpMat*superOp;
%         toc;
%         tic;
%         testX = completeOp*inputData(:);
%         toc;
     infoStruc = ' ';
%         xnew = backprojectionModel*b(:);
    
end