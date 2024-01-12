function [sos,results,varargout] = backProjectTune(fileLoader,dataMeta,Npix,FOV)

    
    [scanSize.Nrun,scanSize.Nz,scanSize.Nrep,scanSize.Nwl,scanSize.Nshot] = size(dataMeta.ScanStructure);
    scanSize.Dims = size(dataMeta.ScanStructure);
    

    %% Build the coordinate system of the data itself.
     Fs = dataMeta.HWDesc.SamplingFrequency;
     Nsamples = double(dataMeta.MeasurementDesc.RecordLength);

     Nxdcr = double(dataMeta.HWDesc.NumDetectors);

    dt = 1 / Fs;
    dataTime = dt * (1:Nsamples);
    xdcrIndex = 1 : Nxdcr;
    xdcrCoordsX = dataMeta.HWDesc.TransducerCoordinatesCartesian_X;
    xdcrCoordsY = dataMeta.HWDesc.TransducerCoordinatesCartesian_Y;
    

%% Build the coordinate system of the target image space. 
    viewRange = [-FOV,FOV]/2;
    
    pixelCorners = linspace(viewRange(1),viewRange(2),Npix+1);
    [dx, dy] = deal(pixelCorners(2)-pixelCorners(1));
    pixelCenters = pixelCorners(1:(end-1)) + dx/2;
    
    [pixelX, pixelY] = meshgrid(pixelCenters);
    
    

    %% Construct the backprojection mapping (parameterized by SoS)
    xDistance = pixelX(:) - (xdcrCoordsX(:).');
    yDistance = pixelY(:) - (xdcrCoordsY(:).');
    rDistance = hypot(xDistance,yDistance);





    %% Loop over a reasonable range of SoS and create a vector of focus scores to manually find the minimum of.
    [~,dataSample]=min(abs([dataMeta.ScanFrames(:).Wavelength]-800));
    
    testData = fileLoader(dataSample,{});
    testSignal = testData.Data;
    testSignal = movmean(testSignal,3);
    
    % Dist = speed * time
    

    sos = 1480:1580;
    nSpeeds = numel(sos);

    [sobelScoreVec, consistentGradientScoreVec,giniVal] = deal(zeros(nSpeeds,1));
    h=waitbar(0, 'testing sound speed 1480-1580');
    for k = 1:nSpeeds
        apoWinStart = (dataMeta.HWDesc.Radius - FOV/2)./sos(k);
        apoWinEnd = (dataMeta.HWDesc.Radius + FOV/2)./sos(k);
        apoWin = ones(size(testSignal));
%         timeStart = find(dataTime<apoWinStart,1,'last');
%         timeEnd = find(dataTime>apoWinEnd,1,'first');
%         Nt = timeEnd-timeStart+1;
%         L = tukeywin(Nt,0.25);
%         apoWin(timeStart:timeEnd,:) = repmat(L,[1 size(testSignal,2)]);
        
        
        apoSig = testSignal.*apoWin;
%         bpImageIter = getBackprojectionImage(apoSig,rDistance,sos(k)+1,dataMeta.HWDesc.SamplingFrequency,dataTime);
        bpImageIter2 = getBackprojectionImage(apoSig,rDistance,sos(k),dataMeta.HWDesc.SamplingFrequency,dataTime);
%         bpImageIter3 = getBackprojectionImage(apoSig,rDistance,sos(k)-1,dataMeta.HWDesc.SamplingFrequency,dataTime);
       testImage = bpImageIter2;
%        testImage = ((bpImageIter -2.*bpImageIter2 + bpImageIter3).^2);
%        testImage = var(cat(3,bpImageIter,bpImageIter2,bpImageIter3),[],3);
%        testImage = (bpImageIter+bpImageIter2-bpImageIter3).^2;
%        testImage = mean(cat(3,bpImageIter,bpImageIter2,bpImageIter3),3);
%        testImage = testImage./sum(testImage(:));
       %figure(1);imagesc(testImage); title(num2str(sos(k))); pause(0.025);
       sobelScoreVec(k) = getFocusScore_sobelVar(testImage);
       consistentGradientScoreVec(k) = getFocusScore_consistentGradient(testImage,0.25);
       giniVal(k) = gini(ones(numel(testImage),1),abs(testImage(:)));
%        giniVal(k) = sum(testImage(:));
       hists(k,:) = histcounts(testImage(:),512,'Normalization','probability');
%        consistentGradientScoreVec(k) = norm(bpImageIter+bpImageIter2);
saveIm{k} = bpImageIter2;
waitbar(k/nSpeeds, h, sprintf('testing sound speed of %d%',sos(k)))
    end
    delete(h)




    sos_SobelVar = sos(sobelScoreVec == min(sobelScoreVec(:)));
    sos_ConsistentGradient = sos(consistentGradientScoreVec == min(consistentGradientScoreVec(:)));

    
    
    results.SosRange = sos;
    results.SobelScores = sobelScoreVec;
    results.ConsistentGradientScores = consistentGradientScoreVec;
    results.SoSSobel = sos_SobelVar;
    results.SoSConsistentGradient = sos_ConsistentGradient;
    results.gini = giniVal;
    results.hist = hists;
    
    sos = mean(sos_SobelVar,sos_ConsistentGradient);

    if nargout >2
        varargout{1} = saveIm;
    end
end













%% Helper functions

% Perform backprojection given structural parameters.

function bpImage = getBackprojectionImage(data,rDistance,speedOfSound,samplingFrequency,samplingTimes)
   
    [Npix,Nxdcr] = size(rDistance);
    
    centeredData = abs(hilbert(data - mean(data)));
    centeredData = data - mean(data);
    p1 = centeredData;
    
    psmooth = p1;
            pderivative = [ diff(psmooth) ; zeros(1,Nxdcr) ] .* samplingFrequency;
            p2 = samplingTimes(:).*pderivative;
            
    b = p1-p2;
    
    q = zeros(Npix,1);
    
    for transducerID = 1:Nxdcr
        q(:) = q(:) + interp1(samplingTimes(:),b(:,transducerID),rDistance(:,transducerID)./speedOfSound);
        
    end
    
    bpImage = reshape(q,[sqrt(Npix),sqrt(Npix)]);
end



% Calculation of the focus score(s)

function focusScore = getFocusScore_sobelVar(imageInput)
    
    gradientImage = abs(imfilter(imageInput,fspecial('sobel')));
    meanG = mean(gradientImage(:));
    enerG = sum(gradientImage(:).^2);
    Npix = numel(imageInput);
    
    focusScore = -1./(enerG).*sum((gradientImage(:)-meanG).^2);
    
    
end




function focusScore = getFocusScore_consistentGradient(imageInput,horizontalEdgeWeight)
    
    [Ny,Nx] = size(imageInput);
%     W = tukeywin(Ny,0.1);
%     imageInput = abs(imageInput).*(W.*W.');
    
    initialGradient = imfilter(fspecial('sobel'),imageInput);
    absGradient = abs(initialGradient);
    cannyWeight = prctile(absGradient(:),95);
    
    diffusedImage = imageInput;
    
    
        dI_north = zeros(Ny,Nx);
        dI_south = zeros(Ny,Nx);
        dI_east  = zeros(Ny,Nx);
        dI_west  = zeros(Ny,Nx);
        lambda = 0.5;
    
   
        
    for timeIt = 2:35
        
        prevImage = diffusedImage;
        
        diffCoeff = @(X)  1./( 1 + (abs(X)./cannyWeight).^2);
        
        dI_north(2:(end-1),2:(end-1))  =  diffusedImage(1:(end-2),2:(end-1)) - diffusedImage(2:(end-1),2:(end-1));
        dI_south(2:(end-1),2:(end-1))  =  diffusedImage(3:( end ),2:(end-1)) - diffusedImage(2:(end-1),2:(end-1));
        dI_east(2:(end-1),2:(end-1))   =  diffusedImage(1:(end-2),1:(end-2)) - diffusedImage(2:(end-1),2:(end-1));
        dI_west(2:(end-1),2:(end-1))   =  diffusedImage(1:(end-2),3:( end )) - diffusedImage(2:(end-1),2:(end-1));
        
        conduction_north = diffCoeff(dI_north) ;
        conduction_south = diffCoeff(dI_south) ;
        conduction_east  = diffCoeff(dI_east ) ;
        conduction_west  = diffCoeff(dI_west ) ;
        
        diffusedImage = diffusedImage + lambda*(conduction_north.*dI_north +...
                                                conduction_south.*dI_south +...
                                                conduction_east.*dI_east +...
                                                conduction_west.*dI_west );
                                            
        nrm(timeIt) = norm(prevImage(:)-diffusedImage(:));
        
        if abs(nrm(timeIt)-nrm(timeIt-1))./norm(diffusedImage(:)) < 0.01
           break; 
        end
        
    end
    
    
    enerG = sum(imageInput(:).^2);
    
    CG_Operator = [-0.003776, -0.010199, 0, 0.010199, 0.003776;
                   -0.026786, -0.070844, 0, 0.070844, 0.026786;
                   -0.046548, -0.122572, 0, 0.122572, 0.046548;
                   -0.026786, -0.070844, 0, 0.070844, 0.026786;
                   -0.003776, -0.010199, 0, 0.010199, 0.003776];
               
   edges_Horizontal = conv2(diffusedImage,CG_Operator,'same');
   edges_Vertical = conv2(diffusedImage,CG_Operator.','same');
   
   focusScore = -1./(Ny*Nx) .* sum( horizontalEdgeWeight .* edges_Horizontal(:).^2 + ...
                                (1- horizontalEdgeWeight).* edges_Vertical(:).^2);
   

    
    
end

function [g,l,a] = gini(pop,val,makeplot)
    % check arguments
    assert(nargin >= 2, 'gini expects at least two arguments.')
    if nargin < 3
        makeplot = false;
    end
    assert(numel(pop) == numel(val), ...
        'gini expects two equally long vectors (%d ~= %d).', ...
        size(pop,1),size(val,1))
    pop = [0;pop(:)]; val = [0;val(:)];     % pre-append a zero
    isok = all(~isnan([pop,val]'))';        % filter out NaNs
    if sum(isok) < 2
        warning('gini:lacking_data','not enough data');
        g = NaN; l = NaN(1,4);
        return;
    end
    pop = pop(isok); val = val(isok);
    
    assert(all(pop>=0) && all(val>=0), ...
        'gini expects nonnegative vectors (neg elements in pop = %d, in val = %d).', ...
        sum(pop<0),sum(val<0))
    
    % process input
    z = val .* pop;
    [~,ord] = sort(val);
    pop    = pop(ord);     z    = z(ord);
    pop    = cumsum(pop);  z    = cumsum(z);
    relpop = pop/pop(end); relz = z/z(end);
    
    % Gini coefficient
    % We compute the area below the Lorentz curve. We do this by
    % computing the average of the left and right Riemann-like sums.
    % (I say Riemann-'like' because we evaluate not on a uniform grid, but
    % on the points given by the pop data).
    %
    % These are the two Rieman-like sums:
    %    leftsum  = sum(relz(1:end-1) .* diff(relpop));
    %    rightsum = sum(relz(2:end)   .* diff(relpop));
    % The Gini coefficient is one minus twice the average of leftsum and
    % rightsum. We can put all of this into one line.
    g = 1 - sum((relz(1:end-1)+relz(2:end)) .* diff(relpop));
    
    % Lorentz curve
    l = [relpop,relz];
    a = [pop,z];
    if makeplot   % ... plot it?
        area(relpop,relz,'FaceColor',[0.5,0.5,1.0]);    % the Lorentz curve
        hold on
        plot([0,1],[0,1],'--k');                        % 45 degree line
        axis tight      % ranges of abscissa and ordinate are by definition exactly [0,1]
        axis square     % both axes should be equally long
        set(gca,'XTick',get(gca,'YTick'))   % ensure equal ticking
        set(gca,'Layer','top');             % grid above the shaded area
        grid on;
        title(['\bfGini coefficient = ',num2str(g)]);
        xlabel('share of population');
        ylabel('share of value');
    end
    
end