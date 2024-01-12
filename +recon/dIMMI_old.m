

function [modelQ,modelTime,varargout]=dIMMI(reconObj)


    
    % If the inputs include a dump variable, write this code as it's
    % structured at runtime into the 
    if  isprop(reconObj,'dump') && reconObj.dump==true
        s = fileread([mfilename('fullpath'),'.m']);
        dump.codeString = s;
        dump.inputStruct = struct(reconObj);
        dump.inputJSON = jsonencode(struct(reconObj));
        varargout{1} = dump;
    else
        varargout{1} = [];
    end

%% Initialize all needed variables
    
    Nproj = reconObj.N_Projections; % Total number of transducers, 
    Mpts = reconObj.N_DiscreteRays; % TODO: Use this to split the dIMMI code into the constant speed of sound version and
                                       % the version that takes a
                                       % speed-of-sound image.
    time_res = reconObj.PixelOversampling;
    speedOfSound = reconObj.SpeedOfSound;

    % Transducer theta, in both circular and cylindrical coordinates.
%     xdcrtheta=linspace(inputs.signal.startAngle,...
%         inputs.signal.endAngle , Nxdcr);
    xdcrR=mean(reconObj.TransducerCoordinatesCylindrical_R);
%     [xdcrX,xdcrY]=pol2cart(xdcrtheta,xdcrR); % convert to X,Y
    xdcrXStart = reconObj.TransducerCoordinatesCartesian_X;
    xdcrYStart = reconObj.TransducerCoordinatesCartesian_Y;
    
    [theta,rho] = cart2pol(xdcrXStart,xdcrYStart);
    [thetaNew] = interp1(1:numel(theta),unwrap(theta),linspace(1,numel(theta),Nproj));
    [rhoNew] = interp1(1:numel(rho),rho,linspace(1,numel(theta),Nproj));
    [xdcrX,xdcrY] = pol2cart(thetaNew,rhoNew);
    
    
% Derive the image grid from the image itself

    % Number of pixels in the reconstructed image, as well as the span.
    % The span is from the edges of pixels, not from the center of them.
    Nx = reconObj.N_x;
    Ny = reconObj.N_y;
    
    xROI = reconObj.FieldOfView_X;
    xrange = [-xROI xROI]/2;
    
    yROI = reconObj.FieldOfView_Y;
    yrange = [-yROI yROI]/2;
    
% Supplement with 'halo' nodes.
    Nx_ghost = Nx+2; % Adds two nodes on either side of the image.
    Ny_ghost = Ny+2;
    dx = xROI./Nx;
    dy = yROI./Ny; % The dimension of the pixels. Square, in gen.
    
    %Guarantee that there are at least timeres samples present in each
    %pixel.
    dt = dx./(time_res*speedOfSound)*sqrt(2);

    % Radius of the minimum time to reach the image grid from the set of
    % transducers.
    RFOV = xROI * sqrt(2)/2;
    R_close = ( xdcrR - ((Nx+1)*dx./2)*sqrt(2) );
    R_close = ( xdcrR - RFOV );
    t_close = R_close./speedOfSound; % Convert distance to time.
    
    R_far = ( xdcrR + ((Nx+1)*dx./2)*sqrt(2) );
    R_far = ( xdcrR + RFOV );
    t_far = R_far./speedOfSound;
    
    % Get complete coverage over the ROI. This will actually mean that
    % there are slightly more than timeres samples per pixel.
    Nt = ceil((t_far-t_close)./dt);
    
    % Make a vector which spans the time range.
    tvec = linspace(t_close,t_far,Nx.*time_res);
    Nt = numel(tvec);

%% Iterate over the transducers to build the matrix.
for xdcrIndex=1:Nproj
    % Initialize the vectors to build the matrix for each transducer.
    rowIndex = [];
    colIndex = [];
    weights = [];
    
    % The ROI size and the xdcr distance define the angle of coverage
    % needed to guarantee that we're always getting the full arc.
    generalHalfAngle = asin((sqrt(2).*(Nx+1).*dx)./(2*xdcrR)); % Can move out of loop. 
    
    % Get the complex vector which points from the transducer to the center
    % of the ROI. Since the ROI center is at (0,0), just use coordinates.
    vecCent = (-xdcrX(xdcrIndex)-1i*xdcrY(xdcrIndex));
    
    % Points of tangency for the transducer to the circumcircle of the ROI.
    % These correspond to the intersection of the dual line of the
    % transducer with the circumcircle. 
    % TODO: Remove.
    firstPolar=vecCent.*(cos(generalHalfAngle)+1i*sin(generalHalfAngle));
    secondPolar=vecCent.*(cos(generalHalfAngle)-1i*sin(generalHalfAngle));
    
    % Define the vector of angles we need to propagate.
    angleVec = linspace( -generalHalfAngle , generalHalfAngle , Mpts ); % Can move out of loop. 
    
    % Normalize. 
    vecCentN = vecCent./norm(vecCent);
    
    % Iterate for every time point to build up the submatrix of the
    % transducer.
    for tIndex=1:Nt
        
        t = tvec(tIndex);
        
        % The complex location of the points, in distance, along the arc of
        % interest, with origin at the xdcr.
        % TODO: Time-distance conversion.
        aP = (speedOfSound.*t.*(vecCentN.*(cos(angleVec)+1i.*sin(angleVec)))).';
        
        
        % Convert from relative coordinates to absolute image coordinates
        xArcPt = real(aP)+xdcrX(xdcrIndex);
        yArcPt = imag(aP)+xdcrY(xdcrIndex);
        
        % Renormalize the coordinates.
        nodeMinX = xrange(1)+dx/2; % Location of the center of the first pixel
        nodeMinY = yrange(1)+dy/2;
        
        % Pixel counting coordinates, i.e. first pixel is 0,0.
        xArcPt_hat = (xArcPt + dx - nodeMinX)./dx; 
        yArcPt_hat = (yArcPt + dy - nodeMinY)./dy;
        
        
        % Gets the lower left coordinate for each pixel.
        xArcPt_L = floor(xArcPt_hat); 
        yArcPt_L = floor(yArcPt_hat);
        
        % Fractional part of the coordinate.
        xArcPt_frac = xArcPt_hat - xArcPt_L; 
        yArcPt_frac = yArcPt_hat - yArcPt_L;
        
        
        
        % Subscripts of the nodes which frame each interpolating point.
        BL_subs = [ xArcPt_L     , yArcPt_L    ]; 
        TL_subs = [ xArcPt_L     , yArcPt_L + 1];
        BR_subs = [ xArcPt_L + 1 , yArcPt_L    ];
        TR_subs = [ xArcPt_L + 1 , yArcPt_L + 1];
        
        % Bilinear weights for each node.
        BL_weights = ( 1 - xArcPt_frac ) .* ( 1 - yArcPt_frac ); 
        TL_weights = ( 1 - xArcPt_frac ) .* (     yArcPt_frac );
        BR_weights = (     xArcPt_frac ) .* ( 1 - yArcPt_frac );
        TR_weights = (     xArcPt_frac ) .* (     yArcPt_frac );
        
        
        % What is being accomplished here?? This appears to be related to
        % some form of normalization based on how big the FOV is. 
        % nodeMinX = -FOV/2, so the normalization here is basically FOV/4.
        % Done because of 4 interpolating points?
        %  IT seems that this is to improve the numerical stability of the
        %  algorithm overall. 
        
        
%         BL_weights = BL_weights.*abs(nodeMinX/2);
%         TL_weights = TL_weights.*abs(nodeMinX/2);
%         BR_weights = BR_weights.*abs(nodeMinX/2);
%         TR_weights = TR_weights.*abs(nodeMinX/2);
%         BL_weights = BL_weights.*abs(nodeMinX.*2);
%         TL_weights = TL_weights.*abs(nodeMinX.*2);
%         BR_weights = BR_weights.*abs(nodeMinX.*2);
%         TR_weights = TR_weights.*abs(nodeMinX.*2);
        
        % Trim out bad indices and bad weights, convert to linear index
        [BL_indices , BL_mask] = convertIndices( BL_subs ,Nx,Ny);
        [TL_indices , TL_mask] = convertIndices( TL_subs ,Nx,Ny);
        [BR_indices , BR_mask] = convertIndices( BR_subs ,Nx,Ny);
        [TR_indices , TR_mask] = convertIndices( TR_subs ,Nx,Ny);
        
        % Remove bad pixels.
        BL_weights = BL_weights(BL_mask);
        TL_weights = TL_weights(TL_mask);
        BR_weights = BR_weights(BR_mask);
        TR_weights = TR_weights(TR_mask);
        
        % Get the arc length over the imaging area.
        Normalizer=sum(cat(1,BL_weights,TL_weights, BR_weights, TR_weights)).*(t*speedOfSound*1000);
        Normalizer = 1./abs(nodeMinX.*2);
        
        BL_weights = BL_weights./Normalizer;
        TL_weights = TL_weights./Normalizer;
        BR_weights = BR_weights./Normalizer;
        TR_weights = TR_weights./Normalizer;

        
        
        % The indices function as row indices, while the columns are by time
        % and xdcr.
        numVals = numel([ BL_indices ; TL_indices ; BR_indices ; TR_indices ]);
        
        rowIndex(end+(1:numVals)) = (tIndex);
        colIndex(end+(1:numVals)) = [BL_indices;TL_indices;BR_indices;TR_indices];
%         weights(end+(1:numVals)) = [BL_weights;TL_weights;BR_weights;TR_weights];
        weights(end+(1:numVals)) = [BL_weights;TL_weights;BR_weights;TR_weights];
    end
    
    % Once there's a completed panel of time data, perform a derivative
    % operation in time.
    Qmat = sparse(rowIndex,colIndex,weights,Nt,Ny.*Nx);
    
    % TODO: Figure out how to make this more numerically stable. For now,
    % just multiplying dt by 1000. 
    Qdtemp = calculateDerivative(Qmat,dt); % Factor out and just loop over different panels of the original sparse matrix.
    % e.g/
%     for k = 1:Nxdcr
%         Qd((k-1)*Nt + 1:Nt) = calculateDerivative(Qmat((k-1)*Nt + 1:Nt),dt);
%     end
    
    Qdholder{xdcrIndex}=Qdtemp;
    
    
end

modelQ=vertcat(Qdholder{:});
modelTime=tvec;



end

%% Helper functions

    function derivative = calculateDerivative(x,dx)
    derivative = sparse([],[],[],size(x,1),size(x,2),nnz(x)+(2*size(x,2)));
    
    derivative(2:(end-1),:) = (x(3:end,:)-x(1:(end-2),:))./dx;
    derivative(1,:) = (sum([-1.5; 2.0; -0.5].*x(1:3,:)))./dx;
    derivative(end,:) = sum([0.5; -2.0; 1.5].*x((end-2):end,:))./dx;
    end


    function [xPole,yPole]=calculatePolarPoints(xp,yp,xc,yc,rc)
        % From the given point, given center of the circle in question, and
        % the radius, compute the points of tangency from the given point to
        % the circle in question.
        
        
        x1=(rc.^2.*(xp-xc)+rc.*(yp-yc).*sqrt((xp-xc).^2+(yp-yc).^2-rc.^2))./((xp-xc).^2+(yp-yc).^2)+xc;
        x2=(rc.^2.*(xp-xc)-rc.*(yp-yc).*sqrt((xp-xc).^2+(yp-yc).^2-rc.^2))./((xp-xc).^2+(yp-yc).^2)+xc;
        
        y1=(rc.^2.*(yp-yc)-rc.*(xp-xc).*sqrt((xp-xc).^2+(yp-yc).^2-rc.^2))./((xp-xc).^2+(yp-yc).^2)+yc;
        y2=(rc.^2.*(yp-yc)+rc.*(xp-xc).*sqrt((xp-xc).^2+(yp-yc).^2-rc.^2))./((xp-xc).^2+(yp-yc).^2)+yc;
        
        xPole=[x1;x2];
        yPole=[y1;y2];
        
        
        
    end
    
    % Converts a given set of coordinate subscripts into their
    % corresponding indices, removing the bad ones in the process.
    function [subIndices,mask]=convertIndices(subs,Nx_real,Ny_real)
        % If the subs are too low or too high, remove.
        mask = ( subs(:,1) >= 1 & subs(:,1) <= (Nx_real) ) & ...
               ( subs(:,2) >= 1 & subs(:,2) <= (Ny_real) );
        subIndices = ( subs(mask,2) + ( subs(mask,1) - 1 ) * Ny_real );
        
    end