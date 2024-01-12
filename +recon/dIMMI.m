

function [modelQ,varargout] = dIMMI(input_coords,output_coords,varargin)
% function [modelQ,varargout]=dIMMI(input_coords,output_coords,model_args,reconObj,varargin)
% function [Qd, modelTimeList,varargout] = dIMMI(input_coords,output_coords,model_args)
%     varargout = {};
if isstruct(varargin{1})&&(numel(varargin)==1)
        model_args = varargin{1}; % Assume this is a struct of arg pairs.
    else
        model_args = struct(varargin{:}); % Convert to name/value pairs
    end
    
    varargout{1} = [];
    
    if isa(output_coords,'frame.CoordinateSystem')
        tvec = output_coords.getCoordinateValues('Time'); % Directly returns a list of time values. TODO: Add unit selection a la .getCoordinateValues('Time','ms');
        Nt = size(tvec,1);
        xdcrvec = output_coords.getCoordinateValues('Transducer'); % Should return a list of Coordinates.
        xdcrX = xdcrvec.getCoordinateValues('EUCLIDEAN_X'); % xdcrvec.getCoordinateValues('EUCLIDEAN_X');
        xdcrY = xdcrvec.getCoordinateValues('EUCLIDEAN_Y');
        Nxdcr = numel(xdcrX);
        
    else
        tvec = output_coords{1};
        Nt = size(tvec,1);
        xdcr = output_coords{2};
        Nxdcr = size(xdcr,1);
        xdcrX = xdcr(:,1);
        xdcrY = xdcr(:,2);
    end
    
    if isa(input_coords,'frame.CoordinateSystem')
        ypix = input_coords.getCoordinateValues('EUCLIDEAN_Y');
        Ny = numel(ypix);
        dy = mean(diff(ypix));
        nodeMinY = ypix(1);
        nodeMaxY = ypix(end);
        
        xpix = input_coords.getCoordinateValues('EUCLIDEAN_X');
        Nx = numel(xpix);
        dx = mean(diff(xpix));
        nodeMinX = xpix(1);
        nodeMaxX = xpix(end);
        
    else
        ypix = input_coords{1};
        Ny = numel(ypix);
        dy = mean(diff(ypix));
        nodeMinY = ypix(1);
        nodeMaxY = ypix(end);
        
        xpix = input_coords{2};
        Nx = numel(xpix);
        dx = mean(diff(xpix));
        nodeMinX = xpix(1);
        nodeMaxX = xpix(end);
    end
    dt = mean(diff(tvec));
    fieldOfView = range(xpix)+dx;
    
    try % TODO: Factor this out.
    N_DiscreteRays = model_args.NDiscreteRays; 
    catch
    N_DiscreteRays = model_args.N_DiscreteRays; 
    end
    speedOfSound = model_args.SpeedOfSound;
%     try 
%         speedOfSound = reconObj.SpeedOfSound;
%     catch
%     end
    
    
    % If the inputs include a dump variable, write this code as it's
    % structured at runtime. TODO: Try to move this up to ReconSystem. 
%     if  isprop(reconObj,'dump') && reconObj.dump==true
%         s = fileread([mfilename('fullpath'),'.m']);
%         dump.codeString = s;
%         dump.inputStruct = struct(reconObj);
%         dump.inputJSON = jsonencode(struct(reconObj));
%         varargout{1} = dump;
%     else
%         varargout{1} = [];
%     end

%% Initialize all needed variables
%     
%     Nproj = reconObj.N_Projections; % Total number of transducers, 
% %     Mpts = reconObj.N_DiscreteRays; 
%     
%     time_res = reconObj.PixelOversampling;
% %     speedOfSound = reconObj.SpeedOfSound;
% 
%     % Transducer theta, in both circular and cylindrical coordinates.
%     xdcrR=mean(reconObj.TransducerCoordinatesCylindrical_R);
%     xdcrXStart = reconObj.TransducerCoordinatesCartesian_X;
%     xdcrYStart = reconObj.TransducerCoordinatesCartesian_Y;
%     
%     [theta,rho] = cart2pol(xdcrXStart,xdcrYStart);
%     [thetaNew] = interp1(1:numel(theta),unwrap(theta),linspace(1,numel(theta),Nproj));
%     [rhoNew] = interp1(1:numel(rho),rho,linspace(1,numel(theta),Nproj));
%     [xdcrX,xdcrY] = pol2cart(thetaNew,rhoNew);
%     
%     
% % Derive the image grid from the image itself
% 
%     % Number of pixels in the reconstructed image, as well as the span.
%     % The span is from the edges of pixels, not from the center of them.
%     Nx = reconObj.N_x;
%     Ny = reconObj.N_y;
%     
%     xROI = reconObj.FieldOfView_X;
%     xrange = [-xROI xROI]/2;
%     
%     yROI = reconObj.FieldOfView_Y;
%     yrange = [-yROI yROI]/2;
%     
% % Supplement with 'halo' nodes.
%     dx = xROI./Nx;
%     dy = yROI./Ny; % The dimension of the pixels. Square, in gen.
%     
%     %Guarantee that there are at least timeres samples present in each
%     %pixel.
%     dt = dx./(time_res*speedOfSound)*sqrt(2);
% 
%     % Radius of the minimum time to reach the image grid from the set of
%     % transducers.
%     RFOV = xROI * sqrt(2)/2;
%     R_close = ( xdcrR - ((Nx+1)*dx./2)*sqrt(2) );
%     R_close = ( xdcrR - RFOV );
%     t_close = R_close./speedOfSound; % Convert distance to time.
%     
%     R_far = ( xdcrR + ((Nx+1)*dx./2)*sqrt(2) );
%     R_far = ( xdcrR + RFOV );
%     t_far = R_far./speedOfSound;
%     
%     % Get complete coverage over the ROI. This will actually mean that
%     % there are slightly more than timeres samples per pixel.
%     Nt = ceil((t_far-t_close)./dt);
    
    % Make a vector which spans the time range.
%     tvec = linspace(t_close,t_far,Nx.*time_res);
%     Nt = numel(tvec);
%     Nt = size(tvec,1);

DEBUG = false;
if DEBUG
    d_Image_Data = zeros(Ny,Nx);
        d_Fig = figure;
        d_Ax = axes;
        d_Image = imagesc(d_Ax,d_Image_Data,'XData',[min(xpix) max(xpix)],'YData',[min(ypix) max(ypix)]); hold on;
        d_Ax.YDir='normal';
        xdcrScat = scatter(d_Ax,xdcrX(:),xdcrY(:));
        
        q_mid= quiver(d_Ax,xdcrX(1),xdcrY(1),-0.01,-0.01,'AutoScale','off');
        
        scat_mid=scatter(0,0,'r');
        scat_int=scatter(0,0,'m');
        
        axis([min(xdcrX) max(xdcrX) min(xdcrY) max(xdcrY)]);
        pbaspect([1 1 1]);
        
        miX = min(xpix)-dx/2;
        maX = max(xpix)+dx/2;
        miY = min(ypix)-dy/2;
        maY = max(ypix)+dy/2;
        
        plot([miX,miX],[miY,maY],'k');
        plot([maX,maX],[miY,maY],'k');
        plot([miX,maX],[maY,maY],'k');
        plot([miX,maX],[miY,miY],'k');
        
        h_circ = viscircles([xdcrX(1),xdcrY(1)],0.02);
        
        for d_k = 1:(Ny-1)
            plot([miX,miX]+dx.*d_k,[miY,maY],'k');
            plot([miX,maX],[miY,miY]+dy.*d_k,'k');
        end
    
end


%% Iterate over the transducers to build the matrix.
for xdcrIndex = 1:Nxdcr
    
    % Coordinate determination
    thisX = xdcrX(xdcrIndex);
    thisY = xdcrY(xdcrIndex);
    
    % Initialize the vectors to build the matrix for each transducer.
    rowIndex = [];
    colIndex = [];
    weights = [];
    
    % The ROI size and the xdcr distance define the angle of coverage
    % needed to guarantee that we're always getting the full arc.
    generalHalfAngle = asin((sqrt(2).*(Nx+1).*dx)./(2*hypot(thisX,thisY))); % Can move out of loop. 
    angleVec = linspace( -generalHalfAngle , generalHalfAngle , N_DiscreteRays ); % Can move out of loop. 
    incAngle = 2.*generalHalfAngle;
    % Get the complex vector which points from the transducer to the center
    % of the ROI. Since the ROI center is at (0,0), just use coordinates.
    vecCent = ( -thisX - 1i*thisY);
    
    % Points of tangency for the transducer to the circumcircle of the ROI.
    % These correspond to the intersection of the dual line of the
    % transducer with the circumcircle. 
    % TODO: Remove.
%     firstPolar=vecCent.*(cos(generalHalfAngle)+1i*sin(generalHalfAngle));
%     secondPolar=vecCent.*(cos(generalHalfAngle)-1i*sin(generalHalfAngle));
    
    % Define the vector of angles we need to propagate.
    
    % Normalize. 
    vecCentN = vecCent./norm(vecCent);
    
    % Iterate for every time point to build up the submatrix of the
    % transducer.
    for tIndex=1:Nt
        
        t = tvec(tIndex);
        
        % The complex location of the points, in distance, along the arc of
        % interest, with origin at the xdcr.
        R = speedOfSound.*t;
        aP = (speedOfSound.*t.*(vecCentN.*(cos(angleVec)+1i.*sin(angleVec)))).';
        
        
        % Convert from relative coordinates to absolute image coordinates
        xArcPt = real(aP) + thisX;
        yArcPt = imag(aP) + thisY;
        
        isInsideCircle = (xArcPt>=(nodeMinX-dx/2))&...
                         (xArcPt<=(nodeMaxX+dx/2))&...
                         (yArcPt>=(nodeMinY-dy/2))&...
                         (yArcPt<=(nodeMaxY+dy/2));
        dL = incAngle./N_DiscreteRays;
        % Renormalize the coordinates.
%         nodeMinX = xrange(1)+dx/2; % Location of the center of the first pixel
%         nodeMinY = yrange(1)+dy/2;
        
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
            
        % The indices function as row indices, while the columns are by time
        % and xdcr.
        numVals = numel([ BL_indices ; TL_indices ; BR_indices ; TR_indices ]);
        
        % Get the arc length over the imaging area.
        Normalizer = (4*pi.*speedOfSound);
%         Normalizer = (4*pi.*speedOfSound.*N_DiscreteRays);
%         Normalizer = ((1*pi)*speedOfSound.*(incAngle/(2*pi)).*speedOfSound); % One speed of sound to make R, the other to handle normalization. 
        
        BL_weights = BL_weights./Normalizer;
        TL_weights = TL_weights./Normalizer;
        BR_weights = BR_weights./Normalizer;
        TR_weights = TR_weights./Normalizer;

        
    
        rowIndex(end+(1:numVals)) = (tIndex);
        colIndex(end+(1:numVals)) = [BL_indices;TL_indices;BR_indices;TR_indices];
        weights(end+(1:numVals)) = [BL_weights;TL_weights;BR_weights;TR_weights].*dL;
        
        if DEBUG && numVals>0
%                 q_mid.Visible = 'off';
                q_mid.UData = real(aP);
                q_mid.VData = imag(aP);
                q_mid.XData = ones(size(q_mid.UData)).*thisX;
                q_mid.YData = ones(size(q_mid.VData)).*thisY;
                
                
                scat_mid.XData = q_mid.XData+q_mid.UData;
                scat_mid.YData = q_mid.YData+q_mid.VData;
                
                d_Image.CData(:)=0;
                d_Image.CData(:) = accumarray([BL_indices;TL_indices;BR_indices;TR_indices],[BL_weights;TL_weights;BR_weights;TR_weights],[Ny*Nx 1]);
%                 d_Image.CData([BL_indices;TL_indices;BR_indices;TR_indices])=[BL_weights;TL_weights;BR_weights;TR_weights];
                delete(h_circ);
                h_circ = viscircles([thisX,thisY],R,'Color','k','LineWidth',1,'EnhanceVisibility',false);
                pause(0.05);
            end
    end
    
    % Once there's a completed panel of time data, perform a derivative
    % operation in time.
    Qmat = sparse(rowIndex,colIndex,weights,Nt,Ny.*Nx);
    
    % This has a quantitative issue, but we ignore this for now. 
    Qdtemp = calculateDerivative(Qmat,dt); % Factor out and just loop over different panels of the original sparse matrix.

    Qdholder{xdcrIndex}=Qdtemp;
    
    
end

modelQ = vertcat(Qdholder{:});
% modelTime = tvec;



end

%% Helper functions

    function derivative = calculateDerivative(x,dx)
    derivative = sparse([],[],[],size(x,1),size(x,2),nnz(x)+(2*size(x,2)));
    
    derivative(2:(end-1),:) = (x(3:end,:)-x(1:(end-2),:))./(2.*dx);
%     derivative(2:(end-1),:) = ((x(3:end,:)-x(2:(end-1),:))./(dx) + (x(2:(end-1),:)-x(1:(end-2),:))./dx)/2;
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