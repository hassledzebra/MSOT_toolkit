
% function [Qd, modelTimeList,varargout] = CDMMI(reconObj,varargin)  
function [Qd, varargout] = CDMMI(input_coords,output_coords,varargin)
%               where input_coords represents the coordinates of the input space
%               (that which the model acts on), similarly for output_coords.
varargout{1} = [];
    % 
%     tvec = output_coords{1};
%         Nt = size(tvec,1);
%     xdcr = output_coords{2};
%         Nxdcr = size(xdcr,1);
%         sensorx = xdcr(:,1);
%         sensory = xdcr(:,2);
%         
%     ypix = input_coords{1};
%         Ny = numel(ypix);
%         dy = mean(diff(ypix));
%         nodeMinY = ypix(1);
%     xpix = input_coords{2};
%         Nx = numel(xpix);
%         dx = mean(diff(xpix));
%         nodeMinX = xpix(1);
% 
%         dt = mean(diff(tvec));
    
    
    
    %%%%%%%%%%%%%%%%%%%
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
        sensorx = xdcrvec.getCoordinateValues('EUCLIDEAN_X'); % xdcrvec.getCoordinateValues('EUCLIDEAN_X');
        sensory = xdcrvec.getCoordinateValues('EUCLIDEAN_Y');
        
        sensorx = sensorx(:);
        sensory = sensory(:);
        
        Nxdcr = size(sensorx,1);
        
    else
        tvec = output_coords{1};
        Nt = size(tvec,1);
        xdcr = output_coords{2};
        Nxdcr = size(xdcr,1);
        sensorx = xdcr(:,1);
        sensory = xdcr(:,2);
    end
    
    if isa(input_coords,'frame.CoordinateSystem')
        ypix = input_coords.getCoordinateValues('EUCLIDEAN_Y');
        Ny = numel(ypix);
        dy = mean(diff(ypix));
        nodeMinY = ypix(1);
        xpix = input_coords.getCoordinateValues('EUCLIDEAN_X');
        Nx = numel(xpix);
        dx = mean(diff(xpix));
        nodeMinX = xpix(1);
        
    else
        ypix = input_coords{1};
        Ny = numel(ypix);
        dy = mean(diff(ypix));
        nodeMinY = ypix(1);
        xpix = input_coords{2};
        Nx = numel(xpix);
        dx = mean(diff(xpix));
        nodeMinX = xpix(1);
    end
    dt = mean(diff(tvec));
    speedOfSound = model_args.SpeedOfSound;
    fieldOfView = range(xpix)+dx;
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
    
%     time_ = varargin{1};
%     xdcr_ = varargin{2};
%     xpix_ = varargin{3};
%     ypix_ = varargin{4};
%     
%     dy_ = mean(diff(ypix_));
%     dx_ = mean(diff(xpix_));
    
    % Dump step. To be factored out. 
%     if  isprop(reconObj,'dump') && reconObj.dump==true
%         s = fileread([mfilename('fullpath'),'.m']);
%         dump.codeString = s;
%         dump.inputStruct = struct(reconObj);
%         dump.inputJSON = jsonencode(struct(reconObj));
%         varargout{1} = dump;
%     else
%         varargout{1} = [];
%     end
    
    %
    % Input_Time_Coords = input_coords{1};
    % Input_Xdcr_Coords = input_coords{2};
    %   Xdcr_coords_X = Input_Xdcr_Coords(:,1);
    %   Xdcr_coords_Y = Input_Xdcr_Coords(:,2);
    
    % Output_Y_Coords = output_coords{1};
    % Output_X_Coords = output_coords{2};
    
    % SpeedOfSound = model_args.SpeedOfSound;
    
    
%     fieldOfView = reconObj.FieldOfView_X;
%     Nx = reconObj.N_x;
%     transducerCoordinatesR = reconObj.TransducerCoordinatesCylindrical_R;
%     transducerCoordinatesTheta = reconObj.TransducerCoordinatesCylindrical_Theta;
                                     
                     
                    
   
    
    %TODO: Post a warning if this isn't scalar (Math doesn't work for varying SoS)'
%     speedOfSound = reconObj.SpeedOfSound;
%     timeres = reconObj.PixelOversampling;
    
    
    
    % Pixel sizes.
%     px = fieldOfView/(n_x);
%     py = fieldOfView/(n_x);
%     
%     % Radius of the circumcircle of the square field of view. 
%     RFOV = fieldOfView*sqrt(2)/2;
%     Rclose = min(transducerCoordinatesR-RFOV); % Distance between xdcr and the
%                                                % nearest point of the
%                                                % circumcircle.
%     Rfar = max(transducerCoordinatesR+RFOV);   % Distance between xdcr and the
                                               % furthest point of the
                                               % circumcircle.
    
    % Timepoints at which the model is sampled. 
%     modelTimeList = linspace(Rclose/speedOfSound,Rfar/speedOfSound,Nx*timeres);
    modelTimeList = tvec;
    
    %% Define the grid giving us the image area.
    
    % Corners and centers for the marching tests. Allows us to determine what
    % kind of intersection we're encountering. Might be able to delete this. 
    marchX = [-fieldOfView/2 , 0 , fieldOfView/2];
    marchY = [-fieldOfView/2 , 0 , fieldOfView/2];
%     XYCoo = (n_x-1)/2*px;
    
%     [marchXgrid,marchYgrid]=meshgrid(marchX,marchY);
    
    % Actual grid (TODO)
%     z = zeros(n_x);
    
    
    %% Calculate distances to each marching test point for each sensor.
%     [sensorx,sensory] = pol2cart(transducerCoordinatesTheta,transducerCoordinatesR);

%     
%     
%     xdcrXStart = reconObj.TransducerCoordinatesCartesian_X;
%     xdcrYStart = reconObj.TransducerCoordinatesCartesian_Y;
%     
%     [theta,rho] = cart2pol(xdcrXStart,xdcrYStart);
%     [thetaNew] = interp1(1:numel(theta),unwrap(theta),linspace(1,numel(theta),N_proj));
%     [rhoNew] = interp1(1:numel(rho),rho,linspace(1,numel(theta),N_proj));
%     [sensorx,sensory] = pol2cart(thetaNew,rhoNew);
    
    % index 1 is the bottom left corner, index 3 is the upper left, index 9 is the upper right.
    
    % We now know how many points will have been passed at each time point, as
    % well as the first and last time points for each sensor.
    
    % Classify each sensor as either a corner or an edge, within the 3x3
    % neighborhood around the field of view. Transducers within corner regions
    % are guaranteed to intersect with the field of view first, and are
    % guaranteed to do so at the corner pixel. 
    % 
    %       C | E | C 
    %       - - - - -
    %       E | F | E
    %       - - - - -
    %       C | E | C
    
    onLeft = sensorx<=marchX(1);
    onRight = sensorx>=marchX(end);
    
    onBottom = sensory<=marchY(1);
    onTop = sensory>=marchY(end);
    
    edgeCornerParity = onLeft+onRight+onTop+onBottom;
    
    isEdge=(edgeCornerParity==1);
    isCorn=(edgeCornerParity==2);
    
    
    
    
    %% Make arrays of the distance of each corner sensor to its nearest edge.
    % Each row is the distance to the nearest edge. 1/3 are the close/far
    % vertical edges, 2/4 are the close/far horizontal edges.
    distToEdges = zeros(4,Nxdcr);
    
    % Corner Sensors
    
    %%
    distToEdges(1,:) = abs(sensorx-marchX(1));
    distToEdges(2,:) = abs(sensory-marchY(1));
    distToEdges(3,:) = abs(sensorx-marchX(end));
    distToEdges(4,:) = abs(sensory-marchY(end));
    
    %% Classify each sensor by which region it's in.
    %       1 | 2 | 3 
    %       - - - - -
    %       8 | F | 4
    %       - - - - -
    %       7 | 6 | 5
    
    sensorRegion=zeros(1,Nxdcr);
    sensorRegion(onTop&onLeft)=1;
    sensorRegion(onTop&isEdge)=2;
    sensorRegion(onTop&onRight)=3;
    sensorRegion(onRight&isEdge)=4;
    sensorRegion(onBottom&onRight)=5;
    sensorRegion(onBottom&isEdge)=6;
    sensorRegion(onBottom&onLeft)=7;
    sensorRegion(onLeft&isEdge)=8;
    
    
    % Corrections to angles calculated using inverse trig functions. These are
    % magic numbers, but help a lot with the efficiency of the calculation. 
    
    corrAng=pi*[1.5 1.0 1.5 1.0 0.5 1.0 0.5 0.0
                2.0 1.5 1.0 1.5 1.0 0.5 0.0 1.5
                1.5 2.0 1.5 1.0 0.5 0.0 0.5 0.0
                2.0 1.5 1.0 0.5 1.0 0.5 0.0 0.5];
    
    corrAngPlus=logical([1 1 0 1 1 0 0 1
                         0 1 1 0 0 1 1 1
                         1 0 0 1 1 1 0 1
                         0 1 1 1 0 1 1 0]);
    
    corrAngMinus=logical([0 0 1 1 0 1 1 1
                          1 1 0 1 1 1 0 0
                          0 1 1 1 0 0 1 1
                          1 1 0 0 1 1 0 1]);
    
    
    
%     thorta = 0:pi/150:2*pi;
    roInd=[];
    coInd=[];
    valInd=[];
    %% The ve;tor pointing to the closest point from each sensor denotes the point which is 'farthest into the image'.
%     Nt = numel(modelTimeList)
%     DEBUG=true;
     DEBUG=false;
    if DEBUG
        %%
        d_Image_Data = zeros(Ny,Nx);
        d_Fig = figure;
        d_Ax = axes;
        d_Image = imagesc(d_Ax,d_Image_Data,'XData',[min(xpix) max(xpix)],'YData',[min(ypix) max(ypix)]); hold on;
        d_Ax.YDir='normal';
        xdcrScat = scatter(d_Ax,sensorx(:),sensory(:));
        
        q_mid= quiver(d_Ax,sensorx(1),sensory(1),-0.01,-0.01,'AutoScale','off');
        q_int= quiver(d_Ax,sensorx(1),sensory(1),-0.01,-0.01,'AutoScale','off');
        
        scat_mid=scatter(0,0,'r');
        scat_int=scatter(0,0,'m');
        
        axis([min(sensorx) max(sensorx) min(sensory) max(sensory)]);
        pbaspect([1 1 1]);
        
        miX = min(marchX);
        maX = max(marchX);
        miY = min(marchY);
        maY = max(marchY);
        
        plot([miX,miX],[miY,maY],'k');
        plot([maX,maX],[miY,maY],'k');
        plot([miX,maX],[maY,maY],'k');
        plot([miX,maX],[miY,miY],'k');
        
%         h_circ = viscircles([sensorx(1),sensory(1)],0.02);
        
        for d_k = 1:(Ny-1)
            plot([miX,miX]+dx.*d_k,[miY,maY],'k');
            plot([miX,maX],[miY,miY]+dy.*d_k,'k');
        end
        
        
    end
    
    firstFlag=0;
    
    for xdcrInd = 1:Nxdcr
        for tInd = 1:Nt
            
            % Get the transducer's distance to each of the edges of the field of
            % view. 
            xdcrDist = distToEdges(:,xdcrInd);
            
            % Fetch the corrections for that transducer region.
            xdcrType = sensorRegion(xdcrInd); 
            angleCorrection = corrAng(:,xdcrType);
            anglePlus = corrAngPlus(:,xdcrType);
            angleMinus = corrAngMinus(:,xdcrType);
            
            % Sensor coordinates. 
            xdcrLocX = sensorx(xdcrInd);
            xdcrLocY = sensory(xdcrInd);
            
            % Radius of the circle defining the wavefront at this timepoint. 
            R = speedOfSound*modelTimeList(tInd);
            
            % Calculate angles based on edge/corner.
            if isCorn(xdcrInd)
                theseAngles = asin(xdcrDist./R);
            elseif isEdge(xdcrInd)
                theseAngles = acos(xdcrDist./R);
            else
                error('Something went wrong when classifying transducer regions');
            end
            
            % Correct angles so that they're in the right orientation. 
            thetaTemp = [angleCorrection(anglePlus) + theseAngles(anglePlus); 
                         angleCorrection(angleMinus)-theseAngles(angleMinus)];
            
            % We only want real angles, and we want them ascending. 
            theta = sort(thetaTemp(thetaTemp==real(thetaTemp)));
            
            % Remove the angles that match the trig geometry but do not lie
            % inside of the field of view. 
            [xadj,yadj] = pol2cart(theta,R);
            pointPosX = xadj + xdcrLocX;
            pointPosY = yadj + xdcrLocY;
            badAngles = (pointPosX < (marchX(1)   - 1E-16) | ...
                         pointPosX > (marchX(end) + 1E-16)) | ...
                        (pointPosY < (marchY(1)   - 1E-16) | ...
                         pointPosY > (marchY(end) + 1E-16));
            thetaChoice = theta(~badAngles); 
            
            % At this point we should have only valid angles which intersect the
            % geometry of the field of view. This can have either 0, 1, 2, 3, 
            % or 4 points. The zero-point case indicates that this transducer 
            % has not yet reached the FOV, so we continue to the next timept.
            % The single-point case signifies that the wave front is
            % exactly tangent to the viewing region, the two-point case
            % indicates that the wave front is one continuous arc inside the
            % FOV, the three-point case indicates that we have an arc inside the
            % FOV which is just about to leave and is exactly tangent to the
            % interior of the FOV, and the four-point case indicates that we
            % have a wave front which is leaving the FOV and has split into
            % three pieces, with the center piece outside of the FOV but the
            % other two inside. 
            
            if isempty(thetaChoice)
                continue;
            end
%             thetaChoice=sort(thetaChoice);
            
            %Correct for wraparound. Since all transducers are outside of the
            %FOV, the angle between the two intersection angles is guaranteed to
            % be less than pi/2. 
            
            if numel(thetaChoice) == 2
                if abs(diff(abs(thetaChoice)))>(pi/2)
                    thetaChoice(2)=thetaChoice(2)-2*pi;
                end
            elseif numel(thetaChoice)==4
                
            else
                
                error('Something terrible has happened');
            end
            
            %% We now have intersections of each vector with the edge of the ROI.
            %There will only be two or four intersections; split up by case
            %to calculate the arc length inside the ROI.
            if numel(thetaChoice)==2
                arcLength=(thetaChoice(2)-thetaChoice(1));
            elseif numel(thetaChoice)==4
                arcLength=((thetaChoice(2)-thetaChoice(1))+(thetaChoice(4)-thetaChoice(3)));
            else
                error('Something terrible has happened');
            end
            
            % Calculate the intersection points with each gridline.
            th1x=acos(((linspace(-fieldOfView/2,fieldOfView/2,Nx+1)-xdcrLocX))/R);
            th1y=asin(((linspace(-fieldOfView/2,fieldOfView/2,Nx+1)-xdcrLocY))/R);
            
            % This is all the result of a lot of geometry diagrams and debugging
            % magic. Basically just corrects the angles so that they have the
            % same, sensible ordering for both x and y, as well as for each
            % transducer region. 
            switch xdcrType
                case 1
                    th1x=2*pi-th1x;
                    th1y=th1y+2*pi;
                case 2
                    th1x=(pi/2)-th1x+3*pi/2;
                    th1y=2*pi+th1y;
                    th1y=[3*pi/2-(th1y-3*pi/2) th1y];
                case 3
                    th1x=2*(pi)-th1x;
                    th1y=pi-th1y;
                case 4
                    th1x=[2*pi-th1x th1x];
                    th1y=pi-th1y;
                case 5
                    th1y=pi-th1y;
                case 6
                    th1y=[th1y pi-th1y];
                case 7
                case 8
                    thetaChoice(thetaChoice>(pi/2))=thetaChoice(thetaChoice>(pi/2))-2*pi;
                    thetaChoice=sort(thetaChoice);
                    th1x=[-th1x th1x];
            end
            
            % bundle and retain only those arc intersections which lie within
            % the angle range we determined earlier. 
            
            th = [th1x(th1x==real(th1x)) th1y(th1y==real(th1y))];
            th = th(th>min(thetaChoice)&th<max(thetaChoice));
            
            % Ignore any intersections which lay outside the FOV.
            if numel(thetaChoice)==4
                th1=sort([thetaChoice(1:2)' th(th<thetaChoice(2))]);
                th2=sort([thetaChoice(3:4)' th(th>thetaChoice(3))]);
                
                mdth1 = (th1(2:end)+th1(1:end-1))/2;
                mdth2 = (th2(2:end)+th2(1:end-1))/2;
                
                delth1 = (th1(2:end)-th1(1:end-1));
                delth2 = (th2(2:end)-th2(1:end-1));
                
                mdth = [mdth1(:);mdth2(:)];
                delth = [delth1(:);delth2(:)];
                
            else
                th=sort([thetaChoice' th]);

%                 Get the midpoint of each arc, and the size of that arc. 
                mdth = (th(2:end)+th(1:end-1))/2;
                delth = (th(2:end)-th(1:end-1));
                
            end
            
            % Get the X, Y coordinates of the midpoints, in terms of pixel
            % indices.
            xCrossInd = (xdcrLocX+R*cos( mdth )-marchX(1))/dx;
            yCrossInd = (xdcrLocY+R*sin( mdth )-marchY(1))/dy;
            
            
            %% We now have each pixel's value and coordinate. We can now perform what amounts to a reverse interpolation.
            
            % We can't have indices outside the range of (0,N_x) or (0,N_y)
%             goodMask=(xCrossInd<(Nx+1))&(yCrossInd<(Ny+1))&(xCrossInd>-1)&(yCrossInd>-1);
            goodMask=(xCrossInd<(Nx))&(yCrossInd<(Ny))&(xCrossInd>0)&(yCrossInd>0);
            xCrossInd=(xCrossInd(goodMask));
            yCrossInd=(yCrossInd(goodMask));
            delth=abs(delth(goodMask)); % Can't have a negative arclength.
            
            
            
%             [crossInd,vals] = reverseInterpolation(xCrossInd,yCrossInd,delth,Ny,Nx,ypix,xpix);
            
            
            
            
            
            if DEBUG %&& tInd>1380
                %%
                q_mid.Visible = 'off';
                q_mid.UData = R*cos( mdth );
                q_mid.VData = R*sin( mdth );
                q_mid.XData = ones(size(q_mid.UData)).*sensorx(xdcrInd);
                q_mid.YData = ones(size(q_mid.VData)).*sensory(xdcrInd);
                
                q_int.Visible = 'off';
                q_int.UData = R*cos( th );
                q_int.VData = R*sin( th );
                q_int.XData = ones(size(q_int.UData)).*sensorx(xdcrInd);
                q_int.YData = ones(size(q_int.VData)).*sensory(xdcrInd);
                
                scat_mid.XData = q_mid.XData+q_mid.UData;
                scat_mid.YData = q_mid.YData+q_mid.VData;
                scat_int.XData = q_int.XData+q_int.UData;
                scat_int.YData = q_int.YData+q_int.VData;
%                 scat_mid.CData = delth;
            end
            
            [crossInd,vals] = reverseInterpolation(xCrossInd,yCrossInd,delth,Ny,Nx);
            
            if DEBUG
                
                d_Image.CData(:)=0;
%                 d_Image.CData(crossInd)=vals;
                d_Image.CData(:) = accumarray(crossInd(:),vals(:),[Ny*Nx 1]);
%                 delete(h_circ);
%                 h_circ = viscircles([sensorx(xdcrInd),sensory(xdcrInd)],R,'Color','k','LineWidth',1,'EnhanceVisibility',false);
%                 th_circ = linspace(0,2*pi,2000);
%                 [xx,yy]=pol2cart(th_circ,R);
%                 h_circ = plot(xx+sensorx(xdcrInd),yy+sensory(xdcrInd),'Color','k');
                
                pause(0.05);
            end
            

%             z(:) = 0;
%             z(crossInd) = vals;
            
            %store in a matrix.
            numVals=numel(crossInd);
            roInd(end+(1:numVals))=tInd;
            coInd(end+(1:numVals))=crossInd;
            valInd(end+(1:numVals))=vals * (1./(4*pi*speedOfSound));
            
            
            
            
        end
        %% Calculate the time derivatives based on the particular transducer's matrix.
        Qmat=sparse(roInd,coInd,valInd,Nt,Nx*Nx);
        Qdtemp = calculateDerivative(Qmat,dt);
        Qdholder{xdcrInd}=Qdtemp;
        
        roInd=[];
        coInd=[];
        valInd=[];
        
        
        
        
    end
    
    Qd=vertcat(Qdholder{:});
    
    
    
end

function derivative = calculateDerivative(x,dx)
    derivative = sparse([],[],[],size(x,1),size(x,2),nnz(x)+(2*size(x,2)));
    
    derivative(2:(end-1),:) = (x(3:(end),:)-x(1:(end-2),:))./(2.*dx);
%      derivative1 = (x(2:(end-1),:)-x(1:(end-2),:))./(dx);
%      derivative2 = (x(3:(end),:)-x(2:(end-1),:))./(dx);
%      derivative(2:(end-1),:) = (derivative1+derivative2)./2;
    derivative(1,:) = (sum([-1.5; 2.0; -0.5].*x(1:3,:)))./dx;
    derivative(end,:) = sum([0.5; -2.0; 1.5].*x((end-2):end,:))./dx;
end


% THIS IS EXTREMELY INEFFICIENT IN TERMS OF MEMORY USAGE.
function [indout,valout] =reverseInterpolation(xindex,yindex,values,Nx,Ny,ypix,xpix)
    %Returns a reverse-interpolated/ gridded set of indices which
    %correspond to the integer grid which would interpolate to the given values.
    
    %inputs:
    %       xindex, a list of indices (not necessarily integers) on an
    %       image grid, x values.
    % ditto y
    % values are the values at those x,y coordinates.
    % Nx is the number of pixels in the X direction.
    % Ny is the number of pixels in the Y direction.
    
    
    xindex=xindex(:);
    yindex=yindex(:);
    values=values(:);
    
%     xArcPt_hat = xindex+0.5; 
%         yArcPt_hat = yindex+0.5;
        
    xArcPt_hat = [xindex+0.5];%; xindex+0.6; xindex+0.5]; 
        yArcPt_hat = [yindex+0.5];% xindex+0.6; xindex+0.5]; 
        values = [values];%;values./4;values./2];
        
        
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
        
%         BL_weights = exp(-1./(1-BL_weights));
%         TL_weights = exp(-1./(1-TL_weights));
%         BR_weights = exp(-1./(1-BR_weights));
%         TR_weights = exp(-1./(1-TR_weights));
%         
        
        % Trim out bad indices and bad weights, convert to linear index
        [BL_indices , BL_mask] = convertIndices( BL_subs ,Nx,Ny);
        [TL_indices , TL_mask] = convertIndices( TL_subs ,Nx,Ny);
        [BR_indices , BR_mask] = convertIndices( BR_subs ,Nx,Ny);
        [TR_indices , TR_mask] = convertIndices( TR_subs ,Nx,Ny);
        
        % Remove bad pixels.
        BL_weights = values(BL_mask).*BL_weights(BL_mask);
        TL_weights = values(TL_mask).*TL_weights(TL_mask);
        BR_weights = values(BR_mask).*BR_weights(BR_mask);
        TR_weights = values(TR_mask).*TR_weights(TR_mask);
        
        indout = [BL_indices;TL_indices;BR_indices;TR_indices];
        valout = [BL_weights;TL_weights;BR_weights;TR_weights];
% %     
% %     dy = mean(diff(ypix));
% %     dx = mean(diff(xpix));
% %     Y0 = min(ypix-dy);
% %     X0 = min(xpix-dx);
% %     Ym = max(ypix+dy);
% %     Xm = max(xpix+dx);
% %     
% %     leftFrameX = -1.*ones(Ny,1);
% %     rightFrameX = (Nx+1).*ones(Ny,1);
% %     topFrameY = (Ny+1).*ones(Nx,1);
% %     botFrameY = -1.*ones(Nx,1);
% %     
% %     xindex = [xindex;leftFrameX;rightFrameX;[1:Ny].';[1:Ny].'];
% %     yindex = [yindex;[1:Nx].';[1:Nx].';topFrameY;botFrameY];
% %     values = [values;zeros(2.*Nx+2.*Ny,1)];
% %     
% %     [yy,xx] = ndgrid(1:Ny,1:Nx);
% %     
% %     vals = griddata(xindex,yindex,values,xx,yy,'v4');
% %     
% %     indout = find(vals);
% %     valout = vals(indout);
%     
%     
%     
%     
%     
% %     indout = sub2ind([Ny Nx],ceil(yindex),ceil(xindex));
% %     valout = values;
% %     return;
%     
%     %Because we're interpolating to the center of pixels in a grid, the indices
%     %and coordinates are offset by 0.5
%     
%     
%     
% %     xpix2 = 0.5:1:(Nx+0.5);
% %     ypix2 = 0.5:1:(Ny+0.5);
% %     
% %     [xx,yy] = meshgrid(xpix2,ypix2);
% %     [Idx,D] = rangesearch([xx(:) yy(:)],[xindex yindex],1);
% %     scaledD = cellfun(@(x,y) exp(-1./(1-x.^2)),D,mat2cell(values,ones(numel(values),1)),'uni',false);
% %     indout = cat(2,Idx{:})';
% %     valout = cat(2,scaledD{:})';
%     
%     westIndex = floor(xindex-0.5)+0.5;
%     eastIndex = floor(xindex+0.5)+0.5;
%     southIndex = floor(yindex-0.5)+0.5;
%     northIndex = floor(yindex+0.5)+0.5;
%     
%     %% Cleanup
%     % Because accumarray won't take 0s, we have to remove those (this
%     % causes signal leakage out the side of the image, but that's fine).
%     % Also remove indices which are Too High, i.e. larger than Nx or Ny.
% 
%     
%     badMaskNorth = northIndex>(Ny);
%     badMaskSouth = southIndex<0;
%     badMaskEast = eastIndex>(Nx);
%     badMaskWest = westIndex<0;
%     
%     northwestMask = (~badMaskNorth)&(~badMaskWest);
%     northeastMask = (~badMaskNorth)&(~badMaskEast);
%     southwestMask = (~badMaskSouth)&(~badMaskWest);
%     southeastMask = (~badMaskSouth)&(~badMaskEast);
%     
%     
%     % Get the distances to each of the pixel grid points.
%     westDist=xindex-westIndex;
%     eastDist=eastIndex-xindex;
%     southDist=yindex-southIndex;
%     northDist=northIndex-yindex;
%     
%     
%     %% Use these to make the interpolated grid points.
%     
%     %Coordinates
%     %Y coordinate is in the first column, but don't make coordinates that
%     %are out of bounds.
%     NorthwestPointCoordinates=ceil([northIndex(northwestMask),westIndex(northwestMask)]);
%     NortheastPointCoordinates=ceil([northIndex(northeastMask),eastIndex(northeastMask)]);
%     SouthwestPointCoordinates=ceil([southIndex(southwestMask),westIndex(southwestMask)]);
%     SoutheastPointCoordinates=ceil([southIndex(southeastMask),eastIndex(southeastMask)]);
%     
%     %Values in each of those points.
%     NorthwestPointValue=values(northwestMask).*(southDist(northwestMask).*eastDist(northwestMask));
%     NortheastPointValue=values(northeastMask).*(southDist(northeastMask).*westDist(northeastMask));
%     SouthwestPointValue=values(southwestMask).*(northDist(southwestMask).*eastDist(southwestMask));
%     SoutheastPointValue=values(southeastMask).*(northDist(southeastMask).*westDist(southeastMask));
% 
%     
%     
%     
%     
%     %Accumulate those values into an image so there's no redundancy. Use
%     %sparsity.
%     %TODO: MAKE BETTER. VERY TIME AND SPACE CONSUMING
%     coordGroup = [NorthwestPointCoordinates ; NortheastPointCoordinates ; SouthwestPointCoordinates ; SoutheastPointCoordinates];
% %     coordMask = coordGroup(:,1)>0 & coordGroup(:,1)<(Nx+1) & coordGroup(:,2)>0 & coordGroup(:,2)<(Ny+1);
% %     coordGroup = coordGroup(coordMask,:);
%     valGroup = [NorthwestPointValue ; NortheastPointValue ; SouthwestPointValue ; SoutheastPointValue];
% %     valGroup = valGroup(coordMask);
%     indGroup = sub2ind([Ny,Nx],coordGroup(:,1),coordGroup(:,2));
%     
%     
%     sss= sparse(indGroup,1,valGroup);
%     sss = accumarray(coordGroup,valGroup,[Ny,Nx],[],[],true);
%     
%     
%     % Locate the indices of all values in the matrix.
%     indout = find(sss);
%     valout = nonzeros(sss);
% %     
%     
%     
end


% THIS IS EXTREMELY INEFFICIENT IN TERMS OF MEMORY USAGE.
% function [indout,valout] =reverseInterpolation(xindex,yindex,values,Nx,Ny)
%     %Returns a reverse-interpolated/ gridded set of indices which
%     %correspond to the integer grid which would interpolate to the given values.
%     
%     %inputs:
%     %       xindex, a list of indices (not necessarily integers) on an
%     %       image grid, x values.
%     % ditto y
%     % values are the values at those x,y coordinates.
%     % Nx is the number of pixels in the X direction.
%     % Ny is the number of pixels in the Y direction.
%     
%     
%     xindex=xindex(:);
%     yindex=yindex(:);
%     values=values(:);
%     
%     
%     %Because we're interpolating to the center of pixels in a grid, the indices
%     %and coordinates are offset by 0.5
%     westIndex = floor(xindex-0.5)+0.5;
%     eastIndex = floor(xindex+0.5)+0.5;
%     southIndex = floor(yindex-0.5)+0.5;
%     northIndex = floor(yindex+0.5)+0.5;
%     
%     %% Cleanup
%     % Because accumarray won't take 0s, we have to remove those (this
%     % causes signal leakage out the side of the image, but that's fine).
%     % Also remove indices which are Too High, i.e. larger than Nx or Ny.
% 
%     
%     badMaskNorth = northIndex>(Ny);
%     badMaskSouth = southIndex<0;
%     badMaskEast = eastIndex>(Nx);
%     badMaskWest = westIndex<0;
%     
%     northwestMask = (~badMaskNorth)&(~badMaskWest);
%     northeastMask = (~badMaskNorth)&(~badMaskEast);
%     southwestMask = (~badMaskSouth)&(~badMaskWest);
%     southeastMask = (~badMaskSouth)&(~badMaskEast);
%     
%     
%     % Get the distances to each of the pixel grid points.
%     westDist=xindex-westIndex;
%     eastDist=eastIndex-xindex;
%     southDist=yindex-southIndex;
%     northDist=northIndex-yindex;
%     
%     
%     %% Use these to make the interpolated grid points.
%     
%     %Coordinates
%     %Y coordinate is in the first column, but don't make coordinates that
%     %are out of bounds.
%     NorthwestPointCoordinates=ceil([northIndex(northwestMask),westIndex(northwestMask)]);
%     NortheastPointCoordinates=ceil([northIndex(northeastMask),eastIndex(northeastMask)]);
%     SouthwestPointCoordinates=ceil([southIndex(southwestMask),westIndex(southwestMask)]);
%     SoutheastPointCoordinates=ceil([southIndex(southeastMask),eastIndex(southeastMask)]);
%     
%     %Values in each of those points.
%     NorthwestPointValue=values(northwestMask).*(southDist(northwestMask).*eastDist(northwestMask));
%     NortheastPointValue=values(northeastMask).*(southDist(northeastMask).*westDist(northeastMask));
%     SouthwestPointValue=values(southwestMask).*(northDist(southwestMask).*eastDist(southwestMask));
%     SoutheastPointValue=values(southeastMask).*(northDist(southeastMask).*westDist(southeastMask));
%     
%     
%     
%     
%     %Accumulate those values into an image so there's no redundancy. Use
%     %sparsity.
%     %TODO: MAKE BETTER. VERY TIME AND SPACE CONSUMING
%     coordGroup = [NorthwestPointCoordinates ; NortheastPointCoordinates ; SouthwestPointCoordinates ; SoutheastPointCoordinates];
%     coordMask = coordGroup(:,1)>0 & coordGroup(:,1)<(Nx+1) & coordGroup(:,2)>0 & coordGroup(:,2)<(Ny+1);
%     coordGroup = coordGroup(coordMask,:);
%     valGroup = [NorthwestPointValue ; NortheastPointValue ; SouthwestPointValue ; SoutheastPointValue];
%     valGroup = valGroup(coordMask);
%     indGroup = sub2ind([Ny,Nx],coordGroup(:,1),coordGroup(:,2));
%     
%     
%     sss= sparse(indGroup,1,valGroup);
%     
%     
%     % Locate the indices of all values in the matrix.
%     indout = find(sss);
%     valout = nonzeros(sss);
%     
%     
%     
% end
% 
% 
 % Converts a given set of coordinate subscripts into their
    % corresponding indices, removing the bad ones in the process.
    function [subIndices,mask]=convertIndices(subs,Nx_real,Ny_real)
        % If the subs are too low or too high, remove.
        mask = ( subs(:,1) >= 1 & subs(:,1) <= (Nx_real) ) & ...
               ( subs(:,2) >= 1 & subs(:,2) <= (Ny_real) );
        subIndices = ( subs(mask,2) + ( subs(mask,1) - 1 ) * Ny_real );
        
    end