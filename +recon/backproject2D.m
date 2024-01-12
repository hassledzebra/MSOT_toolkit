
function [bpModel,varargout] = backproject2D(input_coords,output_coords,varargin)
% backproject2D returns a lookup table of backprojection distances from
% each transducer to each pixel. This is used as the 'model' for the
% backproject2D solver, but does not actually fit the general description
% of what a model is. 
%
% It is possible to write the code such that the backprojection operator is
% written out as a matrix, but this causes a massive matrix operator to be
% created and kills the efficiency of the method. 




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
        Nxdcr = size(xdcrX,1);
        
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
        
        
    FovX = fieldOfView;
    FovY = fieldOfView;















%     varargout{1} = [];
% 
%     tvec = output_coords{1};
%         Nt = size(tvec,1);
%     xdcr = output_coords{2};
%         Nproj = size(xdcr,1);
%         xdcrX = xdcr(:,1);
%         xdcrY = xdcr(:,2);
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
%     speedOfSound = model_args.SpeedOfSound;
%         
%         
%     FovX = range(xpix)+dx;
%     FovY = range(ypix)+dy;
    
%     Nproj = reconstructor.N_Projections;
%     Nxdcr = reconstructor.N_Transducers;
%     Nx = reconstructor.N_x;
%     Ny = reconstructor.N_y;
%     Nt = reconstructor.N_Samples;
    
%     Fs = reconstructor.SamplingFrequency;
    
    % interpolate the x and y coordinates in the polar domain.
        % This should already be done in the refactored data.
%     newTh = interp1(linspace(1,Nxdcr,Nxdcr),reconstructor.TransducerCoordinatesCylindrical_Theta,linspace(1,Nxdcr,Nproj));
    
%     [xCo,yCo] = pol2cart(newTh, mean(reconstructor.TransducerCoordinatesCylindrical_R));
            
    xCo = xdcrX;
    yCo = xdcrY;
    
            % Create the grids for the X and Y pixel coordinates.
%             xVector = linspace( -FovX/2 , FovX/2 , Nx );
%             yVector = linspace( -FovY/2 , FovY/2 , Ny );
            xVector = xpix;
            yVector = ypix;
            
            [imX,imY] = meshgrid( xVector, yVector);
            % We should have the coordinates, but not the grids.
            
            % Array of sample timepoints.
%             bpTimes = double( ( 1 : Nt)' )./  Fs;
                % We should have this. 
            
            
            % Distances and times to each pixel.
            Xdist = xCo(:)' - imX(:);
            Ydist = yCo(:)' - imY(:);
            Rdist = hypot(Xdist,Ydist);
            bpModel = Rdist ./ speedOfSound;
%             bpProjList = 1:Nproj;


%         lowEdge = arrayfun(@(x) find(x>=tvec,1,'last'),bpModel(:));
%         highEdge = arrayfun(@(x) find(x<=tvec,1,'first'),bpModel(:));
% 
%         lowDiff = abs(tvec(lowEdge)-bpModel(:));
%         highDiff = abs(tvec(highEdge)-bpModel(:));
%         
% %         interpMat = sparse(10000,1033*50);
%         
%         pixCo = repmat((1:(Ny.*Nx))',[1 Nproj]);
%         transCo = repmat((1:Nproj),[Ny.*Nx,1]);
%         
%         highJ = sub2ind([Nt Nproj],highEdge,transCo(:));
%         lowJ = sub2ind([Nt Nproj],lowEdge,transCo(:));
%         
% %         highSub = sub2ind([Ny.*Nx,1033,50],pixCo(:),highEdge,transCo(:));
% %         lowSub = sub2ind([10000,1033,50],pixCo(:),lowEdge,transCo(:));
% %         interpMat(highSub) = lowDiff./(lowDiff+highDiff);
% %         interpMat(lowSub) = highDiff./(lowDiff+highDiff);
%         highVal = lowDiff./(lowDiff+highDiff);
%         lowVal = highDiff./(lowDiff+highDiff);
%         
%         interpMat = sparse([pixCo(:);pixCo(:)],[highJ(:);lowJ(:)],[highVal(:);lowVal(:)],Ny.*Nx,Nt*Nproj);
%         bpModel = interpMat;
end