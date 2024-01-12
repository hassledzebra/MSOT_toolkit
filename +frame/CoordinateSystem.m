classdef CoordinateSystem
    %CoordinateSystem defines intra-data relationships and interpretations
    %according to a linear algebraic construction of the underlying space.
    %
    
    properties
        Coordinates % The list of Coordinate objects which compose the CoordinateSystem's axes.
        
        SpaceConstruction % String which describes the topology of the particular coordinate system.
                          
                          % Reserved spaces are: 
                          % []  : Denotes the empty set
                          % [W] : Denotes the universal (omega) set.
                          % [U] : Denotes the limited set [-1,0,1]
                          % [N] : Denotes all natural numbers (1,2,3,...)
                          % [Z] : Denotes all integers (...,-2,-1,0,1,2,...)
                          % [R] : Denotes all real numbers 
                          % [Q] : Denotes all rational numbers 
                          % [C] : Denotes all complex numbers 
                          % [H] : Denotes all quaternionic numbers
                          % [O] : Denotes all octionionic numbers
                          % ['X'] : Denotes any user-defined space, where X
                          % corresponds to the Name variable of the
                          % coordinate.
                          % Note: None of these spaces as-is include values
                          % of infinity.
                          
                          % Reserved operators:
                          % (U) : Union of two spaces.
                          % (/) : Quotient of two spaces.
                          % (\) : Set complement; for B(\)A, the set of
                          %             items in B but not in A.
                          % (+) : Direct sum of two coordinates. 
                          % (x) : Tensor product of two coordinates.

                          % (k) : Kronecker product of two coordinates
                          % TODO: Probably want to include the ability to
                          % conjugate a coordinate, so that e.g. a
                          % vector/vector tensor product correctly gives us
                          % a (2,0)-tensor, similarly for a
                          % vector/covector.
                          
                          
                          
                          
                          % Reserved values are:
                          % 0  : Zero.
                          % 1  : One.
                          % PI : Pi.
                          % E  : Euler constant, ~2.718...
                          % INF : Infinity
                          % EPS : Dual unit, where EPS^2 = 0
                          % I  : Complex unit, where I^2 = -1
                          % J  : Split-Complex unit, where J^2 = +1
                          
                          % Values may be modified by prefixing them with +
                          % or - to denote the positive or negative set, or
                          % ^(number) to denote a space of a given
                          % dimension. The exponent may be represented as a
                          % product itself, so e.g. R^(3 x 2) is valid,
                          % representing the space of 3 x 2 matrices. 
                      
                          
                          % Examples:
                          %     - "R(/){2*PI*Z}" denotes the circle group,
                          %     equivalent to all angles measured in
                          %     radians.
                          %     - "R(/){360*Z}" denotes the same circle
                          %     group, but measured in degrees instead. 
                          %     - "R(\)Q" denotes all irrational numbers.
                          %     - "R(+)R" denotes all pairs of real
                          %     numbers. Entries are [2 x 1] values.
                          %     Equivalent to "R^2"
                          %     - "R(x)R" denotes a mapping between all
                          %     real numbers and all real numbers; this is
                          %     equivalent to standard multiplication.
                          %     - "R^2(x)R^2" denotes a mapping between all
                          %     two-entry vectors, equivalent to a 2x2
                          %     matrix.
                          %     - "+R" denotes the set of all positive
                          %     reals.
                          %     - "R(\)-R" denotes the set of all
                          %     non-negative reals.
                          %
                          %
                          % Say a user has defined four coordinates Y,X,R,
                          % and S, where Y has 16 1x1 entries, X has 16
                          % 1x1 entries, R has 32 2x1 entries, and S has 100
                          % 1x1 entries.
                          % 
                          % They may first define a CoordinateSystem I with
                          % "['Y'](x)['X']", and a second CoordinateSystem
                          % D with "['S'](x)['R']". I might stand for the
                          % space of images with coordinates Y,X, and
                          % similarly for D.
                          %
                          % They may then define a CoordinateSystem M with
                          % "['D'](x)['I']" to describe the mapping from I
                          % to D
           
                          
    end
    properties(Dependent)
        Values
        Names
    end
    
    methods
        function Values = get.Values(obj)
           
            
            Values = cellfun(@(x) x.Values,obj.Coordinates,'uni',false);
        end
        
        
        
        function obj = CoordinateSystem(varargin)
            
            nvargs = numel(varargin);
            
            
            % If we're just given coordinates, just join them as a tensor
            % product.
            isCoordinate = cellfun(@(x) isa(x,'frame.Coordinate'),varargin);
            isCoordinateSystem = cellfun(@(x) isa(x,'frame.CoordinateSystem'),varargin);
            
            if all(isCoordinate) && ~isempty(varargin)
                obj.Coordinates = [varargin(:)];
                
                coordinateBind = ["['" , repmat("'(x)'",[1,nvargs-1]) , "']"];
                coordinateNames = cellfun(@(x) x.Alias,varargin,'uni',false);
                obj.SpaceConstruction = strjoin(coordinateBind,coordinateNames);
                
            elseif all(isCoordinate|isCoordinateSystem) && ~isempty(varargin)
                coordinateBind = ["" , repmat("(x)",[1,nvargs-1]) , ""];
                for k = 1:nvargs
                    if isCoordinate(k)
                        coordinateNames(k) = strcat("['",string(varargin{k}.Alias),"']");
                    else
                            
                        coordinateNames(k) = strcat("[",strip(varargin{k}.SpaceConstruction,'"'),"]");
                    end
                    obj.Coordinates{k} = varargin{k};
                end
                obj.SpaceConstruction = strjoin(coordinateBind,coordinateNames);
                
            elseif nvargs>1 % Assume we were given coordinates and a space construction.
                lastCoordinate = find(isCoordinate,1,'last');
                if isstring(varargin{lastCoordinate+1})
                   % TODO: Add validation for this.
                   obj.SpaceConstruction = varargin{lastCoordinate+1};
                end
                
                obj.Coordinates = [varargin{1:lastCoordinate}];
            else
                
            end
            
            % Make sure coordinates are sorted according to their use in
            % the space construction. 
%             fixCoordinateOrder
%           % Collapse any direct sums into a single coordinate.
            x=0;
        end
        
        % Make a blank DataFrame. The default is to make a tensor array.
        function blankFrame = makeBlankFrame(obj,varargin)
            nVargs = numel(varargin);
            doSparse = false;
            % Input checking.
            if nVargs>1
                error('Only one argument permitted to makeBlankFrame');
            elseif nVargs == 1
                if ischar(varargin{1}) || isstring(varargin{1}) 
                   doSparse = strcmp(varargin{1},'sparse');
                elseif islogical(varargin{1})
                    doSparse = varargin{1};
                else
                    error('makeBlankFrame only accepts ''sparse'',"sparse",or logical as additional args');
                end
            end
            
            if doSparse && numel(obj.Coordinates)>2
               error("Sparse frames only supported for 2-D arrays.")
            end
            
            % Assign coordinate system
            blankFrame.CoordinateSystem = obj;
            
            if doSparse
                nFirst = numel(obj.Coordinates(1));
                try
                    nSecond = numel(obj.Coordinates(2));
                catch
                    nSecond = 1;
                end
                blankFrame.Data = sparse([],[],[],numel(nFirst,nSecond),sqrt(nFirst.*nSecond));
            else
                nCoords = arrayfun(@(x) numel(x),obj.Coordinates);
                blankFrame.Data = zeros([nCoords,1]);
            end
            
            
        end
        
        function vectorizeCoordinates = vectorizeCoordinates(obj)
            %vectorizeCoordinates Returns a vector array of the ordered
            %coordinates for each of the entries in a given object.
            %
            % Ex: A "['S'](x)['R']" coordinate system where S has 100 1x1 time
            % entries and R has 32 2x1 transducer coordinate entries would
            % vectorize as a 3200 x 3 coordinate, with entries 
            %  [S_1,Rx_1,Ry_1]
            %  [S_2,Rx_1,Ry_1]
            %  ...
            %  [S_100,Rx_1,Ry_1]
            %  [S_1,Rx_2,Ry_2]
            %  ...
            %  [S_100,Rx_32,Ry_32]
            %
            coordVals = {obj.Coordinates.Values};
            coordSizes = [1 cellfun(@(x) size(x,1),coordVals) 1];
            for k = 2:(numel(coordVals)+1)
                numToRight(1,k-1) = prod(coordSizes((k+1):end));
                numToLeft(1,k-1) = prod(coordSizes(1:(k-1)));
            end
            o = ones(1,numel(coordVals));
            
            elRep = cellfun(@(x,y) repelem(x,y,1),coordVals,mat2cell(numToLeft,1,o),'UniformOutput',false);
            matRep = cellfun(@(x,y) repmat(x,y,1),elRep,mat2cell(numToRight,1,o),'UniformOutput',false);
            vectorizeCoordinates = cat(2,matRep{:});
            x = 0;
        end
        
        function newObj = interpCoords(obj,N)
            newObj = obj;
            for k=1:size(obj.Coordinates,2)
               newObj.Coordinates(k) = interpCoords(obj.Coordinates(k),N); 
            end
        end
        
        
        % Get the values for each coordinate based on the standard name
        % given. If the named coordinate isn't compound, then
        % directly fetch the values, otherwise return the corresponding
        % Coordinate
        %
        function coordVals = getCoordinateValues(obj,coordinateStandardName)
           coordVals = {};
           coordName = cellfun(@(x) x.Names,obj.Coordinates,'uni',false);
           coName = strcmpi(coordName,coordinateStandardName);
           if any(coName) && sum(coName)==1
                coordVals = obj.Coordinates{coName};
                if (iscell(coordVals.Space)&&numel(coordVals.Space)==1) || ~iscell(coordVals.Space)
                    coordVals = coordVals.getChildValues;
                end
           end
           subCoordSystems = cellfun(@(x) isa(x,'frame.CoordinateSystem'),obj.Coordinates);
           if isempty(coordVals) && any(subCoordSystems)
               subCoordInds = find(subCoordSystems);
               for k = 1:numel(subCoordInds)
                   coordVals = obj.Coordinates{subCoordInds(1)}.getCoordinateValues(coordinateStandardName);
               end
           end
        end
        
        function tensSize = getTensorSize(obj) % Todo: Fix for categorical/endmember-type dimensions. 
            tensSize = [];
           for k = 1:numel(obj.Coordinates)
                if iscell(obj.Values{k})
                    tensSize = [tensSize;cellfun(@numel,obj.Values{k})];
                elseif isempty(obj.Values{k})
                    tensSize = [tensSize;numel(obj.Coordinates{k}.Values)];
                else
                    tensSize = [tensSize;numel(obj.Values{k})];
                end
           end
        end
        
        function Names = get.Names(obj)
           Names = cellfun(@(x) x.Names,obj.Coordinates,'uni',false); 
           for k = 1:numel(obj.Coordinates)
               if isa(obj.Coordinates{k},'frame.CoordinateSystem')
                   Names{k} = obj.Coordinates{k}.Names;
               else
                   Names{k} = obj.Coordinates{k}.Alias;
               end
           end
        end
    end
end

