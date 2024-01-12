classdef Coordinate
    %Coordinate Describes the interpretation of a given coordinate when
    %embedded in a coordinate system.
    %
    % For example, one may have a coordinate t_s, where there are N_t
    % entries, with each entry containing a single time point.
    %
    % One may also have a coordinate s, where there are N_s entries, but
    % each entry contains TWO spatial coordinates.
    %
    % The joint use of multiple coordinates simultaneously allows the
    % description of a data tensor according to its coordinate embedding,
    % which further allows its relationship to other data tensors.
    % As a technical note, because of the lack of algebraic structure, the
    % data is actually a holor, not a tensor.
    
    
    properties
        Values

        Alias
        Space % Which of the standard spaces the data exists in. 
        
        Names % Labels the coordinate. e.g. "Sampling Time"
             % Can have multiple values, but requires a similar
             % number of values for the Name_Range.
             % Can have names for each entry's values, as well.
             

          
        Unit % Labels each of the entries. Must have the same number of values as Name.
          Unit_latex % Provides the rendering for representing the unit using LaTeX renderers.
        
        VertexOrEdge; % Describes whether the value is associated with a point on the grid or with an edge betwen points.
                      % An image is the joint of two edge-oriented
                      % coordinates, and so the face is the association
                      % with the value. 
                      %
                      % Can be derived from the number of entries relative
                      % to the N_values, or vice-versa. e.g. the coordinate
                      % system for an image is properly a (Nx+1)x(Ny+1)
                      % system where each pixel's corners coordinates are
                      % defined. Allows for things like integrating over
                      % the pixel.
                      

        Orientation   % Orientation of the axis relative to the values. This is largely only important for edge-
                      % associated coordinates, to interpret the 'direction
                      % of change'.
                      
        
    end
    properties (Access = private)
       IndexCoordinates 
       Values_
       Children_Names ={};
       Children_Alias ={};
    end
    
    methods
        function obj = Coordinate(varargin)
            %Coordinate Constructor
            cBool = cellfun(@(x) isa(x,'frame.Coordinate'),varargin);
            if ~isempty(cBool)&&all(cBool)
                propertyList = properties(varargin{1});
                for k = 1:numel(propertyList)
                    propName = propertyList{k};
                    switch propName
                        case {'Alias','Names'}
                            obj.(propName) = strjoin(cellfun(@(x) x.(propName),varargin,'uni',false),'(+)');
                        case 'Values'
                            valueLists = cellfun(@(x) x.(propName),varargin,'uni',false);
                            valListSizes = cellfun(@size,valueLists,'uni',false);
                            valListNumel = cellfun(@numel,valueLists);
                            assert(all(valListNumel==valListNumel(1)));
                            assert(all(cellfun(@(x) isequaln(x,valListSizes{1}),valListSizes)));
                            
                            NCoordinates = valListNumel(1);
                            obj.IndexCoordinates = 1:NCoordinates;
                            obj.Values_ = valueLists;
                            
                        otherwise
                        obj.(propName) = cellfun(@(x) x.(propName),varargin,'uni',false);
                    end
                    
                end
            elseif ~isempty(cBool)&&any(cBool)
                nCoords = sum(cBool);
                % Join together the coordinates.
                propertyList = properties(varargin{1});
                for k = 1:numel(propertyList)
                    propName = propertyList{k};
                    switch propName
                        case {'Alias','Names'}
                            obj.(propName) = strjoin(cellfun(@(x) x.(propName),varargin(1:nCoords),'uni',false),'(+)');
                            obj.(['Children_',propName]) = cellfun(@(x) x.(propName),varargin(1:nCoords),'uni',false);
                        case 'Values'
                            valueLists = cellfun(@(x) x.(propName),varargin(1:nCoords),'uni',false);
                            valListSizes = cellfun(@size,valueLists,'uni',false);
                            valListNumel = cellfun(@numel,valueLists);
                            assert(all(valListNumel==valListNumel(1)));
                            assert(all(cellfun(@(x) isequaln(x,valListSizes{1}),valListSizes)));
                            
                            NCoordinates = valListNumel(1);
                            obj.IndexCoordinates = 1:NCoordinates;
                            obj.Values_ = valueLists;
                            
                        otherwise
                        obj.(propName) = cellfun(@(x) x.(propName),varargin(1:nCoords),'uni',false);
                    end
                    
                end
                % Add the other props.
                for k = (nCoords+1):2:numel(varargin)
                        obj.(varargin{k}) = varargin{k+1}; 
                end
            elseif (nargin > 0) && ( ~isa( varargin{1} , 'struct' ) )
                
                for k = 1:2:numel(varargin)
                    if strcmp(varargin{k},'Values')
                       obj.Values_ = varargin{k+1}(:); 
                       NCoordinates = numel(obj.Values_);
                       obj.IndexCoordinates = 1:NCoordinates;
                    else
                        obj.(varargin{k}) = varargin{k+1}; 
                    end
                end
            % Otherwise, if there's only one input and it's a struct,
            % assume that it contains the fields we need to create the
            % filter.
            elseif (nargin==1) && (isa(varargin{1},'struct'))
                structFieldNames = fieldnames(varargin{1});
                structFieldValues = struct2cell(varargin{1});
                objProps = properties(obj);
                [~,keepInds] = ismember(objProps,structFieldNames);
                
                concatStruct = [structFieldNames(keepInds(keepInds~=0)),structFieldValues(keepInds(keepInds~=0))].';
                setProperties(obj,numel(concatStruct),concatStruct{:})
            else
                
                % Nothing should happen in other cases. 
            end
        end
        
        function outputArg = numel(obj)
            outputArg = numel(obj.Values);
        end
        
        function obj = interpCoords(obj,N)
            obj.Values_ = cellfun(@(x) interp1(1:numel(x),x,linspace(1,numel(x),N)),obj.Values_,'uni',false);
            if ~isempty(obj.IndexCoordinates)
                Ncoords = numel(obj.IndexCoordinates);
                obj.IndexCoordinates = interp1(1:Ncoords,1:Ncoords,linspace(1,Ncoords,N));
            end
        end
        
        function Values = get.Values(obj)
           if iscell(obj.Values_) 
               Values = obj.IndexCoordinates;
           else
               Values = obj.Values_;
           end
            
        end
        function Values = getChildValues(obj)
            
            Values = obj.Values_;
        end
        
        function coordVals = getCoordinateValues(obj,coordinateStandardName)
            
           coName = strcmp(obj.Children_Names,coordinateStandardName);
           if any(coName) && sum(coName)==1
%                targCoord = obj.Coordinates{coName};
%                if iscell(targCoord.Space)&&numel(targCoord.Space)==1 || ~iscell(targCoord.Space)
%                 coordVals = targCoord.Values;
%                else
%                 coordVals = targCoord.getChildValues;
%                end
                coordVals = obj.Values_{coName};
           end
           
           
        end
    end
end

