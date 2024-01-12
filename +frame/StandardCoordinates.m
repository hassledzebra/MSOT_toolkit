classdef StandardCoordinates
    % Contains references to standard coordinate systems
    
    methods
        function obj = defaults(varargin)
           obj.contents = varargin;
        end
        
        function cels = cell(obj)
           cels = {obj.contents{:}}; 
        end
        
        function c = struct(obj)
           if isstruct(obj.contents{:})
               c = obj.contents{:};
           else
               c = obj.contents{:};
           end
        end
        
        function out = unpack(obj)
            out = obj.contents{:};
        end
        
                
    end
    
    enumeration
        CARTESIAN_X
        CARTESIAN_Y
        CARTESIAN_Z
        
        CYLINDRICAL_R
        CYLINDRICAL_THETA
        CYLINDRICAL_Z
        
        SPHERICAL_R
        SPHERICAL_PHI
        SPHERICAL_PSI
        
        
        
    end
end