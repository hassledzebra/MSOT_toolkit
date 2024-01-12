classdef SIPrefixes
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
        % NAME (codeSymbol,latexSymbol,rawFactor,exponent)
        YOTTA ('Y' ,'Y'  ,1E24  ,24)
        ZETTA ('Z' ,'Z'  ,1E21  ,21)
        EXA   ('E' ,'E'  ,1E18  ,18)
        PETA  ('P' ,'P'  ,1E15  ,15)
        TERA  ('T' ,'T'  ,1E12  ,12)
        GIGA  ('G' ,'G'  ,1E9   ,9)
        MEGA  ('M' ,'M'  ,1E6   ,6)
        KILO  ('k' ,'k'  ,1E3   ,3)
        HECTO ('h' ,'h'  ,1E2   ,2)
        DECA  ('da','da' ,1E1   ,1)
        DECI  ('d' ,'d'  ,1E-1  ,-1)
        CENTI ('c' ,'c'  ,1E-2  ,-2)
        MILLO ('m' ,'m'  ,1E-3  ,-3)
        MICRO ('u' ,'\mu',1E-6  ,-6)
        NANO  ('n' ,'n'  ,1E-9  ,-9)
        PICO  ('p' ,'p'  ,1E-12 ,-12)
        FEMTO ('f' ,'f'  ,1E-15 ,-15)
        ATTO  ('a' ,'a'  ,1E-18 ,-18)
        ZEPTO ('z' ,'z'  ,1E-21 ,-21)
        YOCTO ('y' ,'y'  ,1E-24 ,-24)
    end
end