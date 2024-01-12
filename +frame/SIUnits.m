classdef SIUnits
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
        % SI base units
        % NAME
        % (unitSymbol,latexUnitSymbol,dimensionSymbol,latexDimensionSymbol,quantityName)
        SECOND      ('s'    ,'s'    ,'T','T'    ,'time')
        METRE       ('m'    ,'m'    ,'L','L'    ,'length')
        GRAM        ('g'    ,'g'    ,'M','M'    ,'mass')
        KILOGRAM        ('kg'    ,'kg'    ,'M','M'    ,'mass')
        AMPERE      ('A'    ,'A'    ,'I','I'    ,'electric current')
        KELVIN      ('K'    ,'K'    ,'O','\Theta','thermodynamic temperature')
        MOLE        ('mol'  ,'mol'  ,'N','N'    ,'amount of substance')
        CANDELA     ('cd'   ,'cd'   ,'J','J'    ,'luminous intensity')
        
        % Base unit aliases
        METER (SIUnits.METRE)
        
        % Special derived units
        % NAME (symbol, latexSymbol, baseUnits, latexBaseUnits, otherRep,
        % description) '$\frac{a^2 \cdot c^3}{b}$
        RADIAN          ('rad'  ,   'rad',  'm*m^-1',               '$\frac{m}{m}$',                        'dimensionless',  'plane angle')
        STERADIAN       ('sr'  ,   'sr',    'm^2*m^-2',             '$\frac{m^2}{m^2}$',                    'dimensionless',  'solid angle')
        HERTZ           ('Hz'  ,   'Hz',    's^-1',                 '$\frac{1}{s}$',                        '',           'frequency')
        NEWTON          ('N'  ,   'N',      'kg*m*s^-2',            '$\frac{kg \cdot m}{s^2}$',             '',           'force, weight')
        PASCAL          ('Pa'  ,   'Pa',    'kg*m^-1*s^-2',         '$\frac{kg}{m \cdot s^2}$',             'N/m^2',      'pressure, stress')
        JOULE           ('J'  ,   'J',      'kg*m^2*s^-2',          '$\frac{kg \cdot m^2}{s^2}$',           'N*m',        'energy, work, heat')
        WATT            ('W'  ,   'W',      'kg*m^2*s^-3',          '$\frac{kg \cdot m^2}{s^3}$',           'J/s',        'power, radiant flux')
        COULOMB         ('C'  ,   'C',      's*A',                  '$s \cdot A$',                          'dimensionless',  'electric charge')
        VOLT            ('V'  ,   'V',      'kg*m^2*s^-3*A^-1',     '$\frac{kg \cdot m^2}{s^3 \cdot A}$',   'J/C',        'electric potential')
        FARAD           ('F'  ,   'F',      'kg^-1*m^-2*s^4*A^2',   '$\frac{s^4 \cdot A^2}{kg \cdot m^2}$', 'C/V',        'capacitance')
        OHM             ('ohm'  ,   '\Omega',  'kg*m^2*s^-3*A^-2',  '$\frac{kg \cdot m^2}{s^3 \cdot A^2}$', 'V/A',        'resistance, impedance, reactance')
        SIEMENS         ('S'  ,   'S',      'kg^-1*m^-2*s^3*A^2',   '$\frac{s^3 \cdot A^2}{kg \cdot m^2}$', '\Omega^-1',  'electrical conductance')
        WEBER           ('Wb'  ,   'Wb',    'kg*m^2*s^-2*A^-1',     '$\frac{kg \cdot m^2}{s^2 \cdot A^1}$', 'V*s',        'magnetic flux')
        TESLA           ('T'  ,   'T',      'kg*s^-2*A^-1',         '$\frac{kg}{s^2 \cdot A^1}$',           'Wb/m^2',     'magnetic flux density')
        HENRY           ('H'  ,   'H',      'kg*m^2*s^-2*A^-2',     '$\frac{kg \cdot m^2}{s^2 \cdot A^2}$', 'Wb/A',       'inductance')
        DEGREECELSIUS   ('degrees Celsius'  ,   '\degree C',  'K',  '$K$',                                  '',           'temperature')
        LUMEN           ('lm'  ,   'lm',    'cd*sr',                '$cd \cdot sr',                         'cd*sr',      'luminous flux')
        LUX             ('lx'  ,   'lx',    'm^-2*cd',              '$\frac{cd}{m^2}$',                     'lm/m^2',     'illuminance')
        BECQUEREL       ('Bq'  ,   'Bq',    's^-1',                 '$\frac{1}{s}$',                        '',           'radioactivity (decays per unit time)')
        GRAY            ('Gy'  ,   'Gy',    'm^2*s^-2',             '$\frac{m^2}{s^2}$',                    'J/kg',       'absorbed dose')
        SIEVERT         ('Sv'  ,   'Sv',    'm^2*s^-2',             '$\frac{m^2}{s^2}$',                    'J/kg',       'equivalent dose')
        KATAL           ('kat'  ,   'kat',  'mol*s^-1',             '$\frac{mol}{s}$',                      '',           'catalytic activity')
    end
end