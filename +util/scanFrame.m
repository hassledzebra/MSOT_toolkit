% SCANFRAME class definition for scan frames
%     
% DESCRIPTION: scanFrame is a class which acts as a container for the
% metadata corresponding to a single frame of MSOT data. If the frame was
% not pre-averaged, each frame corresponds to a single laser shot from the
% MSOT. 
%
% USAGE: In its CURRENT incarnation, simply exists as a convenient handle
% class for referencing scan metadata. There is no need to copy when
% manipulating the data, as it should always be read-only. In the future,
% the Frame class will be extended and abstracted to be useful for the
% recon and MSP classes. The end goal will be to make a memory-mapped
% interface into an MSOT study/scan in order to allow better data usage.
% 


classdef scanFrame < handle 
    
    
    properties
        % For locating each frame in the file.
        IDOffset@uint32
        USOffset
        
        % Dimensional information. 
        Wavelength
        ZPos
        ZPosReal
        Timestamp
        RelTime
        
        % Environmental/Engineering information. 
        Temperature
        DiodeReadout
        LaserEnergy
        CorrectionFactor
        
        % Information on organizing the data tree.
        Run
        Repetition
        ZNum@uint16
        WLNum@uint16
        ShotNum@uint16
        
        % Handle to parent metadata.
        meta@util.msotData
        
        
        FrameNum@uint32
        
        % For correction down the line. 
        WaterAbsCoefficient
    end
    
    methods
        
        
        function S = struct(obj)
            nFrame = numel(obj);
            
            if nFrame == 1
                S = builtin('struct',obj);
                S = rmfield(S,'meta'); 
            else
                S = arrayfun(@(X) rmfield(builtin('struct',X),'meta'),obj);
            end
        end
        
    end
    
end

