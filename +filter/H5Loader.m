classdef H5Loader
    %H5LOADER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FilePath;
        FileName;
        
        BinFileMap;
        BinFileSize;
        BinFileName;
        
        MetaFile; 
        MetaFileSize;
        MetaFileName;
        
        FrameDimensions;
        FrameFormat;
        N_Frames;
    end
    
    methods
        function obj = H5Loader(inputArg1,inputArg2)
            %H5LOADER Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

