classdef BinLoader < matlab.System
    % BinLoader loads binary data into DataFrames.

    % Public, tunable properties
    properties
        
        FilePath;
        FileName;
        
        BinFileMap;
        BinFileSize;
        BinFileName;
        
        MetaFile; % Matfile object for writing the unstructured metadata.
        MetaFileSize;
        MetaFileName;
        
        FrameDimensions;
        FrameFormat;
        N_Frames;

    end

    properties(DiscreteState)
        LoadIndex
    end

    % Pre-computed constants
    properties(Access = private)
        isMetaMatFile = false;
    end

    methods 
        function obj = BinLoader(varargin) 
            if nargin>0 && ~isa(varargin{1},'struct') 
                setProperties(obj,nargin,varargin{:},'FilePath','FileName','FrameDimensions','FrameFormat','N_Frames')
                
            % Otherwise, if there's only one input and it's a struct,
            % assume that it contains the fields we need to create the
            % filter.
            elseif nargin==1 && isa(varargin{1},'struct')
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
    end
    methods(Access = protected)
        function setupImpl(obj)
            obj.LoadIndex = 1;
            
            obj.BinFileName = fullfile(obj.FilePath,[obj.FileName,'.bin']);
            obj.MetaFileName = fullfile(obj.FilePath,[obj.FileName,'.meta']);
            datatype = obj.FrameFormat;
            
            obj.BinFileMap = memmapfile(obj.BinFileName,'Format',{datatype,[obj.FrameDimensions(:)'],'Data'},'Writable',false,'Repeat',obj.N_Frames);
             try
                 tempMeta = matfile(obj.MetaFileName,'Writable',false);
                 obj.MetaFile   = struct('Meta',tempMeta.Meta);
             catch
                 obj.isMetaMatFile = false;
                 obj.MetaFile.Meta = jsondecode(fileread(obj.MetaFileName));
             end
           
            % Append the opening bracket on the meta array.
        end

        function [y,loadLog] = stepImpl(obj,k)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            loadLog = {};
            if k==0 % Assume they're trying to do an incremental load. 
               locLoadInd = obj.LoadIndex;
                y.Data = obj.BinFileMap.Data(locLoadInd).Data;
%                 y.Meta = obj.MetaFile(locLoadInd);
                 y.Meta = obj.MetaFile.Meta(locLoadInd);
               
               obj.LoadIndex = obj.LoadIndex + 1;
            else
                 y.Data = obj.BinFileMap.Data(k).Data;
%                  y.Meta = obj.MetaFile(k);
                 y.Meta = obj.MetaFile.Meta(k);
            end
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end

        function releaseImpl(obj)
            % Release resources, such as file handles.
            % Append the closing bracket on the meta array.
        end
    end
end
