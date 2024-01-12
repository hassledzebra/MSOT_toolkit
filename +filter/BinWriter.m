classdef BinWriter < matlab.System
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
        MetaArray@struct
        
        FrameDimensions;
        FrameFormat;
        N_Frames;
        
        timeToWriteArray;
        AutoIncrement = false;
    end
    
    properties(DiscreteState)
        WriteIndex
    end
    
    % Pre-computed constants
    properties(Access = private)
        
        
    end
    
    
    methods
        function obj = BinWriter(varargin)
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
%             class(obj.MetaArray)
            obj.WriteIndex = 1;
            % Create the names for the objects.
            obj.BinFileName = fullfile(obj.FilePath,[obj.FileName,'.bin']);
            obj.MetaFileName = fullfile(obj.FilePath,[obj.FileName,'.meta']);
            
            if ~isdir(obj.FilePath)
                % mkdir(obj.FilePath);
                disp(obj.FilePath)
                mkdir obj.FilePath % to compatible with latest matlab Zheng Han,zhenghan2021@gmail.com
            end
            
            
            
            % Create the handle to the binary file and the matfile object.
            disp(obj.BinFileName)
            testFID = fopen(obj.BinFileName,'r');
            if testFID==-1
            binFileID = fopen(obj.BinFileName,'w+');
            display(binFileID)
            %             obj.MetaFile = matfile(obj.MetaFileName,'Writable',true);
            
            % Write empty values to the binary file to initialize it.
            fprintf('Frame dimensions: %d',obj.FrameDimensions(:)');
            fprintf('NFrames: %d',obj.FrameDimensions(:)');
            writeDims = [obj.FrameDimensions(:)' obj.N_Frames]
            binFileSize = 0;
            MAX_CHUNK_SIZE = 1E6;
            datatype = obj.FrameFormat;
            
            while binFileSize < prod(writeDims)
                to_write = min(prod(writeDims)-binFileSize,MAX_CHUNK_SIZE);
                binFileSize = binFileSize + fwrite(binFileID, zeros(to_write,1), datatype);
            end
            
            fclose(binFileID);
            else
            end
            
            obj.BinFileMap = memmapfile(obj.BinFileName,'Format',{obj.FrameFormat,[obj.FrameDimensions(:)'],'Data'},'Writable',true,'Repeat',obj.N_Frames);
            %             obj.MetaFile   = matfile(obj.MetaFileName,'Writable',true);
            %             obj.MetaArray = {};
            obj.MetaFile = fopen(obj.MetaFileName,'w+');
            fprintf(obj.MetaFile,'[');
        end
        
        function [y,writeLog] = stepImpl(obj,u,k)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            switch nargin
             case 2
                k=0;
                otherwise
            end
            try
                if u.IterationIndex ~=0
                    k = u.vargs;
                end
            catch
            end
            if iscell(k)
                k = k{1};
            end
%             u
%             k
            y = 0;
            writeLog = {};
            if k==0 % Assume they're trying to do an incremental write.
                locWriteInd = obj.WriteIndex
%                 obj
                obj.BinFileMap.Data(locWriteInd).Data = u.Data;
                %                obj.MetaFile.Meta(locWriteInd) = u.Meta;
%                 u
%                 k
%                 obj.MetaArray
                if ~isempty(obj.MetaArray)
                    obj.MetaArray(locWriteInd) = struct(u.Meta);
                else
                    structEntries = fieldnames(struct(u.Meta))';
                    structEntries{2,1} = cell(1,obj.N_Frames);
%                     class(obj.MetaArray)
                    testStruct = struct(structEntries{:})
%                     class(testStruct)
                    obj.MetaArray = testStruct;
                    
                    obj.timeToWriteArray = zeros(obj.N_Frames,1);
                    
                    
                    obj.MetaArray(locWriteInd) = struct(u.Meta);
                end
                obj.timeToWriteArray(locWriteInd) = toc;
                
                obj.WriteIndex = obj.WriteIndex + 1;
            else
                tic;
                obj.BinFileMap.Data(k).Data = u.Data;
                
                % Write the array contents, or create the array. 
                if ~isempty(obj.MetaArray) 
                    obj.MetaArray(k) = struct(u.Meta);
                else
                    structEntries = fieldnames(struct(u.Meta))';
                    structEntries{2,1} = cell(1,obj.N_Frames);
                    obj.MetaArray = struct(structEntries{:});
                    
                    obj.timeToWriteArray = zeros(obj.N_Frames,1);
                    
                    
                    obj.MetaArray(k) = struct(u.Meta);
                end
                obj.timeToWriteArray(k) = toc;
            end
            
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
        
        function releaseImpl(obj)
            % Release resources, such as file handles
            disp('Released in BinWriter');
            fprintf(obj.MetaFile,'%s',jsonencode(obj.MetaArray));
            fprintf(obj.MetaFile,']');
        end
        
        function numIn = getNumInputsImpl(obj)
            if obj.AutoIncrement
                numIn = 1;
            else
                numIn = 2;
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s

            % Set private and protected properties
            % obj.myproperty = s.myproperty; 

            % Set public properties and states
%             s
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end

        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
%             obj
            s = saveObjectImpl@matlab.System(obj);
            s.MetaArray = obj.MetaArray;
            % Set private and protected properties
            %s.myproperty = obj.myproperty;
        end
        
    end
end
