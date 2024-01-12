classdef H5Writer < matlab.System
    %H5LOADER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FilePath;
        FileName;
        
        BinFileMap; % Here this will correspond more or less to the H5D identifier with some write() semantics.
        BinFileSize; % Not actually that useful?
        BinFileName; % The BinFileName here would probably more accurately be represented as the DATA group. 
        
        MetaFile; % H5D identifier into the vlen_string array corresponding to each JSON object
        MetaFileSize; % Again, not that useful, especially with the strings being as polymorphic as they are.
        MetaFileName; % META group, instead
        MetaArray@struct % Might not actually be necessary to have, since we can write to file directly by flushing the buffer.
            % Should probably add something describing the chunking
            % settings, as well. 
            % Further, should add coordinates and mapping histories.
            % Realistically the coordinate and everything could basically
            % be added during setup, since it should be consistent
            % throughout the entire dataset.
            %
            
        
        FrameDimensions;
        FrameFormat;
        N_Frames;
        
        timeToWriteArray;
        AutoIncrement = false;
    end
    
    methods
        function obj = H5Writer(varargin)
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
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end

    methods(Access = protected)
        function setupImpl(obj)
            %% Create the H5 file 
                % Check if the file already exists.
                % If it does, post a warning and overwrite it.
                    % TODO: Add write/w+ stuff. 
            
                    
            %%
               % Create the record structure
               % DATA - ref to the frame's data.
               % META - ref to the frame's meta JSON string.
               % COORDINATE - ref to the frame's coordinate system.
                    % NOTE: Should probably store coordinate systems the
                    % first time they're seen and do an object comparison
                    % thereon. Would probably be slow, however. Since each
                    % coordinate system is a handle, however, MATLAB can
                    % just do handle comparisons to see if it's the same
                    % object in memory, which would be quite fast.
                    % Remember that there are also datasets like the
                    % single-wavelength images where their embedded
                    % dimensions basically become properties (e.g. Z
                    % position, wavelength). Would still probably be pretty
                    % cheap.
                    % 
                    % 
               % HISTORY - (mothballed) Tags from mappings that describe
               % each transformation that the frame underwent.
               % MAPPINGS - (mothballed) Similar to the Coordinate, the
               % 'recipe' for each of the processes should be recorded.
               % Realistically this will need to be the input from the
               % parsed JSON commands which describe the overall process.
               % Could actually just be the baked recipe file itself.
        end
        
        function stepImpl(obj,u,k)
            
            % The chunking of the dataset matters a lot here, and we need
            % to record those internal coordinate systems as well (IF the
            % setting is set)
            
            % Receive frame and index
            % Format meta into a JSON string.
            
            % NOTE: The memory dataspace, datatype, etc. should not change
            % between frames, so we can just cache those as private
            % properties.
            
            % Copy into the next data slot.
            % Copy into meta JSONs.
            % Get reference to dataslot
            % Get reference to meta entry
            % Get reference to coordinate system dump.
            
            % Construct record with pointers to each frame's data, meta,
            % and coordinate system.
            
            % NOTE::: THE DATASET OBJECT ALSO HAS METADATA WHICH CAN BE
            % ATTACHED VIA A META ATTRIBUTE. Somewhat similar to the
            % mapping, only the mapping is less well defined here because
            % of information lack. This should be an unlimited data structure
            % representing a bunch of metadata which were added based on
            % each of the mappings/transitions, so that DATASET.META{1}
            % would represent msotData or acquisition structures,
            % DATASET.META{2} would represent recon information.
            
            % The frames could also have this sort of structure, but
            % their metadata would tend to have trouble stacking with
            % additional information, hence the concept of transformation
            % and having different frames with different metadata in
            % different places. 
            % How to flush buffers so there's not a lot of wasted time?
            
            
        end
    end
end

