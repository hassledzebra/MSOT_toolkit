classdef DataFrame < handle
    %DataFrame describes the general form of a data object used in the
    %processing pipeline, conformed to a structure that allows for descriptions
    %of filtering and transformation. 
    % The mathematical form that this takes is of a tensor with a specific
    % sensibility. For example, a 3-D image is (X x Y x Z), whereas Z stacks of
    % 2-D images is (X x Y) x Z . Based on whether each coordinate is face- or
    % edge-centered (See below), each 'coordinate point' on the resulting tensor
    % can either be on the 
    
    properties
        
        % Data properties
        Data 
            % The actual payload of a DataFrame. As the information contained
            % within a DataFrame requires interpretation, this should have
            % metadata accompanying it in order to describe it. As a good
            % example, we should be able to have two DataFrame objects
            % describing bijected data, one in 3D pixel coordinates and one in
            % WL-DCT conjugate coordinates.
             
        Data_Dimensions
            % A cell array of strings which contains the names for each of the
            % dimensions in a frame. 
            % Data_Dimensions = {"Y","X","Z","C","T"} for example. 
            % THE NUMEL OF THIS PROPERTY DETERMINES THE NUMBER OF DIMENSIONS, NOT
            % THE DATA SIZE ITSELF - MATLAB DOES NOT CONSIDER TRAILNG SINGLETON
            % DIMENSIONS.
            
        Data_Coordinates 
            % A heterogenous cell array describing the coordinates for each of
            % the dimensions in the data. The number of elements of a given Data_Coordinates
            % entry should be the same as the number of elements in the
            % corresponding data dimension, or one higher. In the case of
            % equality, each coordinate corresponds to a label for the given
            % block (i.e., the face-centered coordinate). In the case of N+1,
            % each coordinate is assumed to label the BOUNDARY between
            % sequential blocks (i.e., the edge-centered coordinate). 
            %
            % In the case of categorical data, the entry should be a cell array
            % of strings. This allows for metadata association as a tensor
            % product.
            %
            % UNCLEAR HOW TO HANDLE SCALE SPACES EXCEPT AS CELL ARRAYS. 
            % 
             
        % Meta properties
        Meta 
            % Metadata can basically be anything, but we should never lose
            % information, only add to it. If the frame is changing context,
            % then we should make a copy and modify that instead.
            %
            % This should reflect the chain of transformations and processing
            % that has gone forth to create the current frame, but that is more
            % for completeness than anything else. 
            
        FrameHistory
            % Frames should generally have an idea of their provenance as they
            % continue to be processed. This cell array allows descriptions of
            % the processing procedures to be passed along as the processing
            % continues.
        
        ParentFrames 
            % From what frame(s) was this frame derived?
            % TODO: How to address e.g. transient frames like those
            % internal to the structure of a state filter?
                % Could potentially sidestep the issue by labelling the
                % process transition via the FrameHistory, i.e. the last
                % ParentFrame is the input frame to the last process in the
                % FrameHistory. 
        
    end
    
    methods
        function obj = DataFrame(dataIn,metaIn)
            %DataFrame Constructs a basic frame.
            obj.Data = dataIn;
            obj.Meta = metaIn;
        end
        
        % COPY
        
        % HISTORY (APPEND)
        
        % SAVE
        
        % LOAD
        
    end
end

