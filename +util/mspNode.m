classdef mspNode
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name            % MSP Name. Currently defaults to an indexed description of the method.
        Comment         % Comment. Current defaults to the same as Name.
        Method          % Unmixing method
        CreationTime    %** Timestamp of when the MSP was performed.
        GUID            % GUID of the MSP. Should be almost surely unique.
        ReconGUID       % GUID of the recon the MSP used. 
        Wavelengths     % Array of wavelength values used by the MSP
        InputSpectra    % In original meta, just the names of the spectra. 
                        % In modified, maybe include the spectra
                        % themselves?
        ReconNodeID     % Index of the recon node used for the MSP
        ReconMethod     % Reconstruction method used
        Resolution      % Nx and Ny
        Projections     % Number of projections simulated in the recon.
        ROI             % Size, in m, of the region of interest.
        ZPositions      % Z positions that the images were taken at.
        ZNum            % Indices of the Z positions used, (relative to the Z list in the recon node?
        
        Slices@ struct vector = struct(...
             'HasErrors',{},...         %? Unsure, might be a flag for convergence issues?
             'ReconFrames',{},...       % 0-ordinal list of which recon frames were used in the MSP.
             'ReconFrameNumbers',{},... % 1-ordinal list of same. Might have different purpose?
             'ScanFrames',{},...        % Which scan frames have condensed into this MSP. Used if reconframes averaged.
             'ZPos',{},...              % Z position to 2 decimal places
             'ZPosReal',{},...
             'ZNum',{},...              % Z position in index form.
             'Run',{},...               % Run coordinate of that slice
             'Repetition',{},...        % Rep coordinate of that slice
             'RelTime',{},...           % Time, in seconds, that the first frame in that dataset was taken.
             'Timestamp',{},...         %** Timestamp corresponding to the absolute time of the first frame. 
             'Components',struct('IDOffset',{},...
                                 'ComponName',{},...
                                 'ComponRef',{},...
                                 'Invert',{},...
                                 'SpectrumName',{}));          %? Components structure, effectively just a listing of the components and how to refer to them.
         
        Structure       % Structured array of linear indices into the MSP dataset
        SliceIndex      % Index list corresponding to which Z slice is used. Helps for determining where each component image came from.
        
        %% Things in the .msot but not in the loadMSOT metadata
        DiscardNegativeValues
        Reapplied
        DefaultBackground
        ICAComponents
        
    end
    
    methods
        
        %% Constructor
        function obj = mspNode(parentMeta,mspNode,mspNodeID)
            
            if nargin == 0
                obj.Name='null';
                return;
            end
            
            % Nested functions to read XML stuff.
            function textElement = getTextElement(javaTagString)
                textElement=char(mspNode.getObjectValue.get_store.getElementsByTagName(javaTagString).item(0).fetch_text(0));
            end
            
            function numElement = getNumElement(javaTagString)
                numElement=str2double(getTextElement(javaTagString));
            end
            
            function boolElement = getBoolElement(javaTagString)
                textBool=char(getTextElement(javaTagString));
                
                if strcmpi(textBool,'true')
                    boolElement=true;
                else
                    boolElement=false;
                end
            end
            
            %% Get the top-level data.
            obj.Name=char(mspNode.getName);
            obj.Comment=char(mspNode.getComment);
            obj.CreationTime=char(mspNode.getCreationTime);
            obj.DefaultBackground=(mspNode.getDefaultBG.getIntValue);
            obj.GUID=char(mspNode.getGUID);
            obj.ICAComponents=(mspNode.getIcaComponents);
            obj.ReconGUID=char(mspNode.getReconGUID);
            for tr = 1:numel(parentMeta.ReconNode)
               if strcmp(obj.ReconGUID,parentMeta.ReconNode(tr).GUID)
                  obj.ReconNodeID = tr; 
               end
            end
            obj.Method=char(mspNode.getMethod);
%             obj.ReconNodeID=find(strcmp(obj.ReconGUID,{parentMeta.ReconNode(:).GUID}));
            obj.Reapplied=getBoolElement(java.lang.String('Reapplied'));
            obj.DiscardNegativeValues=getBoolElement(java.lang.String('DiscardNegativeValues'));
            obj.InputSpectra=string(mspNode.getInputSpectra.getStringArray);
            
            
            mspSlices=mspNode.getSlicesFrame.getMspSliceFrameArray;
            sliceArrayLength=size(mspSlices,1);
            
            % Build the msp slice structure
            %% TODO: Make sure that this functions as expected.
            
            for sliceIndex=1:sliceArrayLength
               msIndex=mspSlices(sliceIndex).getReconFrameID.getIntArray;
               sfs=[parentMeta.ReconNode(obj.ReconNodeID).Frames(msIndex+1)];
               tempStruc = struct('HasErrors',{mspSlices(sliceIndex).getHasErrors},...
                                     'ReconFrames',msIndex,...
                                     'ReconFrameNumbers',msIndex+1,...
                                     'ScanFrames',reshape([sfs.ScanFrames],[],numel(msIndex)),... % Have to index into the recon node for this. 
                                     'ZPos',parentMeta.ReconNode(obj.ReconNodeID).Frames(msIndex(1)+1).ZPos,...
                                     'ZPosReal',parentMeta.ReconNode(obj.ReconNodeID).Frames(msIndex(1)+1).ZPosReal,...
                                     'ZNum',parentMeta.ReconNode(obj.ReconNodeID).Frames(msIndex(1)+1).ZNum,...
                                     'Run',parentMeta.ReconNode(obj.ReconNodeID).Frames(msIndex(1)+1).Run,...
                                     'Repetition',parentMeta.ReconNode(obj.ReconNodeID).Frames(msIndex(1)+1).Repetition,...
                                	 'RelTime',parentMeta.ReconNode(obj.ReconNodeID).Frames(msIndex(1)+1).RelTime,...
                                     'Timestamp',parentMeta.ReconNode(obj.ReconNodeID).Frames(msIndex(1)+1).Timestamp,...
                                     'Components',struct('IDOffset',arrayfun(@(x) x.getIDOffset,mspSlices(sliceIndex).getComponentList.getMspSliceComponentArray,'UniformOutput',false),...
                                                         'ComponName',arrayfun(@(x) char(x.getComponName),mspSlices(sliceIndex).getComponentList.getMspSliceComponentArray,'UniformOutput',false),...
                                                         'ComponRef',arrayfun(@(x) x.getComponRef,mspSlices(sliceIndex).getComponentList.getMspSliceComponentArray,'UniformOutput',false),...
                                                         'Invert',arrayfun(@(x) x.getInvert,mspSlices(sliceIndex).getComponentList.getMspSliceComponentArray,'UniformOutput',false),...
                                                         'SpectrumName',arrayfun(@(x) char(x.getSpectraName),mspSlices(sliceIndex).getComponentList.getMspSliceComponentArray,'UniformOutput',false))); %mspSlices(1).getComponentList.getMspSliceComponentArray(1).getComponName
                obj.Slices(sliceIndex)=tempStruc;
            end
            
            nComponents=numel(obj.Slices(1).Components);
            nRep=numel(unique([obj.Slices(:).Repetition]));
            nZ=numel(unique(round([obj.Slices(:).ZPosReal]*1E4)/1E4));
            nRun=numel(unique([obj.Slices(:).Run]));
            %
            obj.Wavelengths=arrayfun(@(x) x.getIntValue,mspNode.getRelatedWavelengths.getWavelengthArray);
            obj.InputSpectra=[]; 
            obj.ReconMethod=parentMeta.ReconNode.Method;
            obj.Resolution=parentMeta.ReconNode.Resolution;
            obj.Projections=parentMeta.ReconNode.Projections;
            obj.ROI=parentMeta.ReconNode.ROI;
            obj.ZPositions=parentMeta.ReconNode.ZPositions;
            obj.ZNum=parentMeta.ReconNode.ZNum; % Does not match loadMSOT order because of reconNode
            obj.Structure=permute(reshape(1:(sliceArrayLength*nComponents),[nComponents nRep nZ nRun]),[4 3 2 1]);
            tempSliceIndex=repmat(1:sliceArrayLength,[3,1]);
            obj.SliceIndex=tempSliceIndex(:)';
            
            clear pmet;
        end
        %%% End Constructor %%%
        
        %% Ancestor processing. Finds the handles which point to the metadata of the reconstruction nodes that went into each given MSP slice.
        function ancestorMeta=getReconAncestorMeta(mspIndex)
            
        end
        
        function ancestorMeta=getSignalAncestorMeta(mspIndex)
            
        end
        
        function S = struct(obj)
            nNode = numel(obj);
            if nNode == 1
                S = builtin('struct',obj);
            else
                S = arrayfun(@(X) builtin('struct',X),obj);
            end
        end
        
    end %methods
    
end %classdef

