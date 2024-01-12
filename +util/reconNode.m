classdef reconNode
    %RECONNODE defines the reconstruction node structure in the MSOT data.
    
    properties
        Name@char vector        
        Comment@char vector
        GUID@char vector
        Method@char vector
        Resolution@double scalar
        ZResolution@double scalar
        Projections@uint16 scalar
        ROI@double scalar
        
        ROILow@double scalar
        ROIHigh@double scalar
        FilterLow@double scalar
        FilterHigh@double scalar
        SignalFilterType@char vector 
        TrimSpeedOfSound@uint16 scalar
        NumOfFrames@uint16 scalar
        FrameThickness@double scalar
        ROIZLow@double scalar
        ROIZHigh@double scalar
        DepthCorrection@logical scalar
        BackgroundAbsorption@double scalar
        BackgroundOxygenation@double scalar
        AxialOffset@double scalar
        LightSpotSize@double scalar
        CouplantCorrection@logical scalar
        CouplantSpeedOfSound@uint16 scalar
        CouplantCenterOfRotation@double scalar
        CouplantCurvatureRadius@double scalar
        ImpulseResponse@logical scalar
        SensitivityMap@logical scalar
        
        ZPositions@single vector
        ZNum@uint16 vector
        Runs@uint16 vector
        Repetitions@uint16 vector
        Wavelengths@double vector
        
        Frames@struct vector = struct('IDOffset',{},...
                                     'ScanFrames',{},...
                                     'Wavelength',{},...
                                     'ZPos',{},...
                                     'ZPosReal',{},...
                                     'ZNum',{},...
                                     'Run',{},...
                                     'Repetition',{},...
                                     'Timestamp',{},...
                                	 'RelTime',{},...
                                     'LaserEnergy',{},...
                                     'HasErrors',{},...
                                     'ZOffset',{});
        
        ReconStructure@double 
        RelTime@double
        Timestamps@double vector
        IDLookup@uint16 vector
    end
    
    properties (SetAccess = private)
        isLoaded = false;
    end
    
    methods
        
        %% Constructor
        function obj = reconNode(parentMeta,recNode,reconNodeID)
           
            if nargin == 0
                obj.Name='null';
                return;
            end
            
            % Nesteds for reading XML stuff.
            function textElement = getTextElement(javaTagString)
                textElement=char(xmlObject.getElementsByTagName(javaTagString).item(0).fetch_text(0));
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
            
            %% Get the top-level data
            % First the stuff that's got an msotbeans implementation
            obj.Name=char(recNode.getName);
            obj.Comment=char(recNode.getComment);
            obj.GUID=char(recNode.getGUID);
            obj.Method=char(recNode.getMethod);
            
            obj.Resolution=recNode.getResolution; % Number of pixels.
            obj.ROI=double(recNode.getRoi); %ROI field size. 
            
            obj.Projections=recNode.getProjections; 
                
            % Next the stuff that's only accessible through XML
            xmlObject=recNode.getObjectValue.get_store;
                %Text
                    obj.SignalFilterType=getTextElement(java.lang.String('SingalFilterType'));
                %Numbers
                    %doubles
                    
                    obj.ZResolution=getNumElement(java.lang.String('ZResolution'));
                    obj.ROILow=getNumElement(java.lang.String('RoiLow'));
                    obj.ROIHigh=getNumElement(java.lang.String('RoiHigh'));
                    obj.FilterLow=getNumElement(java.lang.String('FilterLow'));
                    obj.FilterHigh=getNumElement(java.lang.String('FilterHigh'));
                    obj.FrameThickness=getNumElement(java.lang.String('FrameThickness'));
                    obj.ROIZLow=getNumElement(java.lang.String('RoiZlow'));
                    obj.ROIZHigh=getNumElement(java.lang.String('RoiZhigh'));
                    obj.BackgroundAbsorption=getNumElement(java.lang.String('BackgroundAbsorption'));
                    obj.BackgroundOxygenation=getNumElement(java.lang.String('BackgroundOxygenation'));
                    obj.AxialOffset=getNumElement(java.lang.String('AxialOffset'));
                    obj.LightSpotSize=getNumElement(java.lang.String('LightSpotSize'));
                    obj.CouplantCenterOfRotation=getNumElement(java.lang.String('CouplantCenterOfRotation'));
                    obj.CouplantCurvatureRadius=getNumElement(java.lang.String('CouplantCurvatureRadius'));
                    %uint16s
                    obj.CouplantSpeedOfSound=uint16(getNumElement(java.lang.String('CouplantSpeedOfSound')));
                    obj.TrimSpeedOfSound=uint16(getNumElement(java.lang.String('TrimSpeedOfSound')));
                    obj.NumOfFrames=uint16(getNumElement(java.lang.String('NumOfFrames')));
                %Booleans
                    obj.CouplantCorrection=getBoolElement(java.lang.String('CouplantCorrection'));
                    obj.ImpulseResponse=getBoolElement(java.lang.String('ImpulseResponse'));
                    obj.SensitivityMap=getBoolElement(java.lang.String('SensitivityMap'));
                    obj.DepthCorrection=getBoolElement(java.lang.String('DepthCorrection'));
            
            %% Frames
            reconFrames=recNode.getReconFrames.getDataModelReconstructionFrameArray;
            frameArrayLength=recNode.getReconFrames.sizeOfDataModelReconstructionFrameArray;
            frameCursor=recNode.getReconFrames.newCursor;
            frameCursor.toFirstChild;
            
            % Build the recon frame node structure.
            %% TODO: Make sure that this actually works for averaged recons.
            
            obj.ReconStructure=[];  
            % This brings the scan frames structure into local scope
            % persistently, which makes the parsing a lot faster. Must
            % ensure that the pmet is different for different scans. 
            persistent pmet pmetString;
            if isempty(pmet)
                pmet=parentMeta.ScanFrames;
                pmetString=parentMeta.XMLFileName;
            else
                if ~strcmp(pmetString,parentMeta.XMLFileName)
                    pmet=parentMeta.ScanFrames;
                end
            end
            
            
            offsetHold = [pmet(:).IDOffset];
            las=[pmet.LaserEnergy];
            obj.Frames(frameArrayLength).IDOffset = frameArrayLength-1;
            frameHolder=obj.Frames;
            
            for frameIndex=1:frameArrayLength
                %% Parse the recon frame.
                
                % First parse out the scanFramesID
                frameCursor.toFirstChild;
                frameCursor.toFirstChild;
                
                frameIDindex=1;
                partialStruc.tempScanFrames(frameIDindex)=getIntValue(getObject(frameCursor));
                while(frameCursor.toNextSibling)
                    frameIDindex=frameIDindex+1;
                    partialStruc.tempScanFrames(frameIDindex)=getIntValue(getObject(frameCursor));
                end
                frameCursor.toParent;
                % Next the IDOffset
                
                frameCursor.toNextSibling;
                partialStruc.IDOffset=getIntValue(getObject(frameCursor));
                
                % Next zOffset
                frameCursor.toNextSibling;
                partialStruc.ZOffset=getIntValue(getObject(frameCursor));
                
                % Next HasErrors
                partialStruc.HasErrors=getIntValue(getObject(frameCursor));
                
                % Get out of the frame
                frameCursor.toParent;
                %%
%                 rfIndex=partialStruc.IDOffset+1;
                rfIndex=partialStruc.tempScanFrames;
                [~,accessIndex]=ismember(partialStruc.tempScanFrames,offsetHold);
%                 [~,accessIndex_]=ismember(partialStruc.tempScanFrames,[pmet.IDOffset]);
%                 assert(isequal(accessIndex,accessIndex_));
% recNode.getReconFrames.getDataModelReconstructionFrameArray(0).getScanFramesID.get_store().find_element_user(javaObject('javax.xml.namespace.QName',"","int"), 1)

                laserArray=[las(accessIndex)];
                
                %% The ZNum does not align in the same way as loadMSOT does. 
                tempFrame= pmet(accessIndex);
                tempStruc=    struct('IDOffset',partialStruc.IDOffset,...
                                     'ScanFrames',partialStruc.tempScanFrames,...
                                     'Wavelength',[tempFrame.Wavelength],...
                                     'ZPos',[tempFrame.ZPos],...
                                     'ZPosReal',[tempFrame.ZPosReal],...
                                     'ZNum',[tempFrame.ZNum],...
                                     'Run',[tempFrame.Run],...
                                     'Repetition',[tempFrame.Repetition],...
                                     'Timestamp',[tempFrame.Timestamp],...
                                	 'RelTime',[tempFrame.RelTime],...
                                     'LaserEnergy',laserArray,...
                                     'HasErrors', partialStruc.HasErrors,...
                                     'ZOffset',partialStruc.ZOffset);
                                 
                frameHolder(frameIndex)=tempStruc;
                frameCursor.toNextSibling;
%                 disp(['frame',num2str(frameIndex)])
            end
            obj.Frames=frameHolder;
            %TODO: Pare down the recon meta so that if we selected only
            %certain frames, those are the only ones reflected in the meta.
            %TODO: Allow for missing/misaligned data.
            uniqueRuns=unique([obj.Frames.Run]);
            %Timestamps of runs
            obj.Timestamps=[pmet(uniqueRuns).Timestamp];%parentMeta.Timestamps(uniqueRuns);            
            
            %TODO: Assign runs, numerical Z positions, wavelengths, ZNum, 
            obj.ZPositions=unique([obj.Frames(:).ZPos],'stable');
            obj.ZNum=unique([obj.Frames(:).ZNum],'stable'); % Does this correspond to the index of the Z relative to the original frames?
            obj.Runs=unique([obj.Frames(:).Run]); % Same question
            obj.Repetitions=unique([obj.Frames(:).Repetition]);
            obj.Wavelengths=unique([obj.Frames(:).Wavelength],'stable');
%             obj.Frames=[];
            obj.ReconStructure=[]; %bit more complicated. 
            
            
            nRun=numel(uniqueRuns);
            nZ=numel(unique([obj.Frames.ZNum]));
            nRep=numel(unique([obj.Frames.Repetition]));
            nWL=numel(unique([obj.Frames.Wavelength]));
            nFrames=numel(partialStruc.tempScanFrames);
            
            relTimeVec=[pmet([obj.Frames.IDOffset]+1).RelTime];
            %Relative time
            
            Nanvec=NaN*ones(1,prod([nRun,nZ,nRep,nWL,nFrames])-numel(relTimeVec));
            
            if nFrames>1
                obj.RelTime=reshape([relTimeVec Nanvec],[nRun,nZ,nRep,nWL,nFrames]);
                obj.ReconStructure=permute(reshape(1:prod([nRun,nZ,nRep,nWL]),[nWL,nRep,nZ,nRun]),[4 3 2 1]);
            else
%                 nf = numel(relTimeVec)./prod([nRun,nZ,nRep,nWL]);
                nf = parentMeta.ShotNum; % CHANGED
                if ~isinteger(nf)
                    nf = round(nf);
                    warning('Number of frames not an integer. May be an implementation issue. Rounding to nearest integer.');
                end
                obj.RelTime = reshape([relTimeVec Nanvec],[nRun,nZ,nRep,nWL,nf]);
                obj.ReconStructure = permute(reshape([1:prod([nRun,nZ,nRep,nWL,nf])],[nf,nWL,nRep,nZ,nRun]),[5 4 3 2 1]);
            end
            
            obj.IDLookup=uint16([obj.Frames.IDOffset]+1);
            obj.isLoaded=true;
            clear pmet;
        end
        
        
        function S = struct(obj)
            nNode = numel(obj);
            if nNode == 1
                S = builtin('struct',obj);
            else
                S = arrayfun(@(X) builtin('struct',X),obj);
            end
        end
        
        
    end
    
end

