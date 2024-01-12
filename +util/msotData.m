classdef msotData < handle
    properties
        verbose@logical = false;
        progress@logical = false;
    end
    
    properties ( SetAccess = private )
        
        %% Filesystem properties. These will be rendered as absolute paths.
        XMLFileName@char vector 
        FileName@char vector
        RealPath@char vector
        StudyName@char vector
        StudyPath@char vector
        
        
        DOM % Document object model, the main XML object that we load. This seems to cause a number of problems with the heap memory in larger scans. 
        
        Name@char vector
        Comment@char vector
        DataModelVersion@char vector
        FriendlyName@char vector
        FolderName@char vector
        Scientist@char vector
        CreationTimeTxt@char vector
        CreationTime@double scalar
        Complete@logical scalar = true;        % legacy
        is3D@logical scalar = false;
        is2D@logical scalar = false;
        TransducerType@char vector;
        
        MeasurementDesc@struct scalar = struct();
        HWDesc@struct scalar = struct();
        USpresent@logical scalar = false;
        US@struct scalar;
        OAMPreset@struct scalar = struct();
        RepNum@uint16 scalar;
        RunNum@uint16 scalar;
        ShotNum@uint16 scalar;
        AverageTemperature@single scalar;
        
        ScanFrames@util.scanFrame vector;
        ScanStructure@single;
        RelTime@double;
        %         LaserEnergy@single;
        LaserEnergy@double;
        DiodeReadout@single;
        %         CorrectionFactor@single;
        CorrectionFactor@double;
        Temperature@single;
        Wavelengths@double vector;
        ZPositions@single vector;
        Timestamps@double vector;
        idoffset@uint32;        % lookup matrix
        
        ImgFileName@char vector;
        ImgSession@struct vector;
        
        ReconNode@util.reconNode vector;
        MSPNode@util.mspNode vector;
        
        structureTime@single = single(nan);
    end
    
    properties ( Access = private )
        structureLoaded@logical = false;
    end
    
    properties ( Dependent = true )
        RunString@cell;
        ZPositionString@cell;
        RepString@cell;
        WavelengthString@cell;
    end
    
    methods
        % *****************************************************************
        %% CONSTRUCTOR
        function obj = msotData(fileName,progress)
            % for array initialisation
            if nargin == 0
                return
            end
            if nargin >= 2 % If there's more than 1 input, the second 
                           %    should be the progress bool. 
                obj.progress = progress;
            end
            
            if isstring(fileName) % Make sure it's handled as a char.
                fileName = char(fileName);
            end
            
            % check if file exists
            if ~exist(fileName,'file')
                error('File %s does not exist - please specify .msot file',fileName);
            end
            
            % check file ending
            tok = regexpi(fileName,'(.+)\.msot(.+)?$','tokens');
            if isempty(tok)
                error('Please specify a .msot file (%s)', fileName);
            end
            binFile = [ tok{1}{1} '.bin' ];
            
            % find directory and study name
            sepind = strfind(fileName, filesep);
            if isempty(sepind)
                dirName = [ '.' filesep ];
                studyPath = ['..' filesep];
                studyName = '<unknown>';
            else
                dirName = fileName(1:sepind(end));
                try
                studyPath = fileName(1:sepind(end-1)-1);
                if numel(sepind) < 3
                    ss = 0;
                else
                    ss = sepind(end-2);
                end
                studyName = fileName(ss+1:sepind(end-1)-1);
                catch
                    studyPath = '<N/A>';
                    studyName = '<N/A>';
                end
            end
            
            
            % read XML
            if obj.progress
                wbar = waitbar(0,'Parsing MSOT File');
            end
            
            try
                patdoc = javaMethod( 'parse', 'com.itheramedical.msotbeans.DataModelMsotProjectDocument$Factory', java.io.File( fileName ) ) ; % JAVA
            catch ex
                error( ['Selected MSOT file could not be loaded. Check well-formedness.\n' ex.message]) ;
            end
            
            % store file information
            if obj.progress
                waitbar(0.5,wbar,'Reading Information...');
            end
            
            % annotate these. Maybe modify the parsing code (pull to function)
            obj.XMLFileName = fileName;
            obj.FileName = binFile;
            obj.RealPath = dirName;
            obj.StudyName = studyName;
            obj.StudyPath = studyPath;
            
            obj.DOM = patdoc.getDataModelMsotProject();
            obj.Name = char(obj.DOM.getScanNode.getName);
            obj.Comment = char(obj.DOM.getScanNode.getComment);
            try
                obj.DataModelVersion = char(obj.DOM.getDomNode.getElementsByTagName('DataModelVersion').item(0).getFirstChild.getNodeValue);
            catch
                obj.DataModelVersion = '2.0';
            end
            
            % Pull the wavelengths information from the DOM.
            localWavelengths = zeros(length(obj.DOM.getScanNode.getWavelengths.getWavelengthArray),1,'double');
            for i = 1:length(obj.DOM.getScanNode.getWavelengths.getWavelengthArray)
                localWavelengths(i) = double(obj.DOM.getScanNode.getWavelengths.getWavelengthArray(i-1).getIntValue);
            end
            obj.Wavelengths=localWavelengths;
            
            obj.FriendlyName = char(obj.DOM.getFriendlyName);
            obj.FolderName = char(obj.DOM.getFolderName);
            obj.Scientist = char(obj.DOM.getScientist);
            obj.CreationTimeTxt = char(obj.DOM.getCreationTime.toString);
            obj.CreationTime = datenum(obj.CreationTimeTxt,'yyyy-mm-ddTHH:MM:SS.FFF'); % Includes time zone.
            
            obj.is3D = strcmp(char(obj.DOM.getHARDWAREDESC.getTRANSDUCER),'msot3');
            obj.is2D = strcmp(char(obj.DOM.getHARDWAREDESC.getTRANSDUCER),'msot2');
            if (obj.is3D)
                obj.TransducerType = '3D';
            else
                obj.TransducerType = '2D';
            end
            
            
            % presence of US data. The metadata might not even have it.
            obj.USpresent = false;
            try obj.USpresent = logical(obj.DOM.getMEASUREMENTDESC.getULTRASOUNDPRESENT); catch, end
            
            
            
            % Imaging Session / View FileName
            xfn = strrep(obj.FileName,'.bin','.img');
            if exist(xfn,'file')
                obj.ImgFileName = xfn;
            end
            
            if obj.progress
                close(wbar);
            end
        end
        
        function tdom = loadImpulseResponse(obj)
            irpath = strrep((obj.XMLFileName),'msot','irf');
            if exist(irpath,'file')
                FID = fopen(irpath);
                tdom = fread(FID,'double')';
                fclose(FID);
                if obj.verbose, fprintf('Electrical Impulse Response loaded (%s)\n',irpath); end
            else
                par.impresp = [];
                warning('Cannot load impulse response (%s), file does not exist',irpath);
            end
        end
        
        
        % *****************************************************************
        %% Operators
        
        
        % == should compare acquisition dates (likely unique)
        function bool = eq(obj1,obj2)
            bool = (obj1.CreationTime == obj2.CreationTime);
        end
        
        % Allows treatment of the metadata as a vectorized object. Probably
        % should be rewritten and pulled out.
        function varargout = subsref(obj,s)
            switch s(1).type
                case '()'
                    varargout={obj.ScanFrames(s.subs{:})};
                case '.'
                    if length(s) == 1
                        
                        varargout = {obj.(s(1).subs)};
                        
                    elseif length(s) == 2 && strcmp(s(2).type,'()')
                        varargout = {obj.(s(1).subs)(s(2).subs{:})};
                    elseif length(s) == 3
                        try
                            varargout = {[obj.(s(1).subs)(s(2).subs{:}).(s(3).subs)]};
                        catch
                            varargout = {builtin('subsref',obj,s)};
                        end
                        
                    else
                        try
                            varargout = {builtin('subsref',obj,s)};
                        catch
                            %                       disp(obj)
                            %                       disp(s.type)
                            %                       disp(s.subs)
                        end
                        
                    end
                otherwise
                    if length(s)==1 && ismethod(obj,s.subs)
                        try
                            varargout={obj.(s.subs)};
                        catch
                            obj.(s.subs);
                            varargout={};
                        end
                    elseif s(1).subs
                        varargout={builtin('subsref',obj,s)};
                    end
                    
            end
        end
        
        function S = struct(obj)
            
            S = builtin('struct',obj);
            if obj.verbose
                S.DOM = char(obj.DOM.toString); % TODO: Modify this so that the system will only dump this in debug. 
            else
                obj.DOM = [];
            end
            S.ScanFrames = struct(S.ScanFrames);
            
            try
                S.ReconNode = struct(S.ReconNode);
            catch
                S.ReconNode = [];
            end
            
            try
                S.MSPNode = struct(S.MSPNode);
            catch
                S.MSPNode = [];
            end
        end
        
        
        % *****************************************************************
        %% Measurement Description
        function MD = get.MeasurementDesc(obj)
            if isempty(fieldnames(obj.MeasurementDesc))
                md = obj.DOM.getMEASUREMENTDESC;
                
                obj.MeasurementDesc.Averages = uint16(md.getAVERAGESPERPROJECTION); % number of shots averaged into one frame
                obj.MeasurementDesc.NumShots = uint16(md.getNUMBEROFFRAMES);       % number of shots per wavelength
                obj.MeasurementDesc.NumFrames = obj.MeasurementDesc.NumShots / obj.MeasurementDesc.Averages; % Number of frames stored in raw data
                obj.MeasurementDesc.RecordLength = uint16(md.getRECORDEDLENGTH);
                obj.MeasurementDesc.RepRate = single(md.getREPETITIONRATE.getDoubleValue);
                
                
                % wavelengths
                obj.MeasurementDesc.Wavelengths = zeros(length(md.getWAVELENGTHS.getWAVELENGTHArray),1,'double');
                for i = 1:length(md.getWAVELENGTHS.getWAVELENGTHArray)
                    obj.MeasurementDesc.Wavelengths(i) = md.getWAVELENGTHS.getWAVELENGTHArray(i-1).getIntValue;
                end
                obj.MeasurementDesc.Sequence = char(md.getSEQUENCE);
                
                % *** Water Absorption Coefficients
                try nwa = length(md.getWATERABSORPTIONCOEFFICIENTS.getWATERABSORPTIONCOEFFICIENTArray);
                catch, nwa = 0; end
                obj.MeasurementDesc.WaterAbsorptionCoeff = zeros(nwa,1,'single');
                for i = 1:nwa
                    obj.MeasurementDesc.WaterAbsorptionCoeff(i) = double(md.getWATERABSORPTIONCOEFFICIENTS.getWATERABSORPTIONCOEFFICIENTArray(i-1).getCoefficient);
                end
                
                % Path Length in Water (for Water correction)
                obj.MeasurementDesc.PathLengthInWater = nan;
                try obj.MeasurementDesc.PathLengthInWater = md.getPATHLENGTHINWATER.getDoubleValue; catch, end
                
                % *** Average Laser Energy Table (used if Average Energy Correction applied)
                try nfle = length(md.getAVERAGEENERGYTABLE.getAVRAGEENERGYArray);
                catch, nfle = 0; end
                %               obj.MeasurementDesc.FactoryLaserEnergyTable = zeros(nfle,1,'single');
                obj.MeasurementDesc.FactoryLaserEnergyTable = zeros(nfle,1,'double');
                for i = 1:nfle
                    obj.MeasurementDesc.FactoryLaserEnergyTable(i) = double(md.getAVERAGEENERGYTABLE.getAVRAGEENERGYArray(i-1).getDoubleValue);
                end
                
                % Type of Laser Energy Correction
                obj.MeasurementDesc.EnergyNormalisation = char(md.getLEnergyNormalization);
                
                sn = obj.DOM.getScanNode; % JAVA 
                obj.MeasurementDesc.SWVersion = '';
                obj.MeasurementDesc.TrimSOS = 0;
                obj.MeasurementDesc.CouplantCorr = false;
                try
                    obj.MeasurementDesc.SWVersion = char(sn.getSWVersion);
                    obj.MeasurementDesc.TrimSOS = double(sn.getTrimSpeedOfSound);
                    obj.MeasurementDesc.CouplantCorr = logical(sn.getCouplantCorrection);
                catch
                end
                
            end
            
            MD = obj.MeasurementDesc;
        end      % MeasurementDesc
        
        % ******************************************************************
        %% Hardware Description
        function HW = get.HWDesc(obj)
            if isempty(fieldnames(obj.HWDesc))
                hw = obj.DOM.getHARDWAREDESC; % JAVA
                
                obj.HWDesc.SamplingFrequency =  hw.getSAMPLINGFREQUENCY.getDoubleValue*1e6;
                obj.HWDesc.TransducerType = char(hw.getTRANSDUCER);
                
                % parametric description. Number of transducers is actually
                % derived from the angles.
                obj.HWDesc.StartAngle = []; % Angle of first transducer.
                obj.HWDesc.StepAngle = []; % Angle between transducer normals.
                obj.HWDesc.EndAngle = []; % Angle of the last transducer.
                obj.HWDesc.Radius = [];
                obj.HWDesc.RadiusZ = []; % Can probably get from the Dima/Burton paper.
                eq = hw.getFRAMEDESC.getEQUAL;
                if (~isempty(eq))
                    obj.HWDesc.StartAngle = double(eq.getSTART);
                    obj.HWDesc.StepAngle = double(eq.getSTEP);
                    obj.HWDesc.EndAngle = double(eq.getEND);
                    obj.HWDesc.NumDetectors = round((obj.HWDesc.EndAngle - obj.HWDesc.StartAngle) / obj.HWDesc.StepAngle)+1;
                    for j = 0:numel(eq.getCONSTANTArray)-1
                        cst = eq.getCONSTANTArray(j);
                        if cst.getAxisRef == 2
                            obj.HWDesc.Radius = cst.getDoubleValue;
                        elseif cst.getAxisRef == 3
                            obj.HWDesc.RadiusZ = cst.getDoubleValue; %Might be reasonable to hardcode this.
                        end
                    end
                end
                
                % sensor coordinate file, if available, gives coordinates
                % of every transducer.
                obj.HWDesc.ProjectionData = [ ];
                try % If it has projection data...
                    pa = hw.getFRAMEDESC.getPROJECTIONArray;
                    if (~isempty(pa)) %i.e. if it's >v3.6
                        obj.HWDesc.NumDetectors = numel(pa);
                        obj.HWDesc.ProjectionData = zeros(numel(pa),3);
                        for j = 1:numel(pa)
                            va = pa(j).getVALUEArray;
                            % 3 coordinates
                            obj.HWDesc.ProjectionData(j,:) = [va(1).doubleValue va(2).doubleValue va(3).doubleValue ];
                            obj.HWDesc.Radius = mean(sqrt(sum(obj.HWDesc.ProjectionData.^2,2)));
                        end
                        switch obj.DataModelVersion
                            case {'2.2'} % v3.8
                                % Get the minimum and maximum angles involved.
                                [thetaPol,rhoPol,zPol] = cart2pol(obj.HWDesc.ProjectionData(:,1),obj.HWDesc.ProjectionData(:,2),obj.HWDesc.ProjectionData(:,3));
                                obj.HWDesc.StartAngle = min(unwrap(thetaPol));
                                obj.HWDesc.EndAngle = max(unwrap(thetaPol));
                                
                                
                                obj.HWDesc.TransducerCoordinatesCartesian_X = obj.HWDesc.ProjectionData(:,1);
                                obj.HWDesc.TransducerCoordinatesCartesian_Y = obj.HWDesc.ProjectionData(:,2);
                                obj.HWDesc.TransducerCoordinatesCartesian_Z = obj.HWDesc.ProjectionData(:,3);
                                
                                obj.HWDesc.TransducerCoordinatesCylindrical_Theta = thetaPol;
                                obj.HWDesc.TransducerCoordinatesCylindrical_Rho   = rhoPol;
                                obj.HWDesc.TransducerCoordinatesCylindrical_Z     = zPol;
                            case {'2.3'} %v4.0
                                % Get the minimum and maximum angles involved.
                                [thetaPol,rhoPol,zPol] = cart2pol(obj.HWDesc.ProjectionData(:,1),obj.HWDesc.ProjectionData(:,3),obj.HWDesc.ProjectionData(:,2));
                                obj.HWDesc.StartAngle = min(unwrap(thetaPol));
                                obj.HWDesc.EndAngle = max(unwrap(thetaPol));
                                
                                
                                obj.HWDesc.TransducerCoordinatesCartesian_X = obj.HWDesc.ProjectionData(:,1);
                                obj.HWDesc.TransducerCoordinatesCartesian_Y = obj.HWDesc.ProjectionData(:,3);
                                obj.HWDesc.TransducerCoordinatesCartesian_Z = obj.HWDesc.ProjectionData(:,2);
                                
                                obj.HWDesc.TransducerCoordinatesCylindrical_Theta = thetaPol;
                                obj.HWDesc.TransducerCoordinatesCylindrical_Rho   = rhoPol;
                                obj.HWDesc.TransducerCoordinatesCylindrical_Z     = zPol;
                            otherwise
                        end
                    else % it has the other coordinate form.
                        
                        thetaPol = obj.HWDesc.StartAngle:obj.HWDesc.StepAngle:obj.HWDesc.EndAngle;
                        rhoPol   = obj.HWDesc.Radius .* ones(size(thetaPol));
                        zPol     = obj.HWDesc.RadiusZ .* ones(size(thetaPol));
                        
                        [xCart,yCart,zCart] = pol2cart(thetaPol,rhoPol,zPol);
                        
                        obj.HWDesc.TransducerCoordinatesCartesian_X = xCart;
                        obj.HWDesc.TransducerCoordinatesCartesian_Y = yCart;
                        obj.HWDesc.TransducerCoordinatesCartesian_Z = zCart;
                        
                        obj.HWDesc.TransducerCoordinatesCylindrical_Theta = thetaPol;
                        obj.HWDesc.TransducerCoordinatesCylindrical_Rho   = rhoPol;
                        obj.HWDesc.TransducerCoordinatesCylindrical_Z     = zPol;
                        
                    end
                catch
                end
                
                obj.HWDesc.SpeedOfSoundBase = 0;
                try
                    obj.HWDesc.SpeedOfSoundBase = double(hw.getSPEEDOFSOUNDBASE);
                catch
                end
                
                obj.HWDesc.LightSpotSize = 0;
                obj.HWDesc.AxialOffset = 0;
                obj.HWDesc.CouplantSOS = 0;
                obj.HWDesc.CouplantCOR = 0;
                obj.HWDesc.CouplantCurvature = 0;
                try
                    obj.HWDesc.LightSpotSize = double(hw.getLIGHTSPOTSIZE);
                    obj.HWDesc.AxialOffset = double(hw.getAXIALOFFSET);
                    obj.HWDesc.CouplantSOS = double(hw.getCOUPLANTSPEEDOFSOUND);
                    obj.HWDesc.CouplantCOR = double(hw.getCOUPLANTCENTEROFROTATION);
                    obj.HWDesc.CouplantCurvature = double(hw.getCOUPLANTCURVATURERADIUS);
                catch
                end
                
                
                sn = obj.DOM.getScanNode; % JAVA
                obj.HWDesc.DeviceSN = '';
                obj.HWDesc.TransducerSN = '';
                obj.HWDesc.DAQAddress = '';
                try
                    obj.HWDesc.DeviceSN = char(sn.getDeviceSN);
                    obj.HWDesc.TransducerSN = char(sn.getTransducerSN);
                    obj.HWDesc.DAQAddress = char(sn.getDAQMACaddress);
                catch
                end
            end
            
            
            obj.HWDesc.ImpulseResponseFunctionFile = strrep(obj.XMLFileName,'.msot','.irf');
            fID=fopen(obj.HWDesc.ImpulseResponseFunctionFile ,'r','l');
%             if fID ~= -1
                obj.HWDesc.ImpulseResponseFunction = fread(fID,Inf,'double');
%             end
            switch obj.DataModelVersion
                case {'2.2'} %v3.8
                    obj.HWDesc.ImpulseResponseFunction = fftshift(obj.HWDesc.ImpulseResponseFunction,1);
                case {'2.3'} %v3.9
                    obj.HWDesc.ImpulseResponseFunction = fftshift(obj.HWDesc.ImpulseResponseFunction,1);
                otherwise
            end
            fclose(fID);
            
            
            HW = obj.HWDesc;
        end      % hw
        
        
        
        
        % ******************************************************************
        
        %% OAMPreset
        % This is the state at time of acquisition, not necessarily the study
        % preset that is predefined
        function OA = get.OAMPreset(obj)
            if isempty(fieldnames(obj.OAMPreset))
                % ** Should go to OAM preset getter
                obj.OAMPreset = struct;
                obj.OAMPreset.Name = '';
                obj.OAMPreset.roi = [];
                try
                    oap = obj.DOM.getOAMPreset; % JAVA
                    obj.OAMPreset.Name = char(oap.getName);
                    %                     obj.OAMPreset.Ident = char(oap.getPresetIdentifier);
                    
                    rpa = obj.DOM.getOAMPreset.getReconPresets.getDataModelReconPreset.getSystemReconPresets; % JAVA
                    rp = rpa(1).getReconstructionPreset;
                    obj.OAMPreset.RoiHigh = double(rp.getRoiHigh);
                    obj.OAMPreset.RoiLow = double(rp.getRoiLow);
                    obj.OAMPreset.roi = double(rp.getRoiHigh) - double(rp.getRoiLow);
                    obj.OAMPreset.TimeRes = rp.getTimeRes;
                    obj.OAMPreset.Projections = rp.getProjections;
                    obj.OAMPreset.n = rp.getResolution;
                    
                    rpm = obj.DOM.getOAMPreset.getReconPresets.getDataModelReconPreset; % JAVA
                    obj.OAMPreset.FilterLow = rpm.getFilterLow;
                    obj.OAMPreset.FilterHigh = rpm.getFilterHigh;
                    obj.OAMPreset.DepthCorrection = logical(rpm.getDepthCorrection);
                    obj.OAMPreset.BackgroundAbsorption = double(rpm.getBackgroundAbsorption);
                    obj.OAMPreset.BackgroundOxygenation = rpm.getBackgroundOxygenation;
                    obj.OAMPreset.UserSoundTrim = rpm.getUserSoundTrim;
                    
                    mp = obj.DOM.getOAMPreset.getMspPresets.getDataModelMspPreset; % JAVA
                    obj.OAMPreset.MSP.DiscardNegatives = mp.getDiscardNegativeValues;
                    obj.OAMPreset.MSP.Method = char(mp.getMethod);
                    obj.OAMPreset.MSP.BGWavelength = mp.getBgWavelength.getIntValue;
                    obj.OAMPreset.MSP.Spectra = cell(mp.getUserSelectedSpectra.getStringArray);
                    
                    sett = obj.DOM.getOAMPreset.getImagingSettingsPreset; % JAVA
                    obj.OAMPreset.ViewSettings.UltrasoundMinimumScaling = double(sett.getUltrasoundMinimumScaling);
                    obj.OAMPreset.ViewSettings.UltrasoundMaximumScaling = double(sett.getUltrasoundMaximumScaling);
                    obj.OAMPreset.ViewSettings.BackgroundMinimumScaling = double(sett.getBackgroundMinimumScaling);
                    obj.OAMPreset.ViewSettings.BackgroundMaximumScaling = double(sett.getBackgroundMaximumScaling);
                    obj.OAMPreset.ViewSettings.ForegroundMinimumScaling = double(sett.getForegroundMinimumScaling);
                    obj.OAMPreset.ViewSettings.ForegroundMaximumScaling = double(sett.getForegroundMaximumScaling);
                    obj.OAMPreset.ViewSettings.Visible3DGridPlanesTypes = char(sett.getVisible3DGridPlanesTypes);
                    
                catch
                end
            end
            
            OA = obj.OAMPreset;
        end
        
        % ******************************************************************
        
        %% US Information
        function US = get.US(obj)
            if ~obj.USpresent, US = []; return; end
            if ~obj.structureLoaded, calculateStructure(obj); end
            
            if isempty(fieldnames(obj.US))
                % US Info (might not be present)
                obj.US.timestamps = single([]);
                obj.US.OAframes = uint32([]);
                obj.US.OArep = uint16([]);
                
                % resolution
                obj.US.n = 0;
                try obj.US.n = uint16(obj.DOM.getMEASUREMENTDESC.getULTRASOUNDRESOLUTION); catch, end
                
                for i = 1:numel(obj.ScanFrames)
                    uoff = obj.ScanFrames(i).USOffset;
                    if numel(obj.US.timestamps) < uoff+1
                        obj.US.timestamps(uoff+1) = obj.ScanFrames(i).RelTime;
                        obj.US.OAframes(uoff+1) = uint32(i);
                        obj.US.OArep(uoff+1) = uint16(obj.ScanFrames(i).Repetition);
                    end
                end
            end
            
            US = obj.US;
        end
        
        
        
        
        
        
        
        
        % ******************************************************************
        
        %% Structural Information (requires complete parse)
        function SF = get.ScanFrames(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            SF = obj.ScanFrames;
        end
        
        function SF = get.ZPositions(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            SF = obj.ZPositions;
        end
        
        function RN = get.RepNum(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            RN = obj.RepNum;
        end
        
        function RN = get.RunNum(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            RN = obj.RunNum;
        end
        
        function RN = get.ShotNum(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            RN = obj.ShotNum;
        end
        
        function AT = get.AverageTemperature(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            AT = obj.AverageTemperature;
        end
        
        function SS = get.ScanStructure(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            SS = obj.ScanStructure;
        end
        
        function SS = get.RelTime(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            SS = obj.RelTime;
        end
        
        function SS = get.LaserEnergy(obj)
            %            if ~obj.structureLoaded, calculateStructure(obj); end
            if isempty(obj.LaserEnergy)
                nFrames=obj.DOM.getScanNode.getScanFrames.sizeOfDataModelScanFrameArray; % JAVA
                tempLas=zeros(nFrames,1);
                sfa=obj.DOM.getScanNode.getScanFrames.getDataModelScanFrameArray; % JAVA
                for k=1:nFrames
                    tempLas(k)=sfa(k).getFrame.getLASERENERGY.getDoubleValue;
                end
                obj.LaserEnergy=double(tempLas);
                SS = double(tempLas);
                return;
            else
                SS= obj.LaserEnergy;
                return;
            end
            
        end
        
        function SS = get.DiodeReadout(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            SS = obj.DiodeReadout;
        end
        
        function SS = get.CorrectionFactor(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            SS = obj.CorrectionFactor;
        end
        
        function SS = get.Temperature(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            SS = obj.Temperature;
        end
        
        function SS = get.Timestamps(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            SS = squeeze(obj.RelTime(:,1,1,1,1));
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        % ******************************************************************
        
        %% Structure Parsing
        
        function calculateStructure(obj)
            tic %delete
            
            if obj.progress
                wbar = waitbar(0,'Reading Structural Information...');
            end
            
            obj.structureLoaded = true; % do this now to avoid recursion
            
            try
                sn = getScanNode(obj.DOM); %rename variable to scanNode
                fa = sn.getScanFrames.getDataModelScanFrameArray; %frameArray
                zpos = single([]); %zPositions. This doesn't even really do anything.
                ts0 = double(0); %timeStampZero
                nfr = length(fa); %nFrames
                locScanFrames(nfr) = util.scanFrame();  % initialise array at the end.
                
                for k = 1:nfr
                    if obj.progress && (mod(k,25) == 1), wbar = waitbar(k/nfr*0.75,wbar,'Reading Frame Information...'); end
                    frshell=fa(k);
                    fr = getFrame(frshell);
                    tempLocFrame=util.scanFrame();
                    
                    tempLocFrame.IDOffset = uint32(getIDOffset(frshell));
                    tempLocFrame.meta = obj; % As msotData is a handle class, this is fine.
                    tempLocFrame.FrameNum = uint32(k); %FrameNum is the index in the listing.
                    
                    try
                        wl1= intValue(getWAVELENGTH(fr));
                    catch
                        wl1 = nan;
                    end
                    
                    wl=double(wl1);
                    tempLocFrame.Wavelength = wl;
                    try % Need to figure out how to get X- and Y-pos.
                        zr= single(getDoubleValue(getZPOS(fr)));
                    catch
                        zr = single(0);
                    end
                    tempLocFrame.ZPosReal = zr;
                    
                    % Note: the correction factor is just the power reading
                    %                     try tempLocFrame.CorrectionFactor = single(doubleValue(getPOWER(fr))); catch, tempLocFrame.CorrectionFactor = single(nan); end
                    try tempLocFrame.CorrectionFactor = double(doubleValue(getPOWER(fr))); catch, tempLocFrame.CorrectionFactor = double(nan); end
                    try tempLocFrame.Temperature = single(doubleValue(getTEMPERATURE(fr))); catch,  tempLocFrame.Temperature = single(0); end
                    
                    ts = doubleValue(getTimestamp(fr)) * 1e-7;
                    if (k == 1), ts0 = ts; end
                    tempLocFrame.RelTime = double(ts - ts0);
                    tempLocFrame.Timestamp = double(ts);
                    
                    % Later Additions (might not be present)
                    %                     tempLocFrame.DiodeReadout = single(nan);
                    %                     tempLocFrame.LaserEnergy = single(nan);
                    %                     try tempLocFrame.LaserEnergy = single(getDoubleValue(getLASERENERGY(fr))); catch, end
                    %                     try tempLocFrame.DiodeReadout = single(getDoubleValue(getDIODEREADOUT(fr))); catch, end
                    tempLocFrame.DiodeReadout = double(nan);
                    tempLocFrame.LaserEnergy = double(nan);
                    try tempLocFrame.LaserEnergy = double(getDoubleValue(getLASERENERGY(fr))); catch, end
                    try tempLocFrame.DiodeReadout = double(getDoubleValue(getDIODEREADOUT(fr))); catch, end
                    
                    tempLocFrame.USOffset = uint32(nan);
                    
                    %%% for US machines
                    %                     try
                    %                         utest =fr.getUltraSoundFrameOffset;
                    %                         uoff = uint32(utest);
                    %                         locScanFrames(k).USOffset = uoff;
                    %                     catch
                    %                     end
                    
                    ru = uint16(getRUN(fr)); %run
                    if ru == 0, ru = uint16(1); end            % if no RUNs, then it will be zero
                    tempLocFrame.Run = ru;
                    tempLocFrame.Repetition = uint16(getREPETITION(fr));
                    
                    % Find ZPositions (with Rounding Problems) % Clarify
                    % this error. Issues with the pos. encoder.
                    [m, zp] = min(abs(zpos - zr));
                    if isempty(m) || m > 0.049 % Can make this more efficient.
                        zp = numel(zpos)+1;
                        zpos(zp) = zr;
                    end
                    tempLocFrame.ZNum = uint16(zp);
                    
                    
                    
                    % Wavelength index
                    tempLocFrame.WLNum = uint16(find(obj.Wavelengths == wl));
         
                    if isempty(obj.MeasurementDesc.WaterAbsorptionCoeff)
                        tempLocFrame.WaterAbsCoefficient = 1; % handle the situation that waterabscoefficient is not given
                    else
                        tempLocFrame.WaterAbsCoefficient = obj.MeasurementDesc.WaterAbsorptionCoeff(tempLocFrame.WLNum);
                    end
                    
                    locScanFrames(k)=tempLocFrame;
                end
                
                %%
                obj.RepNum = max([locScanFrames.Repetition]); % taken from the scan
                obj.RunNum = max([locScanFrames.Run]); % taken from the scan
                obj.AverageTemperature = mean([locScanFrames.Temperature]); %can probably extend this a little. CoV?
                
                % round and assign ZPositions similar to ViewMSOT
                zr = [locScanFrames.ZPosReal];
                zn = [locScanFrames.ZNum];
                for jj = 1:numel(zpos)% Introduces possible ambiguities.
                    zp = round(mean(zr(zn == jj))*10)/10;
                    zr(zn == jj) = zp;
                    [locScanFrames(zn == jj).ZPos] = deal(zp);
                    zpos(jj) = zp;
                end
                obj.ZPositions = zpos;
                
                % Number of Shots
                numshots = obj.MeasurementDesc.NumFrames;
                if (numshots == 0) % Legacy support for datasets without Number of Frames
                    numshots = round(numel(locScanFrames)./(obj.RunNum*numel(obj.ZPositions)*obj.RepNum*numel(obj.Wavelengths)));
                    obj.MeasurementDesc.NumFrames = numshots;
                end
                obj.ShotNum = numshots;
                
                % Initialise Structure Matrices - Make this more efficient.
                %
                locScanStructure = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),numshots,'single');
                obj.RelTime = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),obj.ShotNum,'double');
                %                 locLaserEnergy = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'single');
                locLaserEnergy = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'double');
                locDiodeReadout = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'single');
                %                 locCorrectionFactor = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'single');
                locCorrectionFactor = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'double');
                obj.Temperature = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'single');
                obj.idoffset = zeros(max([locScanFrames.IDOffset])+1,1,'uint32');
                
                oru = 0; ore = 0; ozp = 0; owl = 0; sh = 1;
                
                obj.ScanFrames=locScanFrames;
                for k = 1:nfr
                    if obj.progress && (mod(k,25) == 1), wbar = waitbar(0.75+k/nfr*0.25,wbar,'Reading Structural Information...'); end
                    %                     fr = obj.ScanFrames(k);
                    fr = locScanFrames(k);
                    
                    % reverse linking for recons and msps
                    obj.idoffset(fr.IDOffset+1) = k;  % map id offset to frame number
                    
                    %shotnumber.Can rewrite this for efficiency and clarity
                    if (oru ~= fr.Run) || (ore ~= fr.Repetition) || (ozp ~= fr.ZNum) || (owl ~= fr.WLNum)
                        sh = 1;
                        obj.RelTime(fr.Run,fr.ZNum,fr.Repetition,fr.WLNum) = fr.RelTime; % Only writes the data of the first shot- Redo!
                        obj.Temperature(fr.Run,fr.ZNum,fr.Repetition,fr.WLNum) = fr.Temperature;
                    else
                        sh = sh + 1;
                    end
                    fr.ShotNum = uint16(sh);
                    
                    runC=fr.Run; % Coordinates along each dimension into the metadata arrays.
                    zC=fr.ZNum;
                    repC=fr.Repetition;
                    wlC=fr.WLNum;
                    shotC=fr.ShotNum;
                    
                    locScanStructure(runC,zC,repC,wlC,shotC) = k;
                    locLaserEnergy(runC,zC,repC,wlC,shotC) = fr.LaserEnergy;
                    locDiodeReadout(runC,zC,repC,wlC,shotC) = fr.DiodeReadout;
                    locCorrectionFactor(runC,zC,repC,wlC,shotC) = fr.CorrectionFactor;
                    
                    % memorise for next time to see if there are more shots
                    % in that wavelength.
                    oru = fr.Run; ore = fr.Repetition; ozp = fr.ZNum; owl = fr.WLNum;
                end
                obj.ScanStructure = locScanStructure;
                obj.LaserEnergy = locLaserEnergy;
                obj.DiodeReadout = locDiodeReadout;
                obj.CorrectionFactor = locCorrectionFactor;
                
                obj.structureTime = single(toc); % Redo this to better report metadata.
                if obj.verbose, fprintf('Structure Information read in %.1fs\n',obj.structureTime); end
                if obj.progress, close(wbar); end
                
                
            catch ex
                if obj.progress, close(wbar); end
                obj.structureLoaded = false;
                rethrow(ex);
            end
        end
        
        
        
        % ******************************************************************
        %% View Information
        %        function IS = get.ImgSession(obj)
        %           if isempty(obj.ImgSession) && ~isempty(obj.ImgFileName)
        %               obj.ImgSession = [];%loadMSOTView(obj.ImgFileName,obj);
        %           end
        %
        %           IS = obj.ImgSession;
        %        end
        
        
        
        % ******************************************************************
        %% Recon Node Information
        function RN = get.ReconNode(obj)
            
            % If we've never looked at the recon node, populate it with the
            % headers of the recon nodes, but nothing else.
            if isempty(obj.ReconNode)
                try
                    rna = obj.DOM.getReconNodes.getDataModelNewReconstructionNodeArray;
                catch
                    RN = NaN;
                    return;
                end
                if obj.progress, wbar = waitbar(0,'Reading Recon Information...'); end
                
                % read all Recons
                for rn = 1:length(rna)
                    if obj.progress && mod(rn,25)==0
                        wbar = waitbar(rn/length(rna),wbar,'Reading Recon Information...');
                    end
                    
                    obj.ReconNode(rn) = util.reconNode(obj,rna(rn),uint16(rn));
                    
                end
                
                if isempty(obj.ReconNode)
                    obj.ReconNode(1)=util.reconNode;
                end
                
                if obj.progress, close(wbar); end
            end
            
            RN = obj.ReconNode;
        end
        
        
        % ******************************************************************
        %% MSP Node Information
        function MN = get.MSPNode(obj)
            if isempty(obj.MSPNode)
                try mna = obj.DOM.getMspNodes.getDataModelNewMspNodeArray;
                catch, MN = NaN; return; end
                
                if obj.progress, wbar = waitbar(0,'Reading MSP Information...'); end
                
                %
                for mn = 1:length(mna)
                    if obj.progress, waitbar(mn/length(mna),wbar) ; end
                    obj.MSPNode(mn) = util.mspNode(obj,mna(mn),uint16(mn));
                end
                
                if isempty(obj.MSPNode)
                    obj.MSPNode(1)=util.mspNode;
                end
                %                 obj.MSPNode='To be implemented';
                if obj.progress, close(wbar); end
            end
            
            MN = obj.MSPNode;
        end
        
        
        % ******************************************************************
        %% String Representations
        function c = get.RunString(obj)
            run = (1:obj.RunNum)';
            ts = squeeze(obj.RelTime(:,1,1,1,1));
            hour = floor(ts/60/60);
            min = floor(mod(ts,60*60)/60);
            sec = floor(mod(ts,60));
            
            if numel(run) ~= numel(ts), warning('Unsupported Dataset, number of timestamps != number of repetitions'); c = {}; end
            str = [num2str(hour,'%02i:') num2str(min,'%02i:') num2str(sec,'%02i') num2str(run,'x(%03i)')];
            c = mat2cell(str,ones(numel(ts),1));
            c = regexprep(c,'x',' ');
        end
        function c = get.ZPositionString(obj)
            if isempty(obj.ZPositions), c = {}; return; end
            c = mat2cell(num2str(reshape(obj.ZPositions,numel(obj.ZPositions),1),'%.1fmm'),ones(numel(obj.ZPositions),1));
        end
        function c = get.RepString(obj)
            reps = (1:obj.RepNum)';
            ts = squeeze(obj.RelTime(1,1,:,1,1));
            min = floor(ts/60);
            sec = floor(mod(ts,60));
            msec = floor(mod(ts,1)*100);
            
            if numel(reps) ~= numel(ts), warning('Unsupported Dataset, number of timestamps != number of repetitions'); c = {}; end
            str = [num2str(min,'%02i:') num2str(sec,'%02i.')  num2str(msec,'%02i') num2str(reps,'x(%03i)')];
            c = mat2cell(str,ones(numel(ts),1));
            c = regexprep(c,'x',' ');
        end
        function c = get.WavelengthString(obj)
            if isempty(obj.Wavelengths), c = {}; return; end
            c = mat2cell(num2str(obj.Wavelengths,'%inm'),ones(numel(obj.Wavelengths),1));
        end
        
        
        
        
        
        % ******************************************************************
        
        
        %% Temperature Plot
        function plotTemperature(obj)
            figure;
            plot([obj.ScanFrames(:).RelTime],[obj.ScanFrames(:).Temperature]);
            xlabel('Time (s)');
            ylabel(['Temperature (',char(176), 'C)']); axis tight;
        end
        
        %% Correction Factor Plot
        %% Diode Readout Plot
        %% Laser Energy Plot
        function plotEnergy(obj,plotDiode,typeString)
            
            if nargin == 1
                plotDiode = false;
                typeString='wavelength-time';
            end
            
            if nargin == 2
                typeString='wavelength-time';
            end
            
            switch typeString
                case 'time'
                    % Easiest case: Just plot all energies in sequential order.
                    times=[obj.ScanFrames(:).RelTime];
                    energy=[obj.ScanFrames(:).LaserEnergy];
                    
                    figure;
                    plot(times,energy);
                    ylabel('Laser Energy (mJ)');
                    xlabel('Time (s)'); axis tight;
                case 'wavelength'
                    % Average over all shots corresponding to each wavelength.
                    % Generates a laser energy spectrum.
                    [WL,~,C]=unique([obj.ScanFrames(:).Wavelength]);
                    WL_energy=accumarray(C,[obj.ScanFrames(:).LaserEnergy],[],@mean);
                    
                    if ~all(ismember(U,obj.Wavelengths))
                        warning('Mismatch in wavelengths');
                    end
                    
                    plot(WL,WL_energy);
                    
                case 'wavelength-run'
                    % Average over all shots within each run corresponding to
                    % each wavelength. Generates a trace for each WL.
                    
                case 'wavelength-rep'
                    % Average over all shots within each rep corresponding to
                    % each wavelength. Generates a trace for each WL.
                case 'wavelength-time'
                    % Plots a trace for each wavelength's energy over time.
                    WLenergyMat=repmat([obj.ScanFrames(:).LaserEnergy]',[1,numel(obj.Wavelengths)]);
                    
                    for k=1:numel(obj.Wavelengths)
                        WLenergyMat([obj.ScanFrames(:).Wavelength]~=obj.Wavelengths(k),k)=NaN;
                    end
                    figure;
                    plot([obj.ScanFrames(:).RelTime],WLenergyMat);
                    legend(num2str(obj.Wavelengths),'Location','northeastoutside')
                    ylabel('Laser Energy (mJ)');
                    xlabel('Time (s)'); axis tight;
            end
            
        end
        
    end       % methods
end      % classdef