% Initial input parsing. This function sorts varargin into various buckets
% based on the input data type. Matched parameters all go into p.Results,
% while all other parameters should be gracefully put into a sensible
% variable in p.Unmatched
function [p] = parsePipelineInputs(msotFile,outputDirectory,varargin)
    
    
    p = inputParser;
    p.KeepUnmatched = true;
    p.StructExpand = false;
    
    %% TODO: Get defaults and valid sets by looking in the json folder for the appropriate file.
    jsonFolder = []; 
    
    %% Set up defaults, valid sets, and validation functions. 
%     default_reconModel = 'backproject2D';
%     default_reconSolver = 'lsqr';
%     default_fieldOfView = '25';
%     default_nPixels     = '200';
    
    %% Set up defaults, valid sets, and validation functions. 
    default_reconModel = unpack(util.defaults.defaultReconModel);
    default_reconSolver = unpack(util.defaults.defaultReconSolver);
    default_fieldOfView = unpack(util.defaults.defaultFieldOfView);
    default_nPixels     = unpack(util.defaults.defaultNumPixels);
    
    %% Create the actual inputParser. 
    % Required.
        addRequired(p,'MsotFile',@ischar); %char, string. Must be file that ends in .msot
        addRequired(p,'outputDirectory',@ischar); %char, string. Should be made if it doesn't exist.
    
    % Name/Value pairs.
        addOptional(p,'reconModel',default_reconModel); %char, string that matches valid set.
        addOptional(p,'reconSolver',default_reconSolver); %char, string that matches valid set.
        addOptional(p,'fieldOfView',default_fieldOfView); %numchar, numstring that matches valid set.
        addOptional(p,'nPixels',default_nPixels); %numchar, numstring
        
        
    %% Dig through varargin to do some initial formatting and cleanup. 
        
    % Get the number of optional entries. This is just to make sure that we
    % allow for a certain number of input values without corresponding
    % input names.
        nOptional = numel(struct(p).Optional);
    
    % Convert double-dashed entries to a single dash.
    for k = 1 : numel(varargin)
       if strncmp(varargin{k},'--',2) 
           varargin{k} = varargin{k}(2:end);
       end
    end
    
    
    % Separate out all inputs which are just structs themselves. 
    % TODO: Make sure that we can pass a struct as part of a name/value
    % pair.
    % a.s = 1; a.t = '1';    f(...,'-filter.settings',a,...) 
    isStruct = cellfun(@(x) isa(x,'struct'),varargin);
    structParams = varargin(isStruct);
    varargin = varargin(~isStruct);
    
    % Flag all parameters (dashed entries) and separate out the optional
    % arguments.
    hasLeadingDash = strncmp(varargin,'-',1);
    if any(hasLeadingDash)
        lastHeadless = find(hasLeadingDash,1,'first') - 1;
    else
        lastHeadless = numel(hasLeadingDash);
    end
    optionalParams = varargin(1:lastHeadless);
    varargin = varargin((lastHeadless+1):end);
    
    % Flag all parameters that have a 'do' at the beginning - these are
    % boolean flags which have an implicit 'true' after them.
    isBooleanParam = strncmp(varargin,'-do',3);
    tbooleanParams = varargin(isBooleanParam);
    booleanParams(1:2:(2*numel(tbooleanParams) - 1)) = cellfun(@(x) x(2:end),tbooleanParams,'UniformOutput',false);
    booleanParams(2:2:(2*numel(tbooleanParams))) = {true};
    varargin = varargin(~isBooleanParam);
    
    
    % Identify all remaining parameters that have a dot in their name - these are
    % structures coordinates that need to be formatted.
    isParamWithDot = cellfun(@(x) any(x == '.'),varargin);
    
    if ~isempty(isParamWithDot)
        isValueWithDot = circshift(isParamWithDot,[0 1]); isValueWithDot(1) = 0;

        structNames = varargin(isParamWithDot);
        structValues = varargin(isValueWithDot);
        outputStructures = makeStructuresFromParameters(structNames,structValues);
        deepStructures = cell(2,numel(outputStructures));
        for k = 1:numel(outputStructures)
           deepStructures{1,k} = outputStructures{k}.name;
           deepStructures{2,k} = outputStructures{k}.contents;        
        end
        remainingParams = varargin(~(isParamWithDot|isValueWithDot));
    
    else
        remainingParams = varargin(:);
        deepStructures = {};
        
    end
    
    
    
    
    % Stitch together all the inputs so that the parser can have a look at
    % them. 
    
    %% Parse everything. Keep unmatched for structure formatting. 
    parse(p,msotFile,outputDirectory,optionalParams{:},booleanParams{:},deepStructures{:},remainingParams{:},'unparsedStructure',structParams);
    
    
    
end
