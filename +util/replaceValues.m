% Fill in any settings which have deferred their value to some other object.
%
% The substitution strings are defined as basic grammars of a couple of
% different operators:
%       - obj: This reflexively indexes into the input object,
%       looking for properties and methods which can be used to replace the
%       entry.
%       - meta: This indexes into the Meta property of the input object.
%       Effectively a shorthand for obj:Meta.<stuff>
%       - defaults/default/def: This invokes the provided default, following the
%       default precedence rules. 
%       - frame: This invokes the frame's properties or methods.
%       
%       - Cell expansion is also supported, by using braces instead of parentheses
%       for the bracketing operation. 
%           - e.g. 'eval:(@solver,frame:Data,{obj:CellArrayOfSolverParameters})'


function completedSettings = replaceValues(obj,varargin)
    
    nVargs = numel(varargin);
    unpackFlag = false;
    if isempty(obj)
        completedSettings = obj;
    else
        % Check if it's a structure or something similar.

        if nVargs == 1
            completedSettings = replaceStructWithMeta(varargin{:},obj);
        elseif mod(nVargs,2)==0
            unpackFlag = true;
            completedSettings = cellfun(@(x) replaceStringWithMeta(x,obj),varargin,'uni',false);
        else
            error('Whoopsie!');
        end
    end
    
    
end


function repString = replaceStringWithMeta(inputString,obj)
    
    
    repString = inputString;
    
    
    if isstring(inputString) || ischar(inputString)
        [isRep,imper,rem] = checkReplaceability(inputString);
        
        if isRep
            switch imper
                case {'obj:'}
                    [splitNameList,splitDotList] = split(rem,'.');
                    
                    fusedLists(1:2:(2*numel(splitNameList)-1)) = splitNameList;
                    fusedLists(2:2:(2*numel(splitDotList)+1)) = splitDotList;
                    
                    refStruct = substruct('.',fusedLists{:});
                    
                    repString = subsref(obj,refStruct);
                    
                case {'defaults:','def:'}
                    
                    repString = util.loadDefault(rem);
                case {'frame:'}
                    
                case {'eval:'}
                    assert(strcmp(rem(1),'('));
                    assert(strcmp(rem(end),')'));
                    reparsed = replaceForEval(rem,obj);
                    repString = eval(reparsed);
                    
                otherwise
            end
            
            
            
        end
    end
end




% Helper function for replaceValues - Performs the recursive
% substitution of the inputStruct by traversing all fields looking for
% 'meta:' as an initial string. This signifies where in the metadata
% structure to look for the corresponding value.
function filledStruct = replaceStructWithMeta(inputStruct,obj)
    
    filledStruct = inputStruct;
    
    values = struct2cell(inputStruct);
    fields = fieldnames(inputStruct);
    subStructs = cellfun(@(x) isa(x,'struct'),values);
    nSubStructs = sum(subStructs);
    
    
    % First, recurse on any structures.
    if nSubStructs > 0
        for k = 1 : nSubStructs
            filledStruct.(fields{k}) = replaceStructWithMeta(values{k},obj);
        end
    end
    
    % do a string comparison over values to find any settings that need to
    % inherit from the metadata.
    for valueIndex = 1 : numel(values)
        if ischar(values{valueIndex}) % If the value is a char, check to see if it's a meta coordinate.
            
            [isRep,imper,rem] = checkReplaceability(values{valueIndex});
            
            if isRep
                switch imper
                    case {'obj:'}
                        [splitNameList,splitDotList] = split(rem,'.');
                        
                        fusedLists(1:2:(2*numel(splitNameList)-1)) = splitNameList;
                        fusedLists(2:2:(2*numel(splitDotList)+1)) = splitDotList;
                        
                        refStruct = substruct('.',fusedLists{:});
                        
                        filledStruct.(fields{valueIndex}) = subsref(obj,refStruct);
                        
                    case {'defaults:','def:'}
                        
                        filledStruct.(fields{valueIndex}) = util.loadDefault(rem);
                    case {'frame:'}
                        
                    case {'eval:'}
                        assert(strcmp(rem(1),'('));
                        assert(strcmp(rem(end),')'));
                        reparsed = replaceForEval(rem,obj);
                        filledStruct.(fields{valueIndex}) = eval(reparsed);
                        
                    otherwise
                end
                
                
                
            end
        end
    end
    
end


function [isReplaceable,imperativeCapture,remainderCapture] = checkReplaceability(recipeString)
    
    isReplaceable = true;
    imperativeCapture = '';
    remainderCapture = '';
    
    
    firstColon = strfind(recipeString,':');
    
    imperativeCapture = recipeString(1:firstColon);
    remainderCapture = recipeString((firstColon+1):end);
    
    switch imperativeCapture
        case {'obj:'}
            
        case {'defaults:','def:'}
            
        case {'frame:'}
            
        case {'eval:'}
            
        otherwise
            isReplaceable = false;
    end
    
end



function [replacedString] = replaceForEval(evalRecipe,obj)
    replacedString = evalRecipe;
    [defMatch,defNoMatch] = regexp(evalRecipe,'(?:defaults|def):\w+','match','split');
    
    [objMatch,objNoMatch] = regexp(evalRecipe,'(?:obj):\w+','match','split');
    for k = 1:numel(objMatch)
        [isRep,imper,rem] = checkReplaceability(objMatch{k});
        [splitNameList,splitDotList] = split(rem,'.');
        
        fusedLists(1:2:(2*numel(splitNameList)-1)) = splitNameList;
        fusedLists(2:2:(2*numel(splitDotList)+1)) = splitDotList;
        
        refStruct = substruct('.',fusedLists{:});
        
        objMatch{k} = char(string(subsref(obj,refStruct)));
        
        
    end
    replacedString = strjoin(objNoMatch,objMatch);
    %     [evalReplaceStart,evalReplaceEnd] = regexp(evalRecipe,'(?:defaults|def):\w+');
    
    [matc,tok] = regexp(evalRecipe,'obj:\w+');
end




function out = evalCommand(evalRecipe,obj)
    
    out = eval(replaceForEval(evalRecipe,obj));
    
end






