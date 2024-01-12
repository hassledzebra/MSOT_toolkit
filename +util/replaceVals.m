
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
%       - feval:(functionHandle,arg1,arg2...)
%       - eval:() 
%       
%       - Cell expansion is also supported, by using braces instead of parentheses
%       for the bracketing operation. 
%           - e.g. 'eval:(@solver,frame:Data,{obj:CellArrayOfSolverParameters})'
%
%       - TODO: Calling 'replaceVals(obj)' with no second argument should do a
%       reflexive substitution on obj. In this case, it is the same as if you
%       had called 'replaceVals(obj,obj)'.
%       - TODO: Encoding the expected number of outputs to accommodate
%       things like evaluating high/low-pass filters. Since this only
%       occurs in the feval:() and eval:() cases, maybe change it to
%       something like feval:3:() or eval:1:() while preserving eval:() and
%       feval:() as shortcuts for eval:1:() and feval:1:() respectively.
%       

function completedSettings = replaceVals(obj,varargin)
    
    nVargs = numel(varargin);
    Q = evalin('caller','properties(obj)');
    
    switch nVargs
        case 0 % Reflexive action; let the obj use its own info for replacement.
            completedSettings = replaceVals(obj,obj);
        case 1
            vInput = varargin{1};
            
            if isstruct(vInput)
                completedSettings = replaceStruct(obj,vInput);
            elseif iscell(vInput) % Need to deploy and re-wrap the outputs.
                completedSettings = replaceCell(obj,vInput);
            elseif isstring(vInput) % Go a single layer deep.
            elseif ischar(vInput) % Treat as a string, but throw a warning.
            end
        otherwise % We assume this is basically just a long parameter list we want to fill in.
            % Consider N/V pairs
            % Consider mixed cells and not.
    end
    
    
    
end

function filledStruct = replaceStruct(obj,inputStruct)
    
    filledStruct = inputStruct;
    
    values = struct2cell(inputStruct);
    fields = fieldnames(inputStruct);
    
    subStructs = find(cellfun(@(x) isa(x,'struct'),values));
    nSubStructs = numel(subStructs);
    
    nonSubStructs = find(cellfun(@(x) ~isa(x,'struct'),values));
    
    % First, recurse on any structures.
    if nSubStructs > 0
        for k = subStructs
            filledStruct.(fields{k}) = replaceStruct(obj,values{k});
        end
    end
    
    % do a string comparison over values to find any settings that need to
    % inherit from the metadata.
    for nonStructIndex = 1 : numel(nonSubStructs)
        valueIndex = nonSubStructs(nonStructIndex);
        exposedString = values{valueIndex};
        if isstring(exposedString) 
            
            % Check to see if it's an eval call.
            newString = replaceString(obj,exposedString,false);
            
            filledStruct.(fields{valueIndex}) = newString;

        end
    end
    
end



function [filledString,varargout] = replaceString(obj,repString,doStringConversion)

    filledString = repString;
    % Check for eval calls.
    [a,b,c,d,e,f] = regexp(repString,'.{0}(?:eval:)\({1}(.*)\){1}.{0}','match','tokens');
    
    
    bEmpty = checkEmpty(b);
    
    if ~bEmpty % There's at least one eval call we need to handle. 
        evalContents = b{1};
        if iscell(evalContents)
            evalContents = evalContents{1};
        end
        parsedContents = replaceString(obj,evalContents,true);
        evalVal = eval(parsedContents);
        if doStringConversion
            filledString = string(evalVal);
        else
            filledString = evalVal;
        end
    else % Check for expansion.
        [aE,bE,cE,dE,eE,fE] = regexp(repString,'.{0}{((?:obj:|frame:|meta:).*)}.{0}','match','tokens');
        bEEmpty = checkEmpty(bE);
        if ~bEEmpty
            parsedContents = expandString(obj,bE{1});
            filledString = parsedContents;
        else
            [aS,bS,cS] = regexp(repString,'(?<source>obj|frame|meta):{1}(?<nav>(?:\(?([ \w:\.\(\)]*){1}\)+?)|([ \w:\.\(\)]*))','match','names');
            
            bSEmpty = checkEmpty(bS);
            if ~bSEmpty
                if strcmp(bS.source,'meta')
                    bS.source = 'obj';
                    bS.nav = strjoin(['Meta' bS.nav],'.');
                end
            [splitNameList,splitDotList] = split(bS(1).nav,'.');
                    
            fusedLists(1:2:(2*numel(splitNameList)-1)) = splitNameList;
            fusedLists(2:2:(2*numel(splitDotList)+1)) = splitDotList;
            
            

            refStruct = substruct('.',fusedLists{:});

            filledString = subsref(obj,refStruct);
            end
        end
    end
    
    function empCheck = checkEmpty(thing)
        if iscell(thing)
            if ~isempty(thing)
                empCheck = isempty(thing{1});
            else
                empCheck = isempty(thing);
            end
        else
            empCheck = isempty(thing);
        end 
    end
    
    
     
end



function [expandedString] = expandString(obj,expString)
   expandedString = replaceString(obj,expString,false); 
end



% Will need to do feval for many of these instead of eval, because eval can only
% handle simple things that can be inlined. Complex arg types like cells in
% objects are basically impossible to treat.
function [fevalOutput] = fevaluate(fevalString)
    
    
    
end



function [replacedCell] = replaceCell(obj,cellArray)

    nCellEl = numel(cellArray);
    replacedCell = cellArray;
    
    for k = 1:nCellEl
        entry = cellArray{k};
        if isstruct(entry)
            replacedCell{k} = replaceStruct(obj,entry);
        elseif isstring(entry) || ischar(entry)
            replacedCell{k} = replaceString(obj,entry,false);
        end
    end
    
end








%% Notes
    % The general parse should more or less break each evaluation string into a
    % tree of actions to undertake. The signature of the evaluation MUST BE 
    % one of:
    %   (<opname>:(<opstring>)) - for operators eval: and feval:
    %   <opname>:(<opstring>) - for operators eval: and feval:
    %   <opname>:<opstring> - for property access
    %   (<opname>:<opstring>) - for property access
    %   
    %   Commas should act to split opstrings, as long as they are within
    %   parentheses. 
    %   Consider:
    %       - command(meta:Prop1,
    %





%%









% Eval: has highest priority
% Expansion follows
% Everything else should be dictated by parentheses or commas.

% .{0}(?:eval:)\({1}(.*)\){1}.{0} - captures eval calls, returning the contents
% of their evalstring in a capture group. 

% .{0}{((?:eval:|obj:|frame:|meta:).*)}.{0} - grabs expansion calls. All other
% brackets are treated literally. NOTE: This might interact with eval in strange
% ways. 

% (?:eval:|obj:|frame:|meta:){1}[\w\(\)\.]+ - Grabs a lot of entries, but has
% trouble with eval:(etc) calls.

% 
% (obj|frame|meta):{1}\(?([ \w:\.\(\)]*){1}\)+? - Grabs the target and command and
% splits them into groups. Captures entire outside parenthetical group, but does
% not handle raw colon exposure.

% (obj|frame|meta):{1}(?:(?:\(?([ \w:\.\(\)]*){1}\)+?)|([ \w:\.\(\)]*)) - Grabs
% the target and command from either parentheticals or raw groups and puts them
% into one of two capture groups. 


% function repString = replaceStringWithMeta(inputString,obj)
%     
%     
%     repString = inputString;
%     
%     
%     if isstring(inputString) || ischar(inputString)
%         [isRep,imper,rem] = checkReplaceability(inputString);
%         
%         if isRep
%             switch imper
%                 case {'obj:'}
%                     [splitNameList,splitDotList] = split(rem,'.');
%                     
%                     fusedLists(1:2:(2*numel(splitNameList)-1)) = splitNameList;
%                     fusedLists(2:2:(2*numel(splitDotList)+1)) = splitDotList;
%                     
%                     refStruct = substruct('.',fusedLists{:});
%                     
%                     repString = subsref(obj,refStruct);
%                     
%                 case {'defaults:','def:'}
%                     
%                     repString = util.loadDefault(rem);
%                 case {'frame:'}
%                     
%                 case {'eval:'}
%                     assert(strcmp(rem(1),'('));
%                     assert(strcmp(rem(end),')'));
%                     reparsed = replaceForEval(rem,obj);
%                     repString = eval(reparsed);
%                     
%                 otherwise
%             end
%             
%             
%             
%         end
%     end
% end
% 
% 
% 
% 
% % Helper function for replaceValues - Performs the recursive
% % substitution of the inputStruct by traversing all fields looking for
% % 'meta:' as an initial string. This signifies where in the metadata
% % structure to look for the corresponding value.
% function filledStruct = replaceStructWithMeta(inputStruct,obj)
%     
%     filledStruct = inputStruct;
%     
%     values = struct2cell(inputStruct);
%     fields = fieldnames(inputStruct);
%     subStructs = cellfun(@(x) isa(x,'struct'),values);
%     nSubStructs = sum(subStructs);
%     
%     
%     % First, recurse on any structures.
%     if nSubStructs > 0
%         for k = 1 : nSubStructs
%             filledStruct.(fields{k}) = replaceStructWithMeta(values{k},obj);
%         end
%     end
%     
%     % do a string comparison over values to find any settings that need to
%     % inherit from the metadata.
%     for valueIndex = 1 : numel(values)
%         if ischar(values{valueIndex}) % If the value is a char, check to see if it's a meta coordinate.
%             
%             [isRep,imper,rem] = checkReplaceability(values{valueIndex});
%             
%             if isRep
%                 switch imper
%                     case {'obj:'}
%                         [splitNameList,splitDotList] = split(rem,'.');
%                         
%                         fusedLists(1:2:(2*numel(splitNameList)-1)) = splitNameList;
%                         fusedLists(2:2:(2*numel(splitDotList)+1)) = splitDotList;
%                         
%                         refStruct = substruct('.',fusedLists{:});
%                         
%                         filledStruct.(fields{valueIndex}) = subsref(obj,refStruct);
%                         
%                     case {'defaults:','def:'}
%                         
%                         filledStruct.(fields{valueIndex}) = util.loadDefault(rem);
%                     case {'frame:'}
%                         
%                     case {'eval:'}
%                         assert(strcmp(rem(1),'('));
%                         assert(strcmp(rem(end),')'));
%                         reparsed = replaceForEval(rem,obj);
%                         filledStruct.(fields{valueIndex}) = eval(reparsed);
%                         
%                     otherwise
%                 end
%                 
%                 
%                 
%             end
%         end
%     end
%     
% end
% 
% 
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
% 
% 
% 
% function [replacedString] = replaceForEval(evalRecipe,obj)
%     replacedString = evalRecipe;
%     [defMatch,defNoMatch] = regexp(evalRecipe,'(?:defaults|def):\w+','match','split');
%     
%     [objMatch,objNoMatch] = regexp(evalRecipe,'(?:obj):\w+','match','split');
%     for k = 1:numel(objMatch)
%         [isRep,imper,rem] = checkReplaceability(objMatch{k});
%         [splitNameList,splitDotList] = split(rem,'.');
%         
%         fusedLists(1:2:(2*numel(splitNameList)-1)) = splitNameList;
%         fusedLists(2:2:(2*numel(splitDotList)+1)) = splitDotList;
%         
%         refStruct = substruct('.',fusedLists{:});
%         
%         objMatch{k} = char(string(subsref(obj,refStruct)));
%         
%         
%     end
%     replacedString = strjoin(objNoMatch,objMatch);
%     %     [evalReplaceStart,evalReplaceEnd] = regexp(evalRecipe,'(?:defaults|def):\w+');
%     
%     [matc,tok] = regexp(evalRecipe,'obj:\w+');
% end
% 
% 
% 
% 
% function out = evalCommand(evalRecipe,obj)
%     
%     out = eval(replaceForEval(evalRecipe,obj));
%     
% end






