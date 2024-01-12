function [safeVarString] = parseJSON(recipeFile,paramFile)
% This function should take as input a JSON-formatted string, or path to an
% accessible JSON file, and parse it out. Preprocessor directives should be
% supported.
%
% The function goes through a set of parsed states:
% -The raw parse, which just contains the nested structure of the JSON,
% -The preprocessor dependency graph, in which various values are set up in
% their evaluation order,
% -The parsed dependency graph, which should have all of the preprocessing values
% resolved.


recipeFile = fileread('/project/apps_database/MSOT/data/recipes/testRecipe2.json');
paramFile = fileread('/project/apps_database/MSOT/data/recipes/params.json');



[parsedObj] = parseObject(recipeFile);     % The parsed object that is returned will have the property that all 'normal' cell arrays have Nx1 entries, while only those that denote a decomposition of a structure or object will have an Nx2 shape, first column for names and second for contents.
try
    [parsedFuse] = parseObject(paramFile);
    
    replacedGlobals = replaceVars(parsedObj,parsedFuse);
catch
end



parsedHead = parsedObj;


% Get the dependency list for each item.
% For LHS:
%       - If it starts with a $, interpret it as a variable. This should
%       lead to it going into the variable Map.
%       - If it does not start with a $, and ends in #<num>, interpret it
%       as a reference.
%
% For RHS:
%       - If everything is enclosed in 'eval' or 'feval', tag it as such.
%       - If it contains $<name>, we know we have a dependency on $<name>.
%
%
% We should semantically end up with two graphs: One depicting all
% relationships and dependencies, even among preprocessor variables, and
% one depicting all relationships and dependencies among pipeline elements.
%
% By the end of preprocessing, there should be no $ variables anymore.

depList0 = getDependencyList(parsedHead,1);
adjGraph0 = makeAdjacencyGraph(parsedHead,depList0);

% First resolve all of the evaluations for the global variables so that
% those can be set.

G = digraph(adjGraph0,parsedHead(:,1));
figure;plot(G,'layout','layered');

% Get an evaluation order for the graph. Forcing a stable order means that
% variables written earlier in the recipe will be evaluated sooner.
[n,H] = toposort(G,'Order','stable');

% Begin evaluating objects. This proceeds in topological order, so that
% objects with fewer dependencies are evaluated earlier. This should chase
% down such things as CoordinateSystems and their Coordinate dependencies.
varMap = containers.Map;

evaledHead = parsedHead;

for s = 1:numel(n)
    
    evaledHead{n(s),2} = evalObject(parsedHead{n(s),2},varMap);
    
    
    varMap(parsedHead{n(s),1}) = evaledHead{n(s),2};
    
end

isGlobalPrec = cellfun(@(x) strncmp(x,'$',1),evaledHead(:,1));
% All remaining objects should now be objects which can actually be
% instantiated.
trimmedList = evaledHead(~isGlobalPrec,:);


% Get the dependency list of each of the objects. This will inform us as to
% whether or not we have an object's construction code in the codebase. 
depList = getDependencyList(trimmedList,1);
adjGraph = makeAdjacencyGraph(trimmedList ,depList);

G1 = digraph(adjGraph,trimmedList(:,1));
figure;plot(G1,'layout','layered');





% Get list of types to check for default replacement.
% typeList = cellfun(@(x) x{1,2},trimmedList(:,2),'uni',false);
% isDefault = strncmpi(typeList,'default:',8);
% 
% typeList = unique(getTypeList(trimmedList));

connectedComponents = conncomp(G1,'type','weak','OutputForm','cell');
connectedComponentsVec = conncomp(G1,'type','weak');


discComps = cellfun(@(x) subgraph(G1,x),connectedComponents,'uni', false);


ancestorList = arrayfun(@(x) isAncestor(discComps{2},1,x),1:numnodes(discComps{1}));
discInds = find(connectedComponentsVec==2);

restrictedGraph = subgraph(discComps{1},find(ancestorList));
restrictedObjList = trimmedList(discInds(ancestorList),:);

figure; plot(restrictedGraph);

[n1,H1] = toposort(restrictedGraph,'Order','stable');

initObj = instantiateObjectTree(restrictedObjList,n1,varMap);

for k = 1:size(initObj,1)
    lockedList(k) = isLocked(initObj{k,2});
end

initObj{1,2}(1,1);

end

% loaderType = trimmedList{5,2}{1,2};
% loaderVargs = trimmedList{5,2}(2:end,:).';
% argStruct = struct(loaderVargs{:});
% testLoader = feval(loaderType,argStruct);
% testFrame1 = testLoader(1);
% 
% % PrefilterInit = initObj{6,2};
% % 
% % test = PrefilterInit(testFrame1);
% 
% % ReconInit = initObj{4,2};
% % 
% % testRec = ReconInit(1);
% %%
% IteratorInit = initObj{1,2};
% 
% testIt = IteratorInit(1);
% release(IteratorInit);
% 
% 
% UnmixInit = initObj{2,2};
% 
% unmixIt = UnmixInit(1);
%%
% Determine if all of the types are included in the package. Should
% probably resolve all default:s before this point, so they can benefit
% from the expansion/evaluation.

% For each pipeline step, in order, create children. We may want to save
% the global vars so that they can resolve certain things (e.g. SoS tuning)
% at runtime. Invoke the constructor using the Type information, with the
% other arguments as N/V pairs from the rest of the object.

% If we're doing this top-down, then the ChildFilters will need to be given
% their instantiation values from their parents.
% 
% objectType = trimmedList{6,2}{1,2};
% objectVargs = trimmedList{6,2}(2:end,:).';
% argStruct = struct(objectVargs{:});
% testLoader = feval(objectType,argStruct);
% 
% objectType = trimmedList{8,2}{1,2};
% objectVargs = trimmedList{8,2}(2:end,:).';
% 
% for k = 1:numel(objectVargs)
%     if iscell(objectVargs{k})
%         objectVargs(k) = {objectVargs(k)};
%     end
% end
% argStruct = struct(objectVargs{:});
% testFilter = feval(objectType,argStruct);
% 
% testFrame1 = testLoader(1);
% testOutputFrame = testFilter(testFrame1);
% 
% 
% objectType = trimmedList{9,2}{1,2};
% objectVargs = trimmedList{9,2}(2:end,:).';
% argStruct = struct(objectVargs{:});
% testDeconvFilter = feval(objectType,argStruct);
% 
% deconvFrame = testDeconvFilter(testOutputFrame);
% 
% 
% objectType = trimmedList{14,2}{1,2};
% objectVargs = trimmedList{14,2}(2:end,:).';
% argStruct = struct(objectVargs{:});
% testInterpFilter = feval(objectType,argStruct);
% 
% interpFrame = testInterpFilter(deconvFrame);
% 
% 
% % Reconstruction
% objectType = trimmedList{15,2}{1,2};
% objectVargs = trimmedList{15,2}(2:end,:).';
% % for k = 1:numel(objectVargs)
% %     if iscell(objectVargs{k})
% %         objectVargs{k}={objectVargs(k)};
% %     end
% % end
% % argStruct = struct(objectVargs{:});
% testReconSystem = feval(objectType,objectVargs{:});
% 
% recFrame = testReconSystem(deconvFrame);
% recFrame.Data = reshape(recFrame.Data,[101 101]);
% figure; imagesc(recFrame.Data);
% 
% % State estimation.
% objectType = trimmedList{21,2}{1,2};
% objectVargs = trimmedList{21,2}(2:end,:).';
% 
% testMultispectralFilter = feval(objectType,objectVargs{:});
% msfiltFrame = testMultispectralFilter(recFrame);
% 
% %% Unmixing solution
% objectType = trimmedList{22,2}{1,2};
% objectVargs = trimmedList{22,2}(2:end,:).';
% 
% testUnmixSolver = feval(objectType,objectVargs{:});
% unmixsolveFrame = testUnmixSolver(msfiltFrame);
% 
% % Combo unmixing.
% objectType = trimmedList{19,2}{1,2};
% objectVargs = trimmedList{19,2}(2:end,:).';
% 
% 
% 
% initObj = initializeObjectTree(trimmedList{19,2},varMap);
% testUnmix = feval(objectType,objectVargs{:});
% 
% 
% 
% 
% unmixFrame = testUnmix(recFrame);
% 
% 
% 
% s = 1:10

%% Helper Functions


%%% replaceVars
% Replace values of variables in listSource with the values of those same
% variables in listOverwrite.
%
% listOverwrite should contain variables which do not begin with $, as we
% assume that these variables will generally come from some serialization
% process.

function replaced = replaceVars(listSource,listOverwrite)
    replaced = listSource;
    NOver = size(listOverwrite,1);
    
    for k = 1:NOver
        testName = strcat('$',listOverwrite{k,1}); % Append variable flag.
        matchedName = strcmp(listSource(:,1),testName);
        if any(matchedName) && ~isequaln(listOverwrite{k,2},listSource{find(matchedName),2}) % Check if any values match.
            if size(listOverwrite{k,2},1)==1
                fprintf('Overwriting %s : %s with %s\n',testName,string(listSource{find(matchedName),2}),string(listOverwrite{k,2}));
            else
                fprintf('Overwriting array contained in %s\n',testName);
            end
            replaced{matchedName,2} = listOverwrite{k,2};
        elseif ~any(matchedName)
            fprintf('No variable %s found in recipe file. Skipping.\n',testName);
        end
    end
end


%%% getTypeList
% Returns a list of input types for all entities contained within the cell
% array. Does not guarantee preservation of the relationships between
% parents and children.
function typeList = getTypeList(lists)
    typeList = {};
    
    if iscell(lists)&&(size(lists,2)==2)
        for k = 1:size(lists,1)
            if strcmp(lists{k,1},'Type') % We are in a terminal object which has its type declared.
                typeList = [typeList lists{k,2}];
            else % Dive into the children to see their types.
                typeList = [typeList getTypeList(lists{k,2})]; 
            end
        end
    end
end


%%% makeAdjacencyGraph
% Creates the adjacency graph matrix for the dependency list. 
function [adjGraph] = makeAdjacencyGraph(evaledHead,depList)
    
    adjGraph = zeros(numel(depList));
    
    
    for k = 1:size(depList,1)
        thisName = depList{k};
        
        if isempty(thisName) % Independent variable.
            continue;
        end
        
        for p = 1:numel(thisName)
            thisSource = thisName{p};
            if contains(thisSource,"#") && strncmp(thisSource,"$",1) % We are looking at a preprocessor-dereferenced object ref.
                sourceFlag = cellfun(@(x) strcmp(x,[extractAfter(thisSource{1},"$")]),evaledHead(:,1));
            else % We are looking at a static object ref.
                sourceFlag = cellfun(@(x) strcmp(x,[thisSource{1}]),evaledHead(:,1));
            end
            sourceInd = find(sourceFlag);
            
            adjGraph(sourceInd,k) = 1;
        end
        
    end
    
end


%%% getDependencyList
% Returns a list of all object references which are descendants of the
% given object. If rootFlag is set, then we construct an array of
% dependency lists, otherwise traverse the tree.
%
% depList is an Nx1 cell array of cell arrays if rootFlag == 1
% depList is an Nx1 cell array of strings if rootFlag ==0

function [depList] = getDependencyList(head,rootFlag)
    depList = {};
    
    if iscell(head) && size(head,2)==2 % We are looking at an object
        if rootFlag
            % Get the dependency list for each one of the objects in the
            % root.
            depList = cell(size(head,1),1);
            for k = 1:size(head,1)
                depList{k,1} = getDependencyList(head{k,2},0);
            end
        else
            for k = 1:size(head,1)
                % Concatenate all dependencies of the children into a
                % single array.
                depList = [depList getDependencyList(head{k,2},0)];
            end
        end
    else
        if iscell(head) % Cell array, traverse.
            for k = 1:size(head,1)
                depList = [depList getDependencyList(head{k},0)];
            end
        elseif isstring(head)||ischar(head)
            if contains(head,'$') % Preprocessor variable
                depList = regexp(head,'(\$\w[\w\d#:_]*)','tokens'); % Get the whole reference, including any object ref and preprocessor flag.
            elseif contains(head,'#')
                depList = regexp(head,'((?:\$?)[\w\d_]*#\d{1,3})','tokens'); 
            else
                depList = {};
            end
        end
    end
end


%%% evalObject
% Substitutes each object's dependencies with any values contained in the
% variable Map. Assumes traversal of the dependency list in a topologically
% sorted order, such that the first objects evaluated will have no
% dependencies, and the last objects evaluated will have the deepest
% dependencies.

function [evaledObj] = evalObject(parsedHead,varStruc)
    evaledObj = parsedHead;
    
    if iscell(parsedHead)
        for k = 1:numel(parsedHead) % Replace every element with the dereferenced value, even field names.
            evaledObj{k} = evalObject(parsedHead{k},varStruc);
        end
    elseif ischar(parsedHead)||isstring(parsedHead)
        strSelect = regexp(parsedHead,'(?<!f)eval:\((.*)\){1}','tokens'); %Match eval:() but not feval:()
        if ~isempty(strSelect)
            repString = regexprep(strSelect{1},'(\$[\w\d#:_]*)','varStruc(''$1'')'); % Replace any variables in the eval: statement with a corresponding key into varStruc
            try
                evaledObj = eval(repString{1}); % In the case of evaluations which return single outputs
            catch
                eval(repString{1}); % In the case of state-modifying variables like addpaths.
                evaledObj = "Done";
            end
        elseif contains(parsedHead,"$")
            % trim out literal quotes.
            if contains(parsedHead,'"')
                parsedHead = strip(parsedHead,'"');
            end
            if contains(parsedHead,'#')
                repString = regexprep(parsedHead,'\$([\w\d:_#]*)','varStruc(''$1'')');
            else
                repString = regexprep(parsedHead,'(\$[\w\d:_]*)','varStruc(''$1'')');
            end
            
            if iscell(repString)
                repString = repString{1};
            end
            
            try
                evaledObj = eval(repString);
            catch
                try
                    eval(repString);
                    evaledObj = "DONE";
                catch
                    evaledObj = "ERROR";
                end
            end
        end
    end
    
end




function [parsedObj] = parseObject(inputStr)
    
    trimStr = inputStr(2:(end-1));
    
    
    % https://www.mathworks.com/matlabcentral/answers/121920-how-do-i-match-nested-parenthesis-brackets-or-braces-with-dynamic-regular-expressions#answer_423254
    
    bracefun = @(s)sprintf('.{%d}',find(cumsum((char(s)==']' )-(char(s)=='['))>0,1,'first'));
    brackfun = @(s)sprintf('.{%d}',find(cumsum((char(s)=='}' )-(char(s)=='{'))>0,1,'first'));
    
    % Match braced values as arrays, bracketed values as objects,
    % double-quoted values as strings, literal 'true' or 'false' as
    % booleans, and numbers (including scientific notation) as numbers.
    regEx = '(?<arr>\[(??@bracefun($'')))|(?<obj>\{(??@brackfun($'')))|"(?<string>(?:[^\\"]|\\.)*?)"|(?<boolean>(true|false))|(?<number>\d*[.]?\d*[eE]?\d*)';
    [in,del,el,spl] = regexp(trimStr,regEx,'match','tokens','names','split');
    
    NItems = numel(del)/2; % In a properly formatted JSON, every entity should have a name.
    
    
    itnames = cell(NItems,1);
    itconts = cell(NItems,1);
    
    for k = 1:NItems
        nameInd = 2.*k - 1;
        valInd = 2.*k;
        itnames{k} = el(nameInd).string;
        
        if ~isempty(el(valInd).obj)
            itconts{k} = parseObject(el(valInd).obj); % We recurse on the object, including its braces
            
        elseif ~isempty(el(valInd).arr)
            itconts{k} = jsondecode(in{valInd}); % Use builtin parser for arrays.
            
        elseif ~isempty(el(valInd).string)
            if contains(el(valInd).string,'\"')
                el(valInd).string = strrep(el(valInd).string,'\"','"');
            end
            itconts{k} = el(valInd).string; % We have a string entry. TODO: Let this handle escaped strings.
            
        else
            itconts{k} = jsondecode(in{valInd}); % Everything else should be handled by the builtin JSON parser.
        end
        
        
    end
    
    parsedObj = [itnames itconts];
    
end



function [instantiatedProps] = instantiateObjectTree(objectTree,topoOrder,varMap)
    instantiatedProps = objectTree;
    
    for k = topoOrder
        thisObjectType = objectTree{k,2}{1,2};
        thisObjectVargs = objectTree{k,2}(2:end,:).';
        
        
        for q = 1:numel(thisObjectVargs)
            subEntry = thisObjectVargs{q};
            % Check to see if we have any references to objects in the
            % entry, or in the cell list.
            if isstring(subEntry) || ischar(subEntry)
                regmatch = ~isempty(regexp(subEntry,".*#\d"));
            elseif iscell(subEntry)
                for p = 1:numel(subEntry)
                    if isstring(subEntry{p}) || ischar(subEntry{p})
                        regmatch(p) = ~isempty(regexp(subEntry{p},".*#\d"));
                    else
                        regmatch(p) = 0; % Nothing else can match
                    end
                end
            end
            
            if any(regmatch)
                matchInds = find(regmatch);
                for p = matchInds
                    thisObjectVargs{q}{p} = varMap(subEntry{p}); % Index and substitute that object's properties.
                end
            end
        end
        isCellArrs = cellfun(@iscell,thisObjectVargs);
        
        
        if any(isCellArrs(:))
            
            thisObjectVargs(isCellArrs) = cellfun(@(x) {{x}},thisObjectVargs(isCellArrs));% Scalarize any cell arrays so they're treated correctly when structifying.
        end
        thisObjectArgStruct = struct(thisObjectVargs{:});
        
        thisEvaledObject = feval(thisObjectType,thisObjectArgStruct); % Create an object using the struct input arguments. 
        instantiatedProps{k,2} = thisEvaledObject;
        varMap(instantiatedProps{k,1})  = thisEvaledObject;
    end
    
end


function [isUpstream] = isAncestor(G,n1,n2)
    
    if n1 == n2
        isUpstream = true;
        return;
    end
    p = predecessors(G,n1);
    
    if isempty(p)
        isUpstream = false;
    elseif ismember(n2,p)
        isUpstream = true;
    else
        isUpstream = false;
        for k = 1:numel(p)
           isUpstream = isUpstream || isAncestor(G,p(k),n2); 
        end
    end
end