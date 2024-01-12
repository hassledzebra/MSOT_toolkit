% Fill in any settings which have deferred their value to the meta object.
% Additionally, create the coordinate systems. 
% TODO: Check where we might want to abstract this coordinate system
% creation a little more. 
function completedSettings = settingsFromMeta(incompleteSettings,meta)
   
    %% Go through any references which ask for meta information. 
    values{1} = struct2cell(incompleteSettings);
    fields{1} = fieldnames(incompleteSettings);

    completedSettings = replaceStructWithMeta(incompleteSettings,meta);
    
    
    
end
% Helper function for settingsFromMeta - Performs the recursive
% substitution of the inputStruct by traversing all fields looking for
% 'meta:' as an initial string. This signifies where in the metadata
% structure to look for the corresponding value. 
function filledStruct = replaceStructWithMeta(inputStruct,meta)
    
    filledStruct = inputStruct;
    
    values = struct2cell(inputStruct);
    fields = fieldnames(inputStruct);
    subStructs = cellfun(@(x) isa(x,'struct'),values);
    nSubStructs = sum(subStructs);
    
    
    % First, recurse on any structures.
    if nSubStructs > 0
        for k = 1 : nSubStructs
            filledStruct.(fields{k}) = replaceStructWithMeta(values{k},meta);
        end
    end
    
    % do a string comparison over values to find any settings that need to
    % inherit from the metadata. 
    for valueIndex = 1 : numel(values) 
       if ischar(values{valueIndex}) % If the value is a char, check to see if it's a meta coordinate.
           stringTest = values{valueIndex};
           if strncmp(stringTest,'meta:',5)
               
               trimmedString = stringTest(6:end);
               [splitNameList,splitDotList] = split(trimmedString,'.');
               
               fusedLists(1:2:(2*numel(splitNameList)-1)) = splitNameList;
               fusedLists(2:2:(2*numel(splitDotList)+1)) = splitDotList;
               
               refStruct = substruct('.',fusedLists{:});
               
               filledStruct.(fields{valueIndex}) = subsref(meta,refStruct);
           end
       end
    end
    
end