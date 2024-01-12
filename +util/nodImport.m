function [nuStruct]=nodImport(filename)
    
    if isstring(filename)
        filename = char(filename);
    end
    docNode=xmlread(filename);
    
    % Begin parsing the child nodes.
    try
        theStruct = parseChildNodes(docNode);
    catch
        error('Unable to parse XML file %s.',filename);
    end
    
    %% Recurse again to clean up.
    
    nuStruct=parseStruct(theStruct);
    
    
end

% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
    % Recurse over node children.
    children = [];
    if theNode.hasChildNodes
        childNodes = theNode.getChildNodes;
        numChildNodes = childNodes.getLength;
        allocCell = cell(1, numChildNodes);
        
        children = struct(             ...
            'Name', allocCell, 'Attributes', allocCell,    ...
            'Data', allocCell, 'Children', allocCell);
        
        for count = 1:numChildNodes
            theChild = childNodes.item(count-1);
            children(count) = makeStructFromNode(theChild);
        end
    end
    
    if ~isempty(children)
        children=children(~cellfun(@isempty,{children.Name}));
    end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
    % Create structure of node info.
    
    nodeStruct = struct(                        ...
        'Name', char(theNode.getNodeName),       ...
        'Attributes', parseAttributes(theNode),  ...
        'Data', '',                              ...
        'Children', parseChildNodes(theNode));
    
    
    if any(strcmp(methods(theNode), 'getData'))
        nodeStruct.Data = char(theNode.getData);
    else
        nodeStruct.Data = '';
    end
    
    if numel(nodeStruct.Children)==1
       if strcmp(nodeStruct.Children(1).Name,'#text')
           nodeStruct.Data=nodeStruct.Children(1).Data;
           nodeStruct.Children=[];
       end
    end
    
    %%If it's a whitespace node, make it blank.
    if strcmp(nodeStruct.Name,'#text')&&...
            (isempty(nodeStruct.Data)||all(isspace(nodeStruct.Data)))
        nodeStruct.Name=[];
    end
    
    
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
    % Create attributes structure.
    
    attributes = [];
    if theNode.hasAttributes
        theAttributes = theNode.getAttributes;
        numAttributes = theAttributes.getLength;
        allocCell = cell(1, numAttributes);
        attributes = struct('Name', allocCell, 'Value', ...
            allocCell);
        
        for count = 1:numAttributes
            attrib = theAttributes.item(count-1);
            attributes(count).Name = char(attrib.getName);
            attributes(count).Value = char(attrib.getValue);
        end
    end
end

% ----- Local function PARSEATTRIBUTES -----
function doneStruct = parseStruct(theStruct)
    
    %We're mostly interested in recovering the information from the
    %children.
    
    childs=theStruct.Children;
    
    nChildren=numel(childs);
    
    for p=1:nChildren
       fieldnomen=childs(p).Name;
       fieldcontents=str2double(childs(p).Data);
       if ~isnan(fieldcontents)
        doneStruct.(fieldnomen)=fieldcontents; 
       else
        doneStruct.(fieldnomen)=childs(p).Data; 
       end
    end
    
end
