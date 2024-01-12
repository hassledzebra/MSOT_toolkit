
function [parsedStruct]=JSON2struct(JSONfile)
    JSONString = fileread(JSONfile);
%Find any dots which are within a quoted string and are closed on the right
%with a colon; this implies that the containing string is a JSON key, which
%we need to preserve for later in case we're trying to using the key as a
%structure reference.
[matched,splitted] = regexp(JSONString,'"(?:[\w-]+\.[\w-]+)":{1}','match','split'); 
    dotSafe = strrep(matched,'.','_DOT_');
    safeString = strjoin(splitted,dotSafe);
    
    parsedStruct=jsondecode(safeString);
end