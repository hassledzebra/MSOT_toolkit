function [defaultJPath] = defaultJSONPath
    %defaultJSONPath Return the path to the JSON folder. For deployed
    %application, this will be the folder which contains the .exe 
    if isdeployed
    defaultJPath = fullfile(ctfroot,'json');
    else
    defaultJPath = fullfile(fileparts(mfilename('fullpath')),'json');
    end
end

