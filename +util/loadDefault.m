function loadedDefault = loadDefault(valueString)
    % Loads defaults from possible defaults folders, with specific precedence
    % rules which depend on the environment in which the pipeline is being 
    % executed. Check to see if we're deployed. If so, we shouldn't use the
    % defaults enumeration file.
    %
    %
    % All defaults will be dominated by any command-line inputs the user
    % provided, but defaults can be found in several places, depending on the
    % context:
    %
    % When running the pipeline as a MATLAB invocation, the precedence rules are
    % fairly straightforward: default.m is the only place that defaults are
    % stored. This inherently prevents any duplicate settings. This is probably
    % behavior that we want to degenerate at some point to ensure that the
    % save/load behavior is the same across all invocations of the pipeline, no
    % matter how it's called. 
    %
    % When running the pipeline as a compiled application, the defaults are
    % instead stored in JSON files, within a folder ./json/defaults which is
    % local to one of three places, in increasing order of priority:
    %   - CTFRoot, where the compiled pipeline unpacks to.
    %   - the .exe root 
    %   - the folder where the .exe's driver script was invoked from
    override = 0;
    try
    if isdeployed
        % Check to see if the pwd (where the .exe is) has a json defaults folder
        exeJSONDefault = fullfile(pwd,'json','defaults');
        exeJSON = fullfile(pwd,'json');
        
        
        disp(pwd);
        
        
        if isfile(fullfile(exeJSON,[valueString,'.json']))
            defaultFile = fullfile(exeJSON,[valueString,'.json']);
            loadedDefault = util.JSON2struct(defaultFile);
             override = 1;
        elseif isdir(exeJSONDefault) %If there is a JSON folder...
            defaultFile = fullfile(exeJSONDefault,[valueString,'.json']);
            if isfile(defaultFile)
                loadedDefault = util.JSON2struct(defaultFile);
                override = 1;
            end
        end
        
        if ~override
            loadedDefault = util.JSON2struct(fullfile(defaultJSONPath,'defaults',[valueString,'.json']));
        end
        
    else % We're in a MATLAB session and can use the defaults structure. Look in the pwd first. 
        if isfile(fullfile('.',[valueString '.json']))
            loadedDefault = util.JSON2struct(fullfile('.',[valueString '.json']));
        else
            loadedDefault = unpack(util.defaults.(valueString));
        end
    end
    catch
        warning(['Default "',valueString,'" not found. Returning empty.']);
        loadedDefault = [];
    end
    
end