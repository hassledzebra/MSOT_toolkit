function [lintInfo] = mlint_check(targetPath)
%mlint_check Performs linting on the specified path
%    [lintInfo]=mlint_check(targetPath) returns a structure of linted
%    information for each .m file found on the target path. Each file on
%    the path is linted separately. Errors and warnings are returned as-is,
%    except for warnings about McCabe complexity from files with McCabe
%    complexity less than a given cutoff (default 10).
%    

    if nargin==0
        targetPath=pwd;
    end

    
    SEVERITY_CUTOFF=2;
    MCCABE_CUTOFF=10;
    
    curPath=pwd;
    cd(targetPath);
    
    % Generates the catalog of warnings and their categories. Used for
    % classification of severity.
    [warnings, categories] = mlintCatalog();
    
    % Get information about the m files on the target path. 
    lintInfo = getAllMFiles();
    
    % For every m file found, perform a lint check.
    max_severity = 0;
    NSevereErrors = 0;
    for k = numel(lintInfo):-1:1
        % Get the initial structure with relevant (non-Mathworks) warnings and
        % errors.
        testStruct = checkcode(fullfile(lintInfo(k).folder,lintInfo(k).name),'-cyc','-fullpath','-id');
        testStruct = trimMccabe(testStruct,MCCABE_CUTOFF);
        
        if ~isempty(testStruct) % Add the severity of each error.
            ids={testStruct.id};
            severity=cellfun(@(X) warnings.(X(:)).severity,ids,'UniformOutput',false);
            
            % If it's a McCabe warning, set severity to 2.
            severity(strcmp({testStruct(:).id},'CABE')) = {2}; 
            severity(strcmp({testStruct(:).id},'SCABE')) = {2}; 
            [testStruct.severity] = severity{:};
            
            isSevereError = [severity{:}]>SEVERITY_CUTOFF;
            if any(isSevereError)
                severeErrorInds = find(isSevereError);
                for q = 1:numel(severeErrorInds)
                    NSevereErrors = NSevereErrors + 1;
                    locSevStruc = testStruct(severeErrorInds(q));
                    locSevStruc.folder = lintInfo(k).folder;
                    locSevStruc.name = lintInfo(k).name;
                    severeErrors(NSevereErrors) = locSevStruc;
                end
            end
            max_severity = max(max([severity{:}]),max_severity); 
        end
        lintInfo(k).lintOutputs = testStruct;
    end
    
    if max_severity > SEVERITY_CUTOFF
        warning('There are serious errors in the code!');
        NSE = numel(severeErrors);
        for k = 1:NSE
           fprintf('Serious error in file: %s\n',severeErrors(k).name);
           fprintf('Severity %d\n',severeErrors(k).severity);
           fprintf('    Folder: %s\n',severeErrors(k).folder);
           fprintf('    Line: %d\n',severeErrors(k).line);
           fprintf('    Column Range: [%d,%d]\n',severeErrors(k).column(1),severeErrors(k).column(2));
           fprintf('    ID: %s\n',severeErrors(k).id);
           fprintf('    message: %s\n',severeErrors(k).message);
           
        end
        exit(2);
    end
    
    cd(curPath);
    exit(0);
end



%% Helper functions

% Adds an 'ispackage' field to the dir structure. Only works one level
% deep.


function dirStruct = getAllMFiles()
    % Recursively finds all m-files on the current path's tree. 
    dirStruct = dir('**/*.m')';
    
    % For every m-file found, locate it within the packaging structure.
    for k=numel(dirStruct):-1:1
        [~,parentFolder] = fileparts(dirStruct(k).folder);
        
        if parentFolder(1) == '+' % If the parent folder is a package folder...
            dirStruct(k).isInPackage = true;
            dirStruct(k).packagePath = getPackageChain(dirStruct(k).folder);
            dirStruct(k).isPrivate = false;
            
        elseif strcmp(parentFolder,'private') % Otherwise, if it's in a private folder...
            dirStruct(k).isInPackage = false;
            dirStruct(k).packagePath = '';
            dirStruct(k).isPrivate = true;
        else                                  % Otherwise it's on the path somewhere, but not directly inside a package.
            dirStruct(k).isInPackage = false;
            dirStruct(k).packagePath = '';
            dirStruct(k).isPrivate = false;
        end
    end
    
end


        function packageChain = getPackageChain(packagePath)
            packageChain = '';
            [parentPath,packageFolder] = fileparts(packagePath);
            [~,parentPackageFolder] = fileparts(parentPath);

            if parentPackageFolder(1) == '+'
                % recurse.
                packageChain = [getPackageChain(parentPath) '.'];
            end
            packageChain = [packageChain strrep(packageFolder,'+','')];

        end

%% Removes the 'McCabe Complexity' line if the McCabe complexity is <10.
function trimmedLint = trimMccabe(lintOutputs,cutoff)
    if isempty(lintOutputs)
        trimmedLint=lintOutputs;
        return;
    end
    
    mccabeWarnings=contains({lintOutputs.message},'McCabe')';
    isHighComplexity = cellfun(@(X) cell2double(regexp(X,'([1-9]+)\.','tokens')) > cutoff,{lintOutputs.message},'UniformOutput',false);
    
    % If it's not a mccabe warning OR if it is and has high complexity,
    % keep it.
    keepers = (~mccabeWarnings) | (mccabeWarnings & [isHighComplexity{:}]');
    
    trimmedLint=lintOutputs(keepers);
end

%% Converts cells to doubles, handling empty values as -Inf
function num = cell2double(val)
    if isempty(val)
        num=-Inf;
    else
        num = str2double(val{:});
    end
end



%% Parses any special flags that might be passed to the parser.
function specialFlags = getSpecialFlags(fileList,tags)
    
    nFiles=numel(fileList);
    
    for k=1:nFiles
        fileContents = fileread(filename);
        toDos = strfind(fileContents, 'TODO');
        fixMes = strfind(fileContents, 'FIXME');
    end
end



%% Gets a list of mLint warnings and categories.
% From
% https://stackoverflow.com/questions/35898444/find-category-of-matlab-mlint-warning-id.
function [warnings, categories] = mlintCatalog()
    % Get a list of all categories, mlint IDs, and severity rankings
    output = evalc('checkcode sum.m -allmsg');
    
    % Break each line into it's components
    lines = regexp(output, '\n', 'split').';
    pattern = '^\s*(?<id>[^\s]*)\s*(?<severity>\d*)\s*(?<message>.*?\s*$)';
    warnings = regexp(lines, pattern, 'names');
    warnings = cat(1, warnings{:});
    
    % Determine which ones are category names
    isCategory = cellfun(@isempty, {warnings.severity});
    categories = warnings(isCategory);
    
    % Fix up the category names
    pattern = '(^\s*=*\s*|\s*=*\s*$)';
    messages = {categories.message};
    categoryNames = cellfun(@(x)regexprep(x, pattern, ''), messages, 'uni', 0);
    [categories.message] = categoryNames{:};
    
    % Now pair each mlint ID with it's category
    comp = bsxfun(@gt, 1:numel(warnings), find(isCategory).');
    [category_id, ~] = find(diff(comp, [], 1) == -1);
    category_id(end+1:numel(warnings)) = numel(categories);
    
    % Assign a category field to each mlint ID
    [warnings.category] = categoryNames{category_id};
    
    category_id = num2cell(category_id);
    [warnings.category_id] = category_id{:};
    
    % Remove the categories from the warnings list
    warnings = warnings(~isCategory);
    
    % Convert warning severity to a number
    severity = num2cell(str2double({warnings.severity}));
    [warnings.severity] = severity{:};
    
    % Save just the categories
    categories = rmfield(categories, 'severity');
    
    % Convert array of structs to a struct where the MLINT ID is the field
    warnings = orderfields(cell2struct(num2cell(warnings), {warnings.id}));
end