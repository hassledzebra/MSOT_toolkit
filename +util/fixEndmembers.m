function fixedEndmembers = fixEndmembers(endmembers)
    
    listClass = class(endmembers);
    
%     endmembers
    
    switch listClass
        
        case {'char'}
            if strcmp(endmembers(1),'[') && strcmp(endmembers(end),']')
                endmembers = endmembers(2:(end-1));
            elseif strcmp(endmembers(1:2),'{[') && strcmp(endmembers((end-1):end),']}')
                endmembers = endmembers(3:(end-2));
            end
            
            splitLines = strsplit(endmembers,'(, )(?=[^\]]*(?:\[|$))','DelimiterType','RegularExpression').';
            
%             EndmemberGlossary = util.defaults.EndmemberIndex;
            endMemRef = util.loadDefault('EndmemberIndex');%EndmemberGlossary.contents{:};
            absorberNames = fieldnames(endMemRef);
            absorberAliases = struct2cell(endMemRef);
            
            for k = 1:numel(splitLines)
                candidateName = strrep(splitLines{k},"'","");
                try
                    [r,c] = ind2sub(size(absorberAliases),find(cellfun(@(x) any(strcmp(candidateName,x)),absorberAliases)));
                    replacementName = absorberNames{r,1};
%                     disp('1')
                    fixedEndmembers{k} = replacementName;
                catch
                    warning(['Unable to find equivalent endmember of ',candidateName,' in database.']);
                    fixedEndmembers{k} = splitLines{k};
                end
            end
        case {'cell'}
            endMemRef = util.loadDefault('EndmemberIndex');
            absorberNames = fieldnames(endMemRef);
            absorberAliases = struct2cell(endMemRef);
            for k = 1:numel(endmembers)
                candidateName = strrep(endmembers{k},"'","");
                try
%                     disp('2');
                    [r,c] = ind2sub(size(absorberAliases),find(cellfun(@(x) any(strcmp(candidateName,x)),absorberAliases)));
                    replacementName = absorberNames{r,1};
                    fixedEndmembers{k} = replacementName;
                catch
                    warning(['Unable to find equivalent endmember of ',candidateName,' in database.']);
                    fixedEndmembers{k} = endmembers{k};
                end
            end
            
            
            
        otherwise
            disp("Unhandled endmember name inputs");
    end
    
    
    x = 0;
end