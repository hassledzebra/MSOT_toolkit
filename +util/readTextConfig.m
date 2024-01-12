clearvars; close all; clc;


testFile = '/home2/s418131/Work/Pipeline/Development/MSOT_Analysis_Pipeline/workflow/scripts/COMPILED/driver/for_testing/TEST/INPUT/params.txt';

inputString = fileread(testFile);







% Get each section
r = '(\w*):{([\s\S]+?)}';
[SectionList,tok] = regexp(inputString,r,'match','tokens') 
groupList = cellfun(@(x) x(1),tok);
varList = cellfun(@(x) x(2),tok);


% Split out variables inside each section. 
r='(\w+)[\s]*=[\s]*([^\#\n]+)(?=\s|\#|\n|\r|\f|$)';


[sections,toks] = regexp(varList,r,'match','tokens');

for k = 1:numel(groupList)
locToks = toks{k};

tokenVals = regexprep(cellfun(@(x) x{2},locToks,'uni',false).','[\n|\r]*','');
    tokenVals = convertNums(tokenVals);
    tokenVals = makeCellArraysFromBrackets(tokenVals);
    tokenVals = convertToCharsAndStrings(tokenVals);
tokenNames = cellfun(@(x) x{1},locToks,'uni',false).';


structy.(groupList{k}) = cell2struct(tokenVals,tokenNames,1); 

end




function convertedToNums = convertNums(inputVal)
    convertedToNums = inputVal;
    
    
    for k = 1:numel(inputVal)
       testConv = str2double(inputVal{k});
       if isnan(testConv)
           if strcmpi(inputVal{k},{'true','yes'})
                convertedToNums{k} = true;
           elseif strcmpi(inputVal{k},{'false','no'})
                convertedToNums{k} = false;
           else
           end
       else
           convertedToNums{k} = testConv;
       end
           
        
    end
    
    
    
end






function newCells = makeCellArraysFromBrackets(inputVal)
   newCells = inputVal;
   
   for k = 1:numel(inputVal)
       if inputVal{k}(1)=='[' && inputVal{k}(end)==']'
            [mat,tok] = regexp(inputVal{k},'[''|"](.+?)[''|"]','match','tokens');
            if ~isempty(tok)
                newCells{k} = cat(1,tok{:});
            end
       end
   end
    
    
end


function newStrings = convertToCharsAndStrings(inputVal)
   newStrings = inputVal;
   
   
   for k = 1:numel(newStrings)
      locVal = inputVal{k};
      
      if ischar(locVal) && numel(locVal)>=2
          if locVal(1)=='''' && locVal(end)==''''
              newStrings{k} = char(locVal(2:(end-1)));
          elseif locVal(1)=='"' && locVal(end)=='"'
              newStrings{k} = string(locVal(2:(end-1)));
          end
          newStrings{k} = strtrim(newStrings{k});
      elseif isstring(locVal)
          if startswith(locVal,'''') && endswith(locVal,'''')
              newStrings{k} = char(extractBetween(locVal,'''',''''));
          elseif startswith(locVal,'"') && endswith(locVal,'"')
              newStrings{k} = extractBetween(locVal,'"','"');
          end
          newStrings{k} = strtrim(newStrings{k});
      end
       
   end
    
    
    
end





