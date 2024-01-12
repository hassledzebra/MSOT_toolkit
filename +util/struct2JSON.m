function [fileflags,writeDestination]=struct2JSON(targetStruct,fileDestination,saveName)
    %% STRUCT2JSON serializes a struct as a JSON file.
    
    if nargin==3
        writeDestination=fullfile(fileDestination,[saveName '.json']);
    elseif nargin==2
        writeDestination=fullfile(fileDestination,[inputname(1),'.json']);
    elseif nargin==1 
        writeDestination=fullfile(pwd,[inputname(1) '.json']);
    end
    
    mkdir(fileparts(writeDestination));
    
    jsonString=jsonencode(targetStruct);
    fid=fopen(writeDestination,'w');
    fileflags=fprintf(fid,'%s',jsonString);
    fclose(fid);
   
end