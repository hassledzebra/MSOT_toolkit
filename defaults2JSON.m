


defaultsInfo = ?util.defaults;

enumList = defaultsInfo.EnumerationMemberList.Name;

jsonPath = fullfile(fileparts(mfilename('fullpath')),'json','defaults');

for k = 1:numel(defaultsInfo.EnumerationMemberList)
    enumName = defaultsInfo.EnumerationMemberList(k).Name;
    
    valStruc = unpack(util.defaults.(enumName));
    

    util.struct2JSON(valStruc,jsonPath,enumName);

end