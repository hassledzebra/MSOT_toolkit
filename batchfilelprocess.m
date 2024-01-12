studylist = [2:52];
scanlist = {[],[]};

msotFile = '/Volumes/Samsung_T5/ZH 2-DG 2-LMP/Scan_6/Scan_6.msot';
% msotFile = '/Volumes/Samsung_T5/breast msot data clinical/MSOT Data/Study_23/Scan_2/Scan_2.msot';
path = fullfile('/Volumes/Samsung_T5/ZH 2-DG 2-LMP/Scan_6');
% path = fullfile('/Volumes/Samsung_T5/breast msot data clinical/MSOT Data/Study_23/Scan_1');
modeltype = 'backproject2D';
loadMSOT(msotFile, path, modeltype)
%%
MSOT_FILENAME = '/Volumes/Samsung_T5/ZH 2-DG 2-LMP/Scan_6/Scan_6.msot';
metaData = util.msotData(MSOT_FILENAME);

calculateStructure(metaData);
loader = util.MSOTSignalLoader(MSOT_FILENAME);