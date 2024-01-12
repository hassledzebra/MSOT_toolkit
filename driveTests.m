function [results] = driveTests(varargin)
    
    nVarg = numel(varargin);
    if nVarg == 1 % Assume that the submitted entry is 
        
    end
    
    
    %% Clean up the entire MATLAB setup.
    clear all; close all; clc;
    clear classes; %#ok<CLCLS>
    clear java; %#ok<CLJAVA>
    
    origPath = pwd;
    
    %% Move to the folder where this file is located, which should be at the
    % toolbox root.
    cd(fileparts(mfilename('fullpath')));
    
    %% Import the relevant classes for testing.
    import matlab.unittest.TestRunner
    import matlab.unittest.TestSuite
    import matlab.unittest.plugins.XMLPlugin
    import matlab.unittest.plugins.DiagnosticsRecordingPlugin
    
    warning('off','all');
    
    
    %% For function rewrite: load the suite for each of the classes that are passed in.
    % This should allow for e.g. having a big array of tests that we can more
    % programmatically call, in a similar way to the driver.m file executing other
    % pipeline steps.
    
    %% Test multiple things by concatenating.
    % suite1 = TestSuite.fromClass(?tests.testRecon.testingBackprojection_Classtests);
    % suite2 = TestSuite.fromClass(?tests.testRecon.testingModels_Classtests);
    % suite = [suite1 suite2];
    
    suite = [TestSuite.fromClass(?tests.testRecon.testingSimulatedReconstruction_Classtests),...
            TestSuite.fromClass(?tests.testRecon.testingBackprojection_Classtests),...
            TestSuite.fromClass(?tests.testRecon.testingModels_Classtests)];
            
    
    %
    runner = TestRunner.withTextOutput;
    
    artFolder = fullfile(pwd,'TESTING_ARTIFACTS');
    mkdir(artFolder);
    runner.ArtifactsRootFolder = artFolder;
    
    
    xmlFile = 'testResults.xml';
    xmlPlugin = XMLPlugin.producingJUnitFormat(xmlFile);
    diagnosticPlugin = tests.DiagnosticRecorderPlugin;
    
    runner.addPlugin(xmlPlugin)
    runner.addPlugin(DiagnosticsRecordingPlugin('IncludingPassingDiagnostics',true,'Verbosity',3));
    
    % Run the battery tests
%     tests.testRecon.testRecons_RandomBatteryData(artFolder);
%     tests.testRecon.testModels_RandomBatteryData(artFolder);
    
    
    % Run the suite
    results = runner.run(suite);
    
end
