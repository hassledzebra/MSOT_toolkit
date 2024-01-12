% Set up cluster and parallel pool based on configuration files. 
%
%
% REQUIRED: The parallel configuration must:
%           - Be defined in a .mlsettings file in the config files folders

function configureParallel(parallelPath)
    
    switch nargin
        case 0
            parallelPath = util.loadDefault('PARALLEL_PROFILE');
        case 1
        otherwise
    end
    
    
    if isdeployed
        setmcruserdata('ParallelProfile',parallelPath);
    end
    
end

