function [stateOutput,stateLog] = slidingFunction(input,varargin)
    %SLIDINGFUNCTION Summary of this function goes here
    %   Detailed explanation goes here
    
    
    persistent slidingState
        if isempty(slidingState) % We assume that this is the first call to slidingFunction.
            msImageSize = varargin{1}; % should be a 3-element vector which describes the x,y, and c dimensions of the buffer.
            assert(all(msImageSize(1:2) == size(input.Data)));
            slidingState = zeros(msImageSize);
        end
        
        inputChannelIndex = input.Meta.WLNum;
        slidingState(:,:,inputChannelIndex) = input.Data;
        
        
        
        
    stateLog = '';
    stateOutput.Data = slidingState;
end

