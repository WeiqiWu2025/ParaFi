function [bxSig] = survey_ParaFi_funcBackscatter(exSig,tagData,numSyms)
% Input
    % exSig: The incoming excitation signal;
    % tagData: The tag data;
    % numSyms: The number of symbols to modulate one bit
% Output
    % bxSig: The backscattered signal after modulated operation
global seqLenForEstChannel;
global numTags;

bxSig = exSig;
numTagData = length(tagData);
pulseLen = 80;
pfo = comm.PhaseFrequencyOffset('PhaseOffset',180);

for idx_1 = 1:numTagData
    if idx_1<=(seqLenForEstChannel*numTags) 
        bxSig(801+(idx_1-1)*pulseLen*numSyms:800+idx_1*pulseLen*numSyms) = tagData(idx_1).*bxSig(801+(idx_1-1)*pulseLen*numSyms:800+idx_1*pulseLen*numSyms);
    else
        if tagData(idx_1) == 0
            bxSig(801+(idx_1-1)*pulseLen*numSyms:800+idx_1*pulseLen*numSyms) = pfo(bxSig(801+(idx_1-1)*pulseLen*numSyms:800+idx_1*pulseLen*numSyms));
        end
    end
end

end    

