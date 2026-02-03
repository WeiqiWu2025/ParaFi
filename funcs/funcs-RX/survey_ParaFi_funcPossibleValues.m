function [possiValue] = survey_ParaFi_funcPossibleValues(numConcurrentTags)

len = power(2,numConcurrentTags)-1;
possiValue = [];
for idx_1 = 0:len
    tmp = dec2bin(idx_1,numConcurrentTags)-48;
    possiValue = [possiValue;tmp];
end

end

