function [demodData] = survey_ParaFi_funcMultiAntennaDecoder(estH,pilots,ofdmDemodPilots)
global numTags;

lenPayload= size(ofdmDemodPilots{1},2);
% demodData = zeros(lenPayload,numTags);

possiValue = survey_ParaFi_funcPossibleValues(numTags);
demodData = survey_ParaFi_funcMultiAntennaMLD(possiValue,estH,pilots,ofdmDemodPilots,lenPayload,numTags);

end

