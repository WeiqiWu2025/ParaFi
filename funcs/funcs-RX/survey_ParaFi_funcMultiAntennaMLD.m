function [demdData] = survey_ParaFi_funcMultiAntennaMLD(possiValue,estH,pilots,ofdmDemodPilots,lenPayload,numTags)

global numMultiAntennas;

demdData = zeros(lenPayload,numTags);
possiValue(possiValue==0)=-1;

len = size(possiValue,1);
len_pilots = size(pilots,1);

for idx_1 = 1:lenPayload
    tmp_data = zeros(len_pilots*numMultiAntennas,numTags);
    for idx_2 = 1:numMultiAntennas
        for idx_3 = 1:len_pilots
            min = inf;
            indx = 0;
            y = ofdmDemodPilots{idx_2}(idx_3,idx_1);
            x = pilots(idx_3,idx_1); 
            for idx_4 = 1:len
                tmp_e = norm(y-(possiValue(idx_4,:)*estH{idx_2}(idx_3,:).').'.*x);
                if tmp_e < min
                    min = tmp_e;
                    indx = idx_4;
                end
            end
            tmp_f = possiValue(indx,:);
            tmp_f(tmp_f==-1)=0;
            tmp_data((idx_2-1)*numMultiAntennas+idx_3,:) =  tmp_f;
        end
    end
    for idx_5 = 1:numTags
        if length(find(tmp_data(:,idx_5)==1)) >= (numMultiAntennas*len_pilots/2)
            demdData(idx_1,idx_5) = 1;
        else
            demdData(idx_1,idx_5) = 0;
        end
    end
end

            
end

