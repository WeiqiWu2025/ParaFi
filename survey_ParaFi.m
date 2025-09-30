clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993); % For repeatable results

%%%%%*** Waveform Configuration ***%%%%%
% Create a format configuration object for a 1-by-1 HT transmission
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 1; % 1 transmit antennas
cfgHT.NumSpaceTimeStreams = 1; % 1 space-time streams
cfgHT.PSDULength = 2000; % PSDU length in bytes % 64
cfgHT.MCS = 0; % 1 spatial streams, BPSK rate-1/2
cfgHT.ChannelCoding = 'BCC'; % BCC channel coding


fs = wlanSampleRate(cfgHT); % Get the baseband sampling rate
ofdmInfo = wlanHTOFDMInfo('HT-Data',cfgHT); % Get the OFDM info
ind = wlanFieldIndices(cfgHT); % Indices for accessing each field within the time-domain packet

%%%%%*** Simulation Parameters ***%%%%%
snr = 0:5:15;


global numTags;
numTags = 3;

global seqLenForEstChannel;
seqLenForEstChannel = 40;
preambleForEstChannel = survey_ParaFi_funcGeneratePreamble(seqLenForEstChannel,numTags);

maxNumPackets = 2000; % The maximum number of packets at an SNR point

S = numel(snr); % 返回数组snr中元素的个数
numBitErrs = zeros(S,numTags); % numBitErrs: The number of bit errors
berEst = zeros(S,numTags);
estSNR = zeros(maxNumPackets,numel(snr));

global numMultiAntennas;
numMultiAntennas = 1;

for i = 1:S
%     disp(['SNR: ',num2str(snr(i)),' dB...']);
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = i;
    RandStream.setGlobalStream(stream);
    
    % Loop to simulate multiple packets
    n = 1; % Index of packet transmitted
    while n<=maxNumPackets
        disp(['SNR: ',num2str(snr(i)),' dB -> ','n: ',num2str(n),'-th packet']);
        %%%%%*** TX side ***%%%%%
        % Generate a packet waveform
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgHT); % generate txWaveform
        tx = [tx; zeros(15,cfgHT.NumTransmitAntennas)];  % Add trailing zeros to allow for channel filter delay
        
        exSig = [];
        H_TX_Tags = [];
        %%%%%*** TX-Tags backscatter channel & AWGN
        for chan_tx_tag_idx1 = 1:numTags
            bxCoeffForTxTag_real = -1+(1+1)*rand(1,1);
            bxCoeffForTxTag_imag = -1+(1+1)*rand(1,1);
            bxCoeffForTxTag_real = bxCoeffForTxTag_real*0.1;
            bxCoeffForTxTag_imag = bxCoeffForTxTag_imag*0.1;
            bxCoeffForTxTag = bxCoeffForTxTag_real+1i*bxCoeffForTxTag_imag;
            tmp_exSig = tx.*bxCoeffForTxTag;
            exSig = [exSig,tmp_exSig];
            H_TX_Tags = [H_TX_Tags,bxCoeffForTxTag];
        end
        
        
        %%%%%*** Tags side ***%%%%%
        % compute the number of bits embeded in a packet
        temp = ceil((cfgHT.PSDULength*8+16+6)/26);
        numSymForPsdu = 0;
        numSymForTailPad = 0;
        if mod(temp,2) == 1
            numSymForPsdu = (numel(tx)-720-15-80-80-80)/80;
            numSymForTailPad = 2;
        else
            numSymForPsdu = (numel(tx)-720-15-80-80)/80;
            numSymForTailPad = 1;
        end
        numTagData = numSymForPsdu; % modulate one tag data per one symbol
        
        % Initial tags data
        tagData = zeros(numTagData,numTags);
        numPayload = numTagData-seqLenForEstChannel*numTags;
        actualPayloadBits = zeros(numPayload,numTags);
        for tag_idx1 = 1:numTags
            tmp_payload = randi([0,1],numPayload,1);
            actualPayloadBits(:,tag_idx1) = tmp_payload;
            numTail = numTagData-seqLenForEstChannel*numTags-length(actualPayloadBits(:,tag_idx1));
            tail = ones(numTail,1);
            tagData(:,tag_idx1) = [preambleForEstChannel(:,tag_idx1);actualPayloadBits(:,tag_idx1);tail];
        end
        
        % perform backscatter operation
        for tag_idx2 = 1:numTags
            bxSig{tag_idx2} = survey_ParaFi_funcBackscatter(exSig(:,tag_idx2),tagData(:,tag_idx2),1);
        end
        
        %%%%%** multi-antennas Tags-RX backscatter channel & AWGN **%%%%%
        H_Tags_RX = zeros(numMultiAntennas,numTags);
        for chan_tag_rx_idx1 = 1:numMultiAntennas
            for chan_tag_rx_idx2 = 1:numTags
                bxCoeffForTagRx_real = -1+(1+1)*rand(1,1);
                bxCoeffForTagRx_imag = -1+(1+1)*rand(1,1);
                bxCoeffForTagRx_real = bxCoeffForTagRx_real*0.01;
                bxCoeffForTagRx_imag = bxCoeffForTagRx_imag*0.01;
                bxCoeffForTagRx = bxCoeffForTagRx_real+1i*bxCoeffForTagRx_imag;
                rxSig{chan_tag_rx_idx1,chan_tag_rx_idx2} = bxSig{chan_tag_rx_idx2}.*bxCoeffForTagRx;
                H_Tags_RX(chan_tag_rx_idx1,chan_tag_rx_idx2) = bxCoeffForTagRx;
            end
            actualH{chan_tag_rx_idx1} = H_TX_Tags.*H_Tags_RX(chan_tag_rx_idx1,:);
        end

        %%%%%*** RX side ***%%%%%
        rx = complex(zeros(length(rxSig{1,1}),numMultiAntennas));
        for rx_idx1 = 1:numMultiAntennas
            for rx_idx2 = 1:numTags
                rx(:,rx_idx1) = rx(:,rx_idx1) + rxSig{rx_idx1,rx_idx2};
            end
            [rx(:,rx_idx1),~,~] = func_awgn(rx(:,rx_idx1),snr(i),'measured');
        end
        
        for rx_idx3 = 1:numMultiAntennas
            tmp_rx = rx(:,rx_idx3);
            ofdmDemodPilots{rx_idx3} = survey_ParaFi_funcReceiver(tmp_rx(ind.HTData(1):ind.HTData(2)),cfgHT,1);
        end
        z = 3;
        pilots = wlan.internal.htPilots(numSymForPsdu+1+numSymForTailPad,z,cfgHT.ChannelBandwidth,1);
        pilots = complex(pilots);
        
        % LS estimator
        for rx_idx4 = 1:numMultiAntennas
            tmp_ofdmDemodPilots = ofdmDemodPilots{rx_idx4};
            tmp_estH = ones(4,numTags);
            for rx_idx5 = 1:numTags
                A = tmp_ofdmDemodPilots(:,1+(1+(rx_idx5-1)*seqLenForEstChannel:rx_idx5*seqLenForEstChannel));
                B = pilots(:,1+(1+(rx_idx5-1)*seqLenForEstChannel:rx_idx5*seqLenForEstChannel));
                tmp_LL = size(A,1);
                for rx_idx6 = 1:tmp_LL
                    tmp_H_real = funcLSEstimator(B(rx_idx6,:)',real(A(rx_idx6,:))');
                    tmp_H_imag = funcLSEstimator(B(rx_idx6,:)',imag(A(rx_idx6,:))');
                    tmp_H = tmp_H_real + 1i*tmp_H_imag;
                    tmp_estH(rx_idx6,rx_idx5) = tmp_H;
                end
            end
            estH{rx_idx4} = tmp_estH;
        end
        
        
        %%%** Decoder **%%%
        for rx_idx7 = 1:numMultiAntennas
            payload_ofdmDemodPilots{rx_idx7} = ofdmDemodPilots{rx_idx7}(:,1+seqLenForEstChannel*numTags+1:end-numTail-numSymForTailPad);
        end
        payload_pilots = pilots(:,1+seqLenForEstChannel*numTags+1:end-numTail-numSymForTailPad);
        demodPayloadBits = survey_ParaFi_funcMultiAntennaDecoder(estH,payload_pilots,payload_ofdmDemodPilots);
          
        % calculate the number of bit errors
        for comm_idx1 = 1:numTags
            numBitErrs(i,comm_idx1) = numBitErrs(i,comm_idx1) + biterr(actualPayloadBits(:,comm_idx1),demodPayloadBits(:,comm_idx1));
        end
        n = n+1;
        
    end
    % calculate bit error rate
    for comm_idx2 = 1:numTags
        berEst(i,comm_idx2) = numBitErrs(i,comm_idx2)/(numPayload*maxNumPackets);
    end
    
end
aaa = 1;
