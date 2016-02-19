%{
Alexander Serrano & Max Howald
ECE 408 - Wireless Communications
Prof. Keene
02/18/16

802.11a Standard
%}


%% SIMULATION PARAMETERS ( source: Mathworks  ) 

%http://www.mathworks.com/help/comm/gs/qpsk-and-ofdm-with-matlab-system-objects-1.html
M = 4;                 % Modulation Order
k = log2(M);           % # of bits per symbol
numSC = 52;           % Number of OFDM subcarriers  (standard -> 52)
cpLen = 16;            % OFDM cyclic prefix length
maxBitErrors = 100;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted


%set convolutional encoder 
%hConEnc = comm.ConvolutionalEncoder;
% see link below for more info on convolutional encoding
% http://www.mathworks.com/help/comm/ref/comm.convolutionalencoder-class.html


%set modulator and demodulator
hQPSKMod = comm.QPSKModulator('BitInput',true);
hQPSKDemod = comm.QPSKDemodulator('BitOutput',true);

hOFDMmod = comm.OFDMModulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);
hOFDMdemod = comm.OFDMDemodulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);


hChan = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');

hError = comm.ErrorRate('ResetInputPort',true);

ofdmInfo = info(hOFDMmod); 

numDC = ofdmInfo.DataInputSize(1) ; 

frameSize = [k*numDC 1];

EbNoVec = (0:10)';
snrVec = EbNoVec + 10*log10(k) + 10*log10(numDC/numSC);

berVec = zeros(length(EbNoVec),3);
errorStats = zeros(1,3);

for m = 1:length(EbNoVec)
    snr = snrVec(m);

    while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits
        dataIn = randi([0,1],frameSize);              % Generate binary data
        qpskTx = step(hQPSKMod,dataIn);               % Apply QPSK modulation
        txSig = step(hOFDMmod,qpskTx);                % Apply OFDM modulation
        powerDB = 10*log10(var(txSig));               % Calculate Tx signal power
        noiseVar = 10.^(0.1*(powerDB-snr));           % Calculate the noise variance
        rxSig = step(hChan,txSig,noiseVar);           % Pass the signal through a noisy channel
        qpskRx = step(hOFDMdemod,rxSig);              % Apply OFDM demodulation
        dataOut = step(hQPSKDemod,qpskRx);            % Apply QPSK demodulation
        errorStats = step(hError,dataIn,dataOut,0);   % Collect error statistics
    end

    berVec(m,:) = errorStats;                         % Save BER data
    errorStats = step(hError,dataIn,dataOut,1);       % Reset the error rate calculator
end

berTheory = berawgn(EbNoVec,'psk',M,'nondiff');

figure
semilogy(EbNoVec,berVec(:,1),'*')
hold on
semilogy(EbNoVec,berTheory)
legend('Simulation','Theory','Location','Best')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
grid on
hold off
