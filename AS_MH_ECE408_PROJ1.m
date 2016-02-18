%{
Alexander Serrano & Max Howald
ECE 408 - Wireless Communications
Prof. Keene
02/18/16

802.11a Standard
%}


%% SIMULATION PARAMETERS ( source: Mathworks) 

M = 4;                 % Modulation Order
k = log2(M);           % # of bits per symbol
numSC = 52;           % Number of OFDM subcarriers  (standard -> 52)
cpLen = 16;            % OFDM cyclic prefix length
maxBitErrors = 100;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted


%set convolutional encoder 
hConEnc = comm.ConvolutionalEncoder;

%set modulator and demodulator
hQPSKMod = comm.QPSKModulator('BitInput',true);
hQPSKDemod = comm.QPSKDemodulator('BitOutput',true);

hOFDMmod = comm.OFDMModulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);
hOFDMdemod = comm.OFDMDemodulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);