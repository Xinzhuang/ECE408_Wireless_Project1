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
cpLen = 32;            % OFDM cyclic prefix length
maxBitErrors = 100;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted