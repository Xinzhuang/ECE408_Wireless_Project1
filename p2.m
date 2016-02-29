clear all;

hMod = comm.BPSKModulator;
hDemod = comm.BPSKDemodulator;
numBits = 1e6;


alpha_0 = .5;
theta_0 = 0;
%h0 = alpha_0 * exp(j * theta_0);
%h0 = .0001 .* raylrnd(numBits, 1) .* exp(j .* 2 * pi * rand(numBits, 1));

alpha_1 = .5;
theta_1 = 0;
%h1 = alpha_1 * exp(j * theta_1);
%h1 = .0001 .* raylrnd(numBits, 1) .* exp(j .* 2 * pi * rand(numBits, 1));

hAWGN = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');

SNRvec = 0:2:50;

fprintf('\n\nNo diversity\n');
berVec = zeros(length(SNRvec), 3);

for i = 1:length(SNRvec)
    
    data = randi([0 1], numBits, 1);
    s0 = step(hMod, data);
    SNR = SNRvec(i);

    hErr = comm.ErrorRate;

    iNoise = 1 * randn(numBits, 1);
    qNoise = 1 * randn(numBits, 1);
    r0 = real(s0) .* iNoise + j * imag(s0) .* qNoise;
    
    
    powerDB = 10*log10(var(s0));
    noiseVar = 10 .^ (0.1 * (powerDB - SNR));
    
    
    r0_n = step(hAWGN, r0, noiseVar);    
    s0_hat = real(r0_n).* conj(iNoise) + j * imag(r0_n) .* conj(qNoise) ;

    demodulated = step(hDemod, s0_hat);
    errorStats = step(hErr, data, demodulated);
    berVec(i,:) = errorStats;
    %    fprintf('SNR: %f\n', snr(r0, n0);
    fprintf(['Bit error rate = %5.2e\nNumber of errors = %d\nTotal ' ...
             'bits = %d\n\n'], errorStats)

end
figure
semilogy(SNRvec,berVec(:,1),'*')
grid on
xlabel('SNR (dB)')
ylabel('Bit Error Rate')

% fprintf('\n\n MRRC (1 Tx, 2 Rx)\n');
% for i = 1:length(SNRvec)
%     data = randi([0 1], numBits, 1);
%     s0 = step(hMod, data);
%     SNR = SNRvec(i);

%     hErr = comm.ErrorRate;

%     noisePower = sqrt(var(s0) * 10.^(-SNR/10));
%     n0 = randn(numBits, 1) .* noisePower;
%     snr(s0, n0)
%     n1 = randn(numBits, 1) .* noisePower;

%     r0 = (h0 .* s0) + n0;
%     r1 = (h1 .* s0) + n1;


%     s0_hat = conj(h0) .* r0  + conj(h1) .* r1;

%     demodulated = step(hDemod, s0_hat);
%     errorStats = step(hErr, data, demodulated);

%     fprintf(['Bit error rate = %5.2e\nNumber of errors = %d\nTotal ' ...
%              'bits = %d\n'], errorStats)

% end


