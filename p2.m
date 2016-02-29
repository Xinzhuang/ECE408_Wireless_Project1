clear all;

hMod = comm.BPSKModulator;
hDemod = comm.BPSKDemodulator;
numBits = 5e6;

SNRvec = 0:2.5:50;
berVec = zeros(length(SNRvec), 3);
fprintf('\n\nNo diversity\n');


for i = 1:length(SNRvec)
    
    data = randi([0 1], numBits, 1);
    s0 = step(hMod, data);
    SNR = SNRvec(i);

    hErr = comm.ErrorRate;

    iNoise = 1 * randn(numBits, 1);
    qNoise = 1 * randn(numBits, 1);
    r0 = real(s0) .* iNoise + j * imag(s0) .* qNoise;
    
       
    r0_n = awgn(r0, SNR * 2);
    s0_hat = real(r0_n).* iNoise - j * imag(r0_n) .* conj(qNoise) ;

    demodulated = step(hDemod, s0_hat);
    errorStats = step(hErr, data, demodulated);
    berVec(i,:) = errorStats;
    %    fprintf('SNR: %f\n', snr(r0, n0);
    fprintf(['Bit error rate = %5.2e\nNumber of errors = %d\nTotal ' ...
             'bits = %d\n\n'], errorStats)

end
figure
semilogy(SNRvec,berVec(:,1),'-o', 'LineWidth',2)
grid on
xlabel('SNR (dB)')
ylabel('Bit Error Rate')
hold on

SNRvec = 0:2.5:25;
berVec = zeros(length(SNRvec), 3);
fprintf('\n\n 1 Tx 2 Rx\n');

for i = 1:length(SNRvec)
    
    data = randi([0 1], numBits, 1);
    s0 = step(hMod, data);
    SNR = SNRvec(i);

    hErr = comm.ErrorRate;

    iNoise0 = 1 * randn(numBits, 1);
    qNoise0 = 1 * randn(numBits, 1);
    r0 = real(s0) .* iNoise0 + j * imag(s0) .* qNoise0;
    
    iNoise1 = 1 * randn(numBits, 1);
    qNoise1 = 1 * randn(numBits, 1);

    r1 = real(s0) .* iNoise1 + j * imag(s0) .* qNoise1;
    
       
    r0_n = awgn(r0, SNR * 2);
    r1_n = awgn(r1, SNR * 2);
    
    r0_hat = real(r0_n).* iNoise0 - j * imag(r0_n) .* qNoise0;

    r1_hat = real(r1_n).* iNoise1  - j * imag(r1_n) .* qNoise1;

    s0_hat = r0_hat + r1_hat;
    demodulated = step(hDemod, s0_hat);
    errorStats = step(hErr, data, demodulated);
    berVec(i,:) = errorStats;
    %    fprintf('SNR: %f\n', snr(r0, n0);
    fprintf(['Bit error rate = %5.2e\nNumber of errors = %d\nTotal ' ...
             'bits = %d\n\n'], errorStats)

end
semilogy(SNRvec,berVec(:,1),'-rv', 'LineWidth',2)


SNRvec = 0:2.5:30;
berVec = zeros(length(SNRvec), 3);

fprintf('\n\n 2 Tx 1 Rx\n');
for i = 1:length(SNRvec)
    
    data = randi([0 1], numBits, 1);
    s0 = step(hMod, data);
    SNR = SNRvec(i);

    hErr = comm.ErrorRate;

    antenna_0 = s0;
    antenna_0(2:2:end) = -1 .* conj(s0(2:2:end));

    antenna_1 = complex(zeros(numBits, 1));
    antenna_1(1:2:end) = s0(2:2:end);
    antenna_1(2:2:end) = conj(s0(1:2:end));
        
    iNoise0 = 1 * randn(numBits/2, 1);
    qNoise0 = 1 * randn(numBits/2, 1);
    
    iInterleaved0 = zeros(numBits, 1);
    iInterleaved0(1:2:end) = iNoise0;
    iInterleaved0(2:2:end) = iNoise0;
    
    qInterleaved0 = zeros(numBits, 1);
    qInterleaved0(1:2:end) = qNoise0;
    qInterleaved0(2:2:end) = qNoise0;
    
    r0 = real(antenna_0) .* iInterleaved0 + j * imag(antenna_0) .* qInterleaved0;
    
    
    iNoise1 = 1 * randn(numBits/2, 1);
    qNoise1 = 1 * randn(numBits/2, 1);
    
    iInterleaved1 = zeros(numBits, 1);
    iInterleaved1(1:2:end) = iNoise1;
    iInterleaved1(2:2:end) = iNoise1;
    
    qInterleaved1 = zeros(numBits, 1);
    qInterleaved1(1:2:end) = qNoise1;
    qInterleaved1(2:2:end) = qNoise1;
    
    r1 = real(antenna_1) .* iInterleaved1 + j * imag(antenna_1) .* qInterleaved1;
    
    r0_n = awgn(r0, SNR * 2);
    r1_n = awgn(r1, SNR * 2);
    r = r0_n + r1_n;
    
    r0_hat = r(1:2:end);
    r1_hat = r(2:2:end);
    
    s0_hat = complex(zeros(numBits, 1));
    
    ch0r0 =  iNoise0 .* real(r0_hat) - j * qNoise0 .* imag(r0_hat);
    
    
    h1cr1 = iNoise1 .* real(conj(r1_hat)) + j * qNoise1 .* ...
            imag(conj(r1_hat));
    

    ch1r0 = iNoise1 .* real(r0_hat) - j * qNoise1 .* imag(r0_hat);

    h0cr1 = iNoise0 .* real(conj(r1_hat)) + j * qNoise0 ...
            .* imag(conj(r1_hat));

    
    s0_hat(1:2:end) = ch0r0 +  h1cr1;
    
    s0_hat(2:2:end) = ch1r0 - h0cr1;
    
    
    
    demodulated = step(hDemod, s0_hat);
    errorStats = step(hErr, data, demodulated);
    berVec(i,:) = errorStats;
    %    fprintf('SNR: %f\n', snr(r0, n0);
    fprintf(['Bit error rate = %5.2e\nNumber of errors = %d\nTotal ' ...
             'bits = %d\n\n'], errorStats)

end
semilogy(SNRvec,berVec(:,1),'-gd', 'LineWidth',2)

SNRvec = 0:2.5:15;
berVec = zeros(length(SNRvec), 3);
fprintf('\n\n 2 Tx 2 Rx\n');
for i = 1:length(SNRvec)
    
    data = randi([0 1], numBits, 1);
    s0 = step(hMod, data);
    SNR = SNRvec(i);

    hErr = comm.ErrorRate;

    antenna_0 = s0;
    antenna_0(2:2:end) = -1 .* conj(s0(2:2:end));

    antenna_1 = complex(zeros(numBits, 1));
    antenna_1(1:2:end) = s0(2:2:end);
    antenna_1(2:2:end) = conj(s0(1:2:end));
        
    iNoise0 = 1 * randn(numBits/2, 1);
    qNoise0 = 1 * randn(numBits/2, 1);
    
    iInterleaved0 = zeros(numBits, 1);
    iInterleaved0(1:2:end) = iNoise0;
    iInterleaved0(2:2:end) = iNoise0;
    
    qInterleaved0 = zeros(numBits, 1);
    qInterleaved0(1:2:end) = qNoise0;
    qInterleaved0(2:2:end) = qNoise0;
    
    r0 = real(antenna_0) .* iInterleaved0 + j * imag(antenna_0) .* qInterleaved0;
    
    
    iNoise1 = 1 * randn(numBits/2, 1);
    qNoise1 = 1 * randn(numBits/2, 1);
    
    iInterleaved1 = zeros(numBits, 1);
    iInterleaved1(1:2:end) = iNoise1;
    iInterleaved1(2:2:end) = iNoise1;
    
    qInterleaved1 = zeros(numBits, 1);
    qInterleaved1(1:2:end) = qNoise1;
    qInterleaved1(2:2:end) = qNoise1;
    
    r1 = real(antenna_1) .* iInterleaved1 + j * imag(antenna_1) .* qInterleaved1;
    
    r0_n = awgn(r0, SNR * 2);
    r1_n = awgn(r1, SNR * 2);
    r = r0_n + r1_n;
    
    r0_hat = r(1:2:end);
    r1_hat = r(2:2:end);
    
    
    
    
    iNoise2 = 1 * randn(numBits/2, 1);
    qNoise2 = 1 * randn(numBits/2, 1);
    
    iInterleaved2 = zeros(numBits, 1);
    iInterleaved2(1:2:end) = iNoise2;
    iInterleaved2(2:2:end) = iNoise2;
    
    qInterleaved2 = zeros(numBits, 1);
    qInterleaved2(1:2:end) = qNoise2;
    qInterleaved2(2:2:end) = qNoise2;
    
    r2 = real(antenna_0) .* iInterleaved2 + j * imag(antenna_0) .* qInterleaved2;
    
    
    iNoise3 = 1 * randn(numBits/2, 1);
    qNoise3 = 1 * randn(numBits/2, 1);
    
    iInterleaved3 = zeros(numBits, 1);
    iInterleaved3(1:2:end) = iNoise3;
    iInterleaved3(2:2:end) = iNoise3;
    
    qInterleaved3 = zeros(numBits, 1);
    qInterleaved3(1:2:end) = qNoise3;
    qInterleaved3(2:2:end) = qNoise3;
    
    r3 = real(antenna_1) .* iInterleaved3 + j * imag(antenna_1) .* qInterleaved3;
    
    r2_n = awgn(r2, SNR * 2);
    r3_n = awgn(r3, SNR * 2);
    r_antenna_1 = r2_n + r3_n;
    
    r2_hat = r_antenna_1(1:2:end);
    r3_hat = r_antenna_1(2:2:end);
    
    
   
    
    
    s0_hat = complex(zeros(numBits, 1));
    
    ch0r0 =  iNoise0 .* real(r0_hat) - j * qNoise0 .* imag(r0_hat);
    
    
    h1cr1 = iNoise1 .* real(conj(r1_hat)) + j * qNoise1 .* ...
            imag(conj(r1_hat));
    

    ch1r0 = iNoise1 .* real(r0_hat) - j * qNoise1 .* imag(r0_hat);

    h0cr1 = iNoise0 .* real(conj(r1_hat)) + j * qNoise0 ...
            .* imag(conj(r1_hat));

    ch2r2 =  iNoise2 .* real(r2_hat) - j * qNoise2 .* imag(r2_hat);
    
    
    h3cr3 = iNoise3 .* real(conj(r3_hat)) + j * qNoise3 .* ...
            imag(conj(r3_hat));
    

    ch3r2 = iNoise3 .* real(r2_hat) - j * qNoise3 .* imag(r2_hat);

    h2cr3 = iNoise2 .* real(conj(r3_hat)) + j * qNoise2 ...
            .* imag(conj(r3_hat));

    
    s0_hat(1:2:end) = ch0r0 +  h1cr1 + ch2r2 + h3cr3;
    
    s0_hat(2:2:end) = ch1r0 - h0cr1 + ch3r2 - h2cr3;
    
    
    
    demodulated = step(hDemod, s0_hat);
    errorStats = step(hErr, data, demodulated);
    berVec(i,:) = errorStats;
    %    fprintf('SNR: %f\n', snr(r0, n0);
    fprintf(['Bit error rate = %5.2e\nNumber of errors = %d\nTotal ' ...
             'bits = %d\n\n'], errorStats)

end
semilogy(SNRvec,berVec(:,1),'-m^', 'LineWidth',2)

SNRvec = 0:2.5:12.5;
berVec = zeros(length(SNRvec), 3);


fprintf('\n\n 1 Tx 4 Rx\n');
for i = 1:length(SNRvec)
    
    data = randi([0 1], numBits, 1);
    s0 = step(hMod, data);
    SNR = SNRvec(i);

    hErr = comm.ErrorRate;

    iNoise0 = 1 * randn(numBits, 1);
    qNoise0 = 1 * randn(numBits, 1);
    r0 = real(s0) .* iNoise0 + j * imag(s0) .* qNoise0;
    
    iNoise1 = 1 * randn(numBits, 1);
    qNoise1 = 1 * randn(numBits, 1);

    r1 = real(s0) .* iNoise1 + j * imag(s0) .* qNoise1;
    
    
    iNoise2 = 1 * randn(numBits, 1);
    qNoise2 = 1 * randn(numBits, 1);

    r2 = real(s0) .* iNoise2 + j * imag(s0) .* qNoise2;
    
    iNoise3 = 1 * randn(numBits, 1);
    qNoise3 = 1 * randn(numBits, 1);

    r3 = real(s0) .* iNoise3 + j * imag(s0) .* qNoise3;
    
       
    r0_n = awgn(r0, SNR * 2);
    r1_n = awgn(r1, SNR * 2);
    r2_n = awgn(r2, SNR * 2);
    r3_n = awgn(r3, SNR * 2);
    
    r0_hat = real(r0_n).* iNoise0 - j * imag(r0_n) .* qNoise0;
    r1_hat = real(r1_n).* iNoise1  - j * imag(r1_n) .* qNoise1;    
    r2_hat = real(r2_n).* iNoise2  - j * imag(r2_n) .* qNoise2;
    r3_hat = real(r3_n).* iNoise3  - j * imag(r3_n) .* qNoise3;

    s0_hat = r0_hat + r1_hat + r2_hat + r3_hat;
    demodulated = step(hDemod, s0_hat);
    errorStats = step(hErr, data, demodulated);
    berVec(i,:) = errorStats;
    %    fprintf('SNR: %f\n', snr(r0, n0);
    fprintf(['Bit error rate = %5.2e\nNumber of errors = %d\nTotal ' ...
             'bits = %d\n\n'], errorStats)

end
semilogy(SNRvec,berVec(:,1),'-ks', 'LineWidth',2)
legend('no diversity (1 Tx, 1 Rx)','MMRC (1 Tx, 2 Rx)', ... 
       'new scheme (2 Tx 1 Rx)', 'new scheme (2 Tx 2 Rx)', ...
       'MMRC (1 Tx, 4 Rx');

