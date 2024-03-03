%==========================================================================
% Authors: Ali Waqar Azim, Ahmad Bazzi, Raed Shubair, and Marwa Chafii 
% Last Modified: March, 2022
% If you use this code or any (modified) part of it in any publication, please cite the paper: Azim, A. W., Bazzi, A., Shubair, R., & Chafii, M. (2022). Dual-Mode Chirp Spread Spectrum Modulation. IEEE Wireless Communications Letters.
% Contact person email: aliwaqarazim@gmail.com
%==========================================================================
clear all; close all; clc;

N = 64; % Number of available frequencies
EbN0_dB = 1:10; % SNR per bit on dB scale
EbN0_lin = 10.^(EbN0_dB/10); % SNR per bit on a linear scale
iter = 2e5; % Number of Monte Carlo iterations
nbits = 2*log2(N/2) + 2*log2(2) + log2(2); % Number of bits per symbol
% 2*log2(N/2): because of possible even and odd cyclic shifts
% log2(2): because of the use of either upchirp or down shirp
% 2*log2(2): because of bpsk alphabets
cu = exp(1i*(2*pi*(0:N-1).^2)./N); % Upchrip
cd = 1./cu; % Downchirp
se = nbits/N; % Spectral efficiency 
%==========================================================================
% Monte Carlo runs
%==========================================================================
for ii = 1:length(EbN0_lin) % Loop for SNR per bit
ii
EbN0 = EbN0_lin(ii);
    parfor kk = 1:iter % Loop for Monte Carlo runs
        % Defining variables to avoid warnings
        idx_rx_even_coh = [];
        idx_rx_odd_coh = [];
        idx_rx_even_nc = [];
        idx_rx_odd_nc = [];
        chirp_use_rx_coh = [];
        chirp_use_rx_nc = [];
        data_psk_even_rx_coh = [];
        data_psk_odd_rx_coh  = [];
        data_psk_even_rx_nc = [];
        data_psk_odd_rx_nc = [];
        RR_even = [];
        RR_odd = [];
        %==================================================================
        % TRANSMITTER
        %==================================================================
        chirp_use = randi([0 1],1,1); 
        % generate a random variable to determine weather upchirp or downchirp
        % will be used
        idx_tx_even = randi([1 N/2],1,1);
        % generate a random variable for even activated frequency shift
        idx_tx_odd = randi([1 N/2],1,1);
        % generate a random variable for even activated frequency shift
        data_psk_even = randi([0 1],1,1);
        % generate a random variable for the phase shift of the even activated frequency shift
        data_psk_odd = randi([0 1],1,1);
        % generate a random variable for the phase shift of the even activated frequency shift
        psk_data_even = real(pskmod(data_psk_even,2));
        % Phase shift alphabet for the even activated frequency shift
        psk_data_odd = real(pskmod(data_psk_odd,2));
        % Phase shift alphabet for the even activated frequency shift
        vec = zeros(1,N);
        vec(2*idx_tx_even) = psk_data_even;
        vec(2*idx_tx_odd-1) = psk_data_odd;
        % Generate a vector with even and odd frequency shift
        x = sqrt(N)*ifft(vec,N); % Unchiped Symbol
        if chirp_use == 0 % Chirped symbol
            s = cu.*x; 
        else
            s = cd.*x;
        end
        Es = sum(abs(s).^2); % Symbol energy
        Eb = Es/nbits; % Energy per bit
        n0 = Eb/EbN0; % Noise variance
        %==================================================================
        % CHANNEL
        %==================================================================
       
        n = sqrt(n0/2)*(randn(1,N) + 1i*randn(1,N)); % Noise vector
        r = s + n; % Received symbol
        %==================================================================
        % RECEIVER
        %==================================================================
        rtilde_1 = cd.*r; % Received waveform multiplied with downchirp
        rtilde_2 = cu.*r; % Received waveform multiplied with upchirp
        R1 = (1/sqrt(N))*fft(rtilde_1,N); % FFT of downchirped signal
        R2 = (1/sqrt(N))*fft(rtilde_2,N); % FFT of upchirped signal
        %=================================================================
        % Non-Coherent Detection
        %==================================================================
        [val1_nc,~] = max(abs(R1).^2);
        [val2_nc,~] = max(abs(R2).^2);
        if val1_nc > val2_nc % Determine chirp rate and the phase shift
            chirp_use_rx_nc = 0;
            RR = real(R1);
            RR_even = RR(2:2:end);
            RR_odd = RR(1:2:end);
        elseif val2_nc > val1_nc
            chirp_use_rx_nc = 1;
            RR = real(R2);
            RR_even = RR(2:2:end);
            RR_odd = RR(1:2:end);
        end
        [~,idx_rx_even_nc] = max(abs(RR_even));
        [~,idx_rx_odd_nc] = max(abs(RR_odd));
        data_psk_even_rx_nc = pskdemod(RR_even(idx_rx_even_nc),2);
        data_psk_odd_rx_nc = pskdemod(RR_odd(idx_rx_odd_nc),2);

        %==================================================================
        % BER Evaluation
        %==================================================================
       
        [bits_in_error_even_nc(kk),~]= biterr(idx_tx_even,idx_rx_even_nc);
        [bits_in_error_odd_nc(kk),~]= biterr(idx_tx_odd,idx_rx_odd_nc);
        [bits_in_error_data_even_nc(kk),~]= biterr(data_psk_even,data_psk_even_rx_nc);
        [bits_in_error_data_odd_nc(kk),~]= biterr(data_psk_odd,data_psk_odd_rx_nc);
        [bits_in_error_cu_nc(kk),~]= biterr(chirp_use,chirp_use_rx_nc);
    end
    
    ber_nc(ii) = (sum(bits_in_error_even_nc) + sum(bits_in_error_odd_nc) +...
        sum(bits_in_error_data_even_nc) + sum(bits_in_error_data_odd_nc) + ...
        sum(bits_in_error_cu_nc))/iter/nbits;
end
%


figure;
semilogy(EbN0_dB,ber_nc,'--o','LineWidth',2,'MarkerSize',8); grid on; hold on;
legend('Non-coherent','Interpreter','latex');
xlabel('$E_b/N_0$, dB','Interpreter','latex')
ylabel('BER','Interpreter','latex')
ylim([10^-5 1]); xlim([EbN0_dB(1) EbN0_dB(end)])

