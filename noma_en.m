clear all
%%%%%%%%%%%%% Introduction %%%%%%%%%%%%%%%%%
% This code was modified by Mr. Wang from Anhui University of Technology based on the source code of "Introduction to MATLAB Communication System Modeling and Simulation" (Electronic Industry Press, Liu Xueyong), in order to realize the design of NOMA communication system based on convolutional codes, QPSK, OFDM.

%%%%%%%%%%%%% Parameter Setting Section %%%%%%%%%%%%%%%%%

Nsp=52; %Number of subcarriers of the system (excluding DC subcarrier)
Nfft=64; % FFT length
Ncp=16; % Cyclic prefix length
Ns=Nfft+Ncp; % Length of one complete OFDM symbol
noc=53; % Total number of subcarriers including DC subcarrier
Nd=6; % Number of OFDM symbols per frame (excluding training symbols)
M1=4; % QPSK modulation
sr=250000; % OFDM symbol rate
EbNo=0:2:30; % Normalized SNR
Nfrm=10000; % Number of simulation frames for each SNR
SF = 7 ;
BW = 5e3 ;
fc = 50e3 ;
Pt = 15 ;
Fs = 150e3 ;
Fc = 55e3 ;
% Doppler effect is not considered for now
ts=1/sr/Ns; % OFDM symbol sampling time interval
t=0:ts:(Ns*(Nd+1)*Nfrm-1)*ts; % Sampling time
fd=10; % Maximum Doppler frequency shift
% Convolutional code parameters setting
L=7; % Constraint length of convolutional code
tblen=6*L; % Traceback depth of Viterbi decoder
stage = 3; % Order of m sequence
ptap1 = [1 3]; % Connections of m sequence register
regi1 = [1 1 1]; % Initial value of m sequence register
% Power parameters for two users
Rp_db=0; % Power allocation ratio in dB between user 1 and 2
Rp=10^(Rp_db/10); % Power allocation ratio between user 1 and 2, Rp=p_1/p_2, assuming user 1 gets more power
p_u1=Rp/(1+Rp); % Power allocation coefficients p_1 and p_2 for user 1 and 2 respectively, assuming p_1 + p_2 = 1
p_u2=1/(1+Rp);
%% Rayleigh channel generation
h=rayleigh(fd,t); % Generate single-path Rayleigh fading channel
h1=sqrt(2/3)*h;
h2=sqrt(1/3)*rayleigh(fd,t);
h2=[zeros(1,4) h2(1:end-4)];

%%
% Frequency domain data of training symbols, using long training symbols data of 802.11a
Preamble=[1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 ...
1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1]; %52-bit data
Preamble1=zeros(1,Nfft);
Preamble1(2:27)=Preamble(27:end); % Rearrange the training symbol data
Preamble1(39:end)=Preamble(1:26);

preamble1=ifft(Preamble1); % Time domain data of training symbol
preamble1=[preamble1(Nfft-Ncp+1:end) preamble1]; % Add cyclic prefix to training symbol

%% %%%%%%%%%%%%%%%%%%% Simulation Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:length(EbNo)
%% Generate random sequences for users
data_tx1=randsrc(1,Nsp*Nd*Nfrm,[0:1]); %User 1
data_tx2=randsrc(1,Nsp*Nd*Nfrm,[0:1]); %User 2
data_tx1yanzhen=reshape(data_tx1,52,[]); %Variable for BER calculation
data_tx2yanzhen=reshape(data_tx2,52,[]);
%% Channel encoding (convolutional code or interleaver)
%Convolutional code: Non-linear forward error correction code
%Interleaver: Maximally disperses burst errors

trellis = poly2trellis(7,[133 171]);             % (2,1,7) convolutional code with rate 1/2  
code_data=convenc(data_tx1,trellis);             % Input (1*3120000), output (1*6240000)   
trellis = poly2trellis(7,[133 171]);    
code_data2=convenc(data_tx2,trellis);    
%% Verify convolutional code? Remove if correct
trellis = poly2trellis(7,[133 171]);  
yanzheng_data = vitdec(code_data,trellis,tblen,'trunc','hard');                 
yanzheng_shuchu=reshape(yanzheng_data,52,[]);
%Verification correct √
%% QPSK modulation
% data1=pskmod(msg1,M1,pi/4); %User 1 QPSK modulation
%data2=pskmod(msg2,M1,pi/4); %User 2 QPSK modulation
%User 3 data conversion and QPSK modulation
data_temp1=reshape(code_data,2,[])'; %
data_temp2=bi2de(data_temp1); %Binary to decimal for QPSK modulation
data_temp3=pskmod(data_temp2,M1,pi/M1); %4PSK modulation
data_temp4=reshape(data_temp3,52,[]) ; %User 3 (convolutional coded) QPSK modulated data
%User 4 data conversion and QPSK modulation
data_temp11=reshape(code_data2,2,[])';            %
    data_temp22=bi2de(data_temp11);                   % 
    data_temp33=pskmod(data_temp22,M1,pi/M1);         % 
    data_temp44=reshape(data_temp33,52,[]) ; 
%Constellation diagram verification 
% data_xingzuotu1=reshape(data1,1,[]);
% data_xingzuotu2=reshape(data2,1,[]);
% data_xingzuotu3=reshape(data_temp4,1,[]);
%% User superimposition
data_diejia=sqrt(p_u1)*data_temp4+sqrt(p_u2)*data_temp44;
%% ifft data rearrangement
%data3=zeros(Nfft,Nd*Nfrm);                   % According to FFT requirements, rearrange data 
%data4=zeros(Nfft,Nd*Nfrm);                   
data5=zeros(Nfft,Nd*Nfrm);
data6=zeros(Nfft,Nd*Nfrm); 
data7=zeros(Nfft,Nd*Nfrm);
%User 1 rearrange data
% data3(2:27,:)=data1(27:end,:); % Rearrange user 1 QPSK data
% data3(39:end,:)=data1(1:26,:);
%User 2 rearrange data
% data4(2:27,:)=data2(27:end,:); % Rearrange user 2 QPSK data
% data4(39:end,:)=data2(1:26,:);

%User 3 (convolutional coded) rearrange data
data5(2:27,:)=data_temp4(27:end,:);          % Rearrange user 3 QPSK data
data5(39:end,:)=data_temp4(1:26,:);          
%User 4 (convolutional coded) rearrange data 
data6(2:27,:)=data_temp44(27:end,:);          % Rearrange user 4 QPSK data
data6(39:end,:)=data_temp44(1:26,:);         
 %Superimposed signal (convolutional coded) rearrange data

 data7(2:27,:)=data_diejia(27:end,:);          % Rearrange superimposed QPSK data
 data7(39:end,:)=data_diejia(1:26,:);          
 
%clear data1 data2 data11;                   % Clear unnecessary temporary variables

%data3=ifft(data3);                           % User 1 IFFT   
%data4=ifft(data4);                           % User 2 IFFT  
data5=ifft(data5);                           % User 3 (convolutional coded) ifft
data6=ifft(data6);                           % User 4 (convolutional coded) ifft
data7=ifft(data7);
%% Add cyclic prefix after ifft
% data3=[data3(Nfft-Ncp+1:end,:);data3]; % Add cyclic prefix, length 16
% data4=[data4(Nfft-Ncp+1:end,:);data4];
data5=[data5(Nfft-Ncp+1:end,:);data5];
data6=[data6(Nfft-Ncp+1:end,:);data6];
data7=[data7(Nfft-Ncp+1:end,:);data7];

%spow1=norm(data3,'fro').^2/(Nsp*Nd*Nfrm);    % Calculate symbols power of user 1
%spow2=norm(data4,'fro').^2/(Nsp*Nd*Nfrm);    % Calculate symbols power of user 2
spow3=norm(data5,'fro').^2/(Nsp*Nd*Nfrm);    % Calculate symbols power of user 3 
spow4=norm(data6,'fro').^2/(Nsp*Nd*Nfrm);    % Calculate symbols power of user 4
spow5=norm(data7,'fro').^2/(Nsp*Nd*Nfrm);    %norm calculates the norm of the QPSK data matrix, calculating the symbol energy

%data10=zeros(Ns,(Nd+1)*Nfrm);
%data11=zeros(Ns,(Nd+1)*Nfrm);
data12=zeros(Ns,(Nd+1)*Nfrm);
data13=zeros(Ns,(Nd+1)*Nfrm);
data14=zeros(Ns,(Nd+1)*Nfrm);
%% Insert training symbol every 7 frames
for indx=1:Nfrm


%data10(:,(indx-1)*(Nd+1)+1)=preamble1.';   %Add training sequence for user 1
%data10(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data3(:,(indx-1)*Nd+1:indx*Nd); %Add symbol info for user 1

%data11(:,(indx-1)*(Nd+1)+1)=preamble1.';    %Add training sequence for user 2
%data11(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data4(:,(indx-1)*Nd+1:indx*Nd); %

data12(:,(indx-1)*(Nd+1)+1)=preamble1.';    %Add training sequence for user 3 (convolutional coded)    
data12(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data5(:,(indx-1)*Nd+1:indx*Nd); %

data13(:,(indx-1)*(Nd+1)+1)=preamble1.';    %Add training sequence for user 4 (convolutional coded)    
data13(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data6(:,(indx-1)*Nd+1:indx*Nd); %

data14(:,(indx-1)*(Nd+1)+1)=preamble1.';    %Add training sequence for superimposed signal (convolutional coded)    
data14(:,(indx-1)*(Nd+1)+2:indx*(Nd+1))=data7(:,(indx-1)*Nd+1:indx*Nd);
end
%% Serial to parallel conversion

% data10=reshape(data10,1,Ns*(Nd+1)Nfrm); % User 1
%data11=reshape(data11,1,Ns(Nd+1)Nfrm); % User 2
data12=reshape(data12,1,Ns*(Nd+1)*Nfrm); % User 3 (convolutional coded)
data13=reshape(data13,1,Ns*(Nd+1)*Nfrm); % User 4 (convolutional coded)
data14=reshape(data14,1,Ns*(Nd+1)*Nfrm); % Superimposed users

%data51=zeros(1,length(data10)); % User 1
%data511=zeros(1,length(data11)); % User 2
data71=zeros(1,length(data12)); % User 3 (convolutional coded)
data81=zeros(1,length(data13)); % User 4 (convolutional coded)
data91=zeros(1,length(data14)); % Superimposed users

%data51(5:end)=data10(1:end-4); % Second path is delayed 4 samples from first path
%data511(5:end)=data11(1:end-4);
data71(5:end)=data12(1:end-4);
data81(5:end)=data13(1:end-4);
data91(5:end)=data14(1:end-4);

%% Calculate noise standard deviation based on Eb/No

%sigma1=sqrt(1/2*spow1/log2(M1)10.^(-EbNo(ii)/10));
%sigma2=sqrt(1/2*spow2/log2(M1)10.^(-EbNo(ii)/10));
sigma3=sqrt(1/2*spow3/log2(M1)*10.^(-EbNo(ii)/10));
sigma4=sqrt(1/2*spow4/log2(M1)*10.^(-EbNo(ii)/10));
sigma5=sqrt(1/2*spow5/log2(M1)*10.^(-EbNo(ii)/10));

%% Channel
for indx=1:Nfrm

%dd_user1=data10((indx-1)*Ns(Nd+1)+1:indx*Ns*(Nd+1));
%dd_user2=data11((indx-1)*Ns(Nd+1)+1:indx*Ns*(Nd+1));
dd_user3=data12((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
dd_user4=data13((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
dd_user5=data14((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
   %dd_user11=data51((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1)); 
   %dd_user22=data511((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
   dd_user33=data71((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
   dd_user44=data81((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
   dd_user55=data91((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
   
   % Current frame's single path channel parameters  
   hh=h((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));      % Not used currently, this system generates dual path channel
   % Current frame's dual path channel parameters
   hh1=h1((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));    
   hh2=h2((indx-1)*Ns*(Nd+1)+1:indx*Ns*(Nd+1));
%  % The signal passes through the dual-path Rayleigh fading channel and adds white Gaussian noise
    %r11=hh1.*dd_user1+hh2.*dd_user11+sigma1*(randn(1,length(dd_user1))+1i*randn(1,length(dd_user1)));   
    %r31=hh1.*dd_user2+hh2.*dd_user22+sigma2*(randn(1,length(dd_user2))+1i*randn(1,length(dd_user2)));
    r41=hh1.*dd_user3+hh2.*dd_user33+sigma3*(randn(1,length(dd_user3))+1i*randn(1,length(dd_user3)));
    r51=hh1.*dd_user4+hh2.*dd_user44+sigma4*(randn(1,length(dd_user4))+1i*randn(1,length(dd_user4)));
    r61=hh1.*dd_user5+hh2.*dd_user55+sigma5*(randn(1,length(dd_user5))+1i*randn(1,length(dd_user5)));
    %% Receiver serial to parallel conversion
    % Dual path propagation 
    %r11=reshape(r11,Ns,Nd+1);
    %r31=reshape(r31,Ns,Nd+1); 
    r41=reshape(r41,Ns,Nd+1);
    r51=reshape(r51,Ns,Nd+1);
    r61=reshape(r61,Ns,Nd+1);

     % Remove cyclic prefix from dual path    
     % r11=r11(Ncp+1:end,:);
     %r31=r31(Ncp+1:end,:);
     r41=r41(Ncp+1:end,:);  
     r51=r51(Ncp+1:end,:);
     r61=r61(Ncp+1:end,:);
     %% Receiver FFT

      % Dual path FFT
      %R11=fft(r11);
      %R31=fft(r31);
      R41=fft(r41);
      R51=fft(r51);
      R61=fft(r61);

      % Dual path data rearrangement  
      %R11=[R11(39:end,:);R11(2:27,:)]; 
      %R31=[R31(39:end,:);R31(2:27,:)];
      R41=[R41(39:end,:);R41(2:27,:)];
      R51=[R51(39:end,:);R51(2:27,:)]; 
      R61=[R61(39:end,:);R61(2:27,:)];
     %% Receiver channel estimation

      % Channel estimation
      % HH11=(Preamble.')./R11(:,1);            
      % HH31=(Preamble.')./R31(:,1);
      HH41=(Preamble.')./R41(:,1);
      HH51=(Preamble.')./R51(:,1);
      HH61=(Preamble.')./R61(:,1);

       %HH11=HH11*ones(1,Nd);                       %User 1 channel estimation
       %HH31=HH31*ones(1,Nd);                       %User 2 channel estimation 
       HH41=HH41*ones(1,Nd);                       %User 3 (convolutional coded)
       HH51=HH51*ones(1,Nd);  
       HH61=HH61*ones(1,Nd);
      %% Channel compensation

       %x1=R11(:,2:end).*HH11;                      %Dual path compensation for user 1       
       %x2=R31(:,2:end).*HH31;                      %Dual path compensation for user 2
       x3=R41(:,2:end).*HH41;                      
       x4=R51(:,2:end).*HH51;
       x5=R61(:,2:end).*HH61;
       %% OMA data demodulation

    %x1=pskdemod(x1,M1,pi/4);  
    %x2=pskdemod(x2,M1,pi/4);  
    x3=pskdemod(x3,M1,pi/4);
    x4=pskdemod(x4,M1,pi/4);
    x5_3=pskdemod(x5,M1,pi/4);

   %% NOMA demodulation for User 3

    % QPSK demodulation for user 3 is complete
  	x5_3= reshape(x5_3,[],1);
    x5_3= de2bi(x5_3);
    x5_3= reshape(x5_3.',1,[]);                                        % Superimposed signal for user 3 complete

    x5_decode = vitdec(x5_3,trellis,tblen,'trunc','hard');             % Bitstream information for user 3

    x5_3_bit=reshape(x5_decode.',Nsp,[]);

     %% NOMA demodulation for User 4

    x5_4=convenc(x5_decode ,trellis);
    x5_4=reshape(x5_4,2,[])';                                            % Reconstruct the encrypted data into a []*2 matrix, m=4
	
	x5_4=bi2de(x5_4);                                                    % Binary to decimal for QPSK modulation

    x5_4=reshape(x5_4,1,[]);
    x5_4_temp=pskmod(x5_4,M1,pi/M1);                                  

   x5_4_temp=reshape(x5_4_temp.',Nsp,[]);

    % Subtract superimposed signal  
    x5_4_temp2=x5-x5_4_temp;
    x5_4_temp2=reshape(x5_4_temp2,1,[]);

    x5_4_temp3=pskdemod(x5_4_temp2,M1,pi/M1);

    x5_4_temp4= reshape(x5_4_temp3,[],1);

    x5_4out= de2bi(x5_4_temp4);

    x5_4out2= reshape(x5_4out.',1,[]);% Superimposed signal for user 4 complete

    x5_4decode = vitdec(x5_4out2,trellis,tblen,'trunc','hard');    

    x5_4_bit=reshape(x5_4decode.',Nsp,[]);

    
    %% OMA demodulation for User 3

    x3 = reshape(x3,[],1);

    x3 = de2bi(x3); 

    x3 = reshape(x3',1,[]);

    % Viterbi decoding

    trellis = poly2trellis(7,[133 171]);

    x3 = vitdec(x3,trellis,tblen,'trunc','hard');                  

    x3=reshape(x3,52,[]);
%% OMA demodulation for User 4
   x4= reshape(x4,[],1);
   x4 = de2bi(x4);
   x4= reshape(x4',1,[]);

   % Viterbi decoding
   trellis = poly2trellis(7,[133 171]);  
   x4 = vitdec(x4,trellis,tblen,'trunc','hard');                  
   x4=reshape(x4,52,[]);

   %% Count the number of error bits in a frame

    % Theoretical case
    [neb1(indx),temp]=biterr(x3, data_tx1yanzhen(:,(indx-1)*Nd+1:indx*Nd),log2(M1));            % Bit errors for user 1 on dual channel

    [neb2(indx),temp]=biterr(x4,data_tx2yanzhen(:,(indx-1)*Nd+1:indx*Nd),log2(M1));            % Bit errors for user 2 on dual channel

    [neb3(indx),temp]=biterr(x5_3_bit,data_tx1yanzhen(:,(indx-1)*Nd+1:indx*Nd),log2(M1));     % Bit errors for user 3 on dual channel 

    [neb4(indx),temp]=biterr(x5_4_bit,data_tx2yanzhen(:,(indx-1)*Nd+1:indx*Nd),log2(M1));           % Bit errors for user 2 on single channel
end

%% BER calculation

ber1(ii)=sum(neb1)/(Nsp*log2(M1)*Nd*Nfrm);                              % BER for user 1
ber2(ii)=sum(neb2)/(Nsp*log2(M1)*Nd*Nfrm);                              % BER for user 2                                 
ber3(ii)=sum(neb3)/(Nsp*log2(M1)*Nd*Nfrm);                              % BER for user 3     
ber4(ii)=sum(neb4)/(Nsp*log2(M1)*Nd*Nfrm);                               % BER for user 3
end

%% Constellation diagram display

%figure(1);stem(msg);
%msg1 = pskmod(msg,4,pi/4); % 4psk modulation with initial phase pi/4
%scatterplot(msg1); axis([-1.2,1.2,-1.2,1.2]);% Draw constellation diagram
%hold on;

%% Plotting

figure(2)
semilogy(EbNo,ber1,'-bo',EbNo,ber2,'-rv',EbNo,ber3,'-y*',EbNo,ber4,'-gd')
grid on
title('BER comparison of NOMA and OFDM systems')
legend('OFDM User 3','OFDM User 4','NOMA User 3','NOMA User 4')
xlabel('Eb/No[dB]')
ylabel('BER')