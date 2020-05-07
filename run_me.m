% Space time trellis code (4 state)
% See Section 21.2.1 in the book "Quadrature Amplitude Modulation : From Basics 
% to Adaptive Trellis-Coded, Turbo-Equalised and Space-Time Coded OFDM, CDMA and 
% MC-CDMA Systems
% 2 X 2 MIMO system (QPSK modulation).

clear all
close all
clc
num_bit = 10^4; % number of data bits transmitted from each antenna
fade_var = 0.5; % fade variance of the channel per dimension
decoding_delay = 20; % decoding delay of the Viterbi algorithm

% SNR parameters
SNR_dB = 30; % SNR per bit (dB)
noise_var = 2*2*fade_var*2*2/(2*10^(0.1*SNR_dB)*1*2); % noise variance per dimension

% QPSK mapping
QPSK_SYM = [1+1i 1-1i -1+1i -1-1i];
QPSK_SYM1 = (1+1i)*ones(1,num_bit-1);
QPSK_SYM2 = (1-1i)*ones(1,num_bit-1);
QPSK_SYM3 = (-1+1i)*ones(1,num_bit-1);
QPSK_SYM4 = (-1-1i)*ones(1,num_bit-1);

tic()
% source
dk1 = randi([0 1],1,num_bit);
dk2 = randi([0 1],1,num_bit);

% Space time encoding
% output symbol from transmit antenna 1
tx_seq1_indices = mod(1*dk1(1:end-1)+2*dk2(1:end-1),4);
tx_seq1 = QPSK_SYM(tx_seq1_indices+1);

% output symbol from transmit antenna 2
tx_seq2_indices = mod(1*dk1(2:end)+2*dk2(2:end),4);
tx_seq2 = QPSK_SYM(tx_seq2_indices+1);

%------------ CHANNEL----------------------------------------------
% fade channel corresp. to transmit antenna 1 and receive antenna 1
h11 = normrnd(0,sqrt(fade_var),1,length(tx_seq1))+1i*normrnd(0,sqrt(fade_var),1,length(tx_seq1));

% fade channel corresp. to transmit antenna 2 and receive antenna 1
h21 = normrnd(0,sqrt(fade_var),1,length(tx_seq1))+1i*normrnd(0,sqrt(fade_var),1,length(tx_seq1));

% fade channel corresp. to transmit antenna 1 and receive antenna 2
h12 = normrnd(0,sqrt(fade_var),1,length(tx_seq1))+1i*normrnd(0,sqrt(fade_var),1,length(tx_seq1));

% fade channel corresp. to transmit antenna 2 and receive antenna 2
h22 = normrnd(0,sqrt(fade_var),1,length(tx_seq1))+1i*normrnd(0,sqrt(fade_var),1,length(tx_seq1));

% AWGN (note that lengths of tx_seq1 and tx_seq2 is same)
noise1 = normrnd(0,sqrt(noise_var),1,length(tx_seq1))+1i*normrnd(0,sqrt(noise_var),1,length(tx_seq1));
noise2 = normrnd(0,sqrt(noise_var),1,length(tx_seq1))+1i*normrnd(0,sqrt(noise_var),1,length(tx_seq1));

% channel output
chan_op = zeros(2,length(tx_seq1));
chan_op(1,:) = h11.*tx_seq1 + h21.*tx_seq2 + noise1;
chan_op(2,:) = h12.*tx_seq1 + h22.*tx_seq2 + noise2;

%--------------- RECEIVER  -----------------------------------
% Viterbi algorithm (VA)
% branch metrics for the VA
BM = zeros(16,num_bit-1);

BM(1,:)=abs(chan_op(1,:)-h11.*QPSK_SYM1-h21.*QPSK_SYM1+chan_op(2,:)-h12.*QPSK_SYM1-h22.*QPSK_SYM1).^2;
BM(2,:)=abs(chan_op(1,:)-h11.*QPSK_SYM1-h21.*QPSK_SYM2+chan_op(2,:)-h12.*QPSK_SYM1-h22.*QPSK_SYM2).^2;
BM(3,:)=abs(chan_op(1,:)-h11.*QPSK_SYM1-h21.*QPSK_SYM3+chan_op(2,:)-h12.*QPSK_SYM1-h22.*QPSK_SYM3).^2;
BM(4,:)=abs(chan_op(1,:)-h11.*QPSK_SYM1-h21.*QPSK_SYM4+chan_op(2,:)-h12.*QPSK_SYM1-h22.*QPSK_SYM4).^2;

BM(5,:)=abs(chan_op(1,:)-h11.*QPSK_SYM2-h21.*QPSK_SYM1+chan_op(2,:)-h12.*QPSK_SYM2-h22.*QPSK_SYM1).^2;
BM(6,:)=abs(chan_op(1,:)-h11.*QPSK_SYM2-h21.*QPSK_SYM2+chan_op(2,:)-h12.*QPSK_SYM2-h22.*QPSK_SYM2).^2;
BM(7,:)=abs(chan_op(1,:)-h11.*QPSK_SYM2-h21.*QPSK_SYM3+chan_op(2,:)-h12.*QPSK_SYM2-h22.*QPSK_SYM3).^2;
BM(8,:)=abs(chan_op(1,:)-h11.*QPSK_SYM2-h21.*QPSK_SYM4+chan_op(2,:)-h12.*QPSK_SYM2-h22.*QPSK_SYM4).^2;

BM(9,:)=abs(chan_op(1,:)-h11.*QPSK_SYM3-h21.*QPSK_SYM1+chan_op(2,:)-h12.*QPSK_SYM3-h22.*QPSK_SYM1).^2;
BM(10,:)=abs(chan_op(1,:)-h11.*QPSK_SYM3-h21.*QPSK_SYM2+chan_op(2,:)-h12.*QPSK_SYM3-h22.*QPSK_SYM2).^2;
BM(11,:)=abs(chan_op(1,:)-h11.*QPSK_SYM3-h21.*QPSK_SYM3+chan_op(2,:)-h12.*QPSK_SYM3-h22.*QPSK_SYM3).^2;
BM(12,:)=abs(chan_op(1,:)-h11.*QPSK_SYM3-h21.*QPSK_SYM4+chan_op(2,:)-h12.*QPSK_SYM3-h22.*QPSK_SYM4).^2;

BM(13,:)=abs(chan_op(1,:)-h11.*QPSK_SYM4-h21.*QPSK_SYM1+chan_op(2,:)-h12.*QPSK_SYM4-h22.*QPSK_SYM1).^2;
BM(14,:)=abs(chan_op(1,:)-h11.*QPSK_SYM4-h21.*QPSK_SYM2+chan_op(2,:)-h12.*QPSK_SYM4-h22.*QPSK_SYM2).^2;
BM(15,:)=abs(chan_op(1,:)-h11.*QPSK_SYM4-h21.*QPSK_SYM3+chan_op(2,:)-h12.*QPSK_SYM4-h22.*QPSK_SYM3).^2;
BM(16,:)=abs(chan_op(1,:)-h11.*QPSK_SYM4-h21.*QPSK_SYM4+chan_op(2,:)-h12.*QPSK_SYM4-h22.*QPSK_SYM4).^2;

% Viterbi algorithm
dec_ip=Viterbi_alg(BM,num_bit-1,decoding_delay);

% demapping indices to bits
dec_dk1 = zeros(1,num_bit-decoding_delay-1);
dec_dk2 = zeros(1,num_bit-decoding_delay-1);

for i1=1:num_bit-decoding_delay-1
temp = de2bi(dec_ip(i1),2);
dec_dk1(i1)=temp(2);
dec_dk2(i1)=temp(1);
end
toc()
% Bit error rate
BER = nnz(dk1(2:num_bit-decoding_delay)-dec_dk1+dk2(2:num_bit-decoding_delay)-dec_dk2)/...
    (2*(num_bit-decoding_delay-1))