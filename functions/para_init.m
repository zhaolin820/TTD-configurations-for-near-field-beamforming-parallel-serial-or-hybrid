function [para] = para_init()
%Construct a struct of the initial values for all the parameters 
%  [values] = para_init()
%Inputs:
%   None
%Outputs:
%   values: a struct
%Date: 04/04/2024
%Author: Zhaolin Wang

para.fc = 1e11; % carrier frequency 100 GHz
B = 1e10; % 10 GHz system bandwidth 
c = 3e8;
para.lambda = c/para.fc;
para.d = para.lambda/2;

para.Pt = 10^(20/10); % dBm overall transmit power
para.K = 4; % user number

para.N = 512; % number of antennas
para.N_RF = para.K; % number of RF chains
para.N_T = 32; % number of TTDs connected to each RF chain
para.t_max = para.N/para.N_T * para.d/c;


para.D = (para.N-1)*para.d;
para.Rayl = 2*para.D^2/para.lambda;

para.M = 10; % number of subcarriers
para.Lcp = 4;

% subcarrier frequency
m = 1:para.M;
para.fm_all =  para.fc + B*(2*m-1-para.M) / (2*para.M);

para.noise_dB = -174; % noise power in dBm/Hz
para.noise_dB = para.noise_dB + 10*log10(B/10); % noise power in dBm

para.Gt = 10; % dBi, transmit antenna gain
para.Gr = 5; % dBi, receive antenna gain


end

