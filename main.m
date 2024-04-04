clc
clear all
close all

addpath("functions\");

% parameters
para = para_init();
user_r = rand(para.K, 1) * 10 + 5; % 5 ~ 15 m
user_theta = sort(rand(para.K, 1) * pi); % 0 ~ 180 degree

% generate channel
[h] = generate_channel(para, user_r, user_theta);

% initialize fully digital beamformers
P_initial = randn(para.N, para.K) + 1i * randn(para.N, para.K);
P_initial = P_initial / norm(P_initial, 'fro') * sqrt(para.Pt);


% algorithms
[R, P] = algorithm_fully_digital(para, h, P_initial); 
[R_s] = algorithm_TTD_serial(para, h, P, user_r, user_theta);
[R_h] = algorithm_TTD_hybrid(para, h, P, user_r, user_theta); 
[R_p] = algorithm_TTD_parallel(para, h, P, user_r, user_theta); 

% convergence behavior
figure; hold on; box on;
plot(R_s, '-r', 'LineWidth', 1.5);
plot(R_h, '-.b', 'LineWidth', 1.5);
plot(R_p, ':k', 'LineWidth', 1.5);
legend('Serial', 'Hybrid', 'Parallel');
xlabel('Number of cumulative BCD iterations');
ylabel('Spectral efficiency (bit/s/Hz)');