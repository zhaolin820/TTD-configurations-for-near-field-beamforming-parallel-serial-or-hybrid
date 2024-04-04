function [h] = generate_channel(para, user_r, user_theta)
%Generate channels
%  [h] = generate_channel(para, user_r, user_theta)
%Inputs:
%   para: structure of the initial parameters
%   user_r: distance for all users
%   user_theta: angle for all users
%Outputs:
%   h: channel for all users
%Date: 04/04/2024
%Author: Zhaolin Wang

HITRANparams = importdata('data_freq_abscoe.txt');

L = 3;
r_NLoS = rand(para.K, L) * 10 + 5; % 5 ~ 15 m
theta_NLoS = rand(para.K, L) * pi; % 0 ~ 180 degree


h = zeros(para.N, para.K, para.M);
for k = 1:para.K
    for m = 1:para.M
        fm = para.fm_all(m);
        
        % LoS channel
        path_loss = getSpreadLoss(fm, user_r(k)) + getAbsLoss(fm, user_r(k), HITRANparams );
        path_loss = 10.^((-path_loss - para.noise_dB + para.Gt + para.Gr)/10);
        beta_k =  path_loss;
        h(:,k,m) = beta_k * beamfocusing(user_r(k), user_theta(k), para.N, para.d, fm);
        
        % NLoS channel
        path_loss = getSpreadLoss(fm, user_r(k)) + getAbsLoss(fm, user_r(k), HITRANparams );
        path_loss = 10.^((-path_loss - para.noise_dB + para.Gt + para.Gr - 15)/10);
        for l = 1:L
            beta_k_l = path_loss*exp(1i*rand(1) * pi);
            h(:,k,m) = h(:,k,m) + beta_k_l * beamfocusing(r_NLoS(k,l), theta_NLoS(k,l), para.N, para.d, fm);
        end
    end
end
h = conj(h);

end