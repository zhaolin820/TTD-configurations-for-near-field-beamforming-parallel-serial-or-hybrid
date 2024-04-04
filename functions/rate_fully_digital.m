function [R_sum, R] = rate_fully_digital(para, P, h)
%Calculate spectral efficiency
%  [R_sum, R] = rate_fully_digital(para, P, h)
%Inputs:
%   para: structure of the initial parameters
%   P: beamformer
%   h: channel for all users
%Outputs:
%   R_sum: spectral efficiency
%   R: rate of each user at each subcarrier
%Date: 04/04/2024
%Author: Zhaolin Wang


R = zeros(para.K, para.M);
for m = 1:para.M
    Pm = P(:,:,m);
    for k = 1:para.K
        hmk = h(:,k,m);
        pmk = Pm(:,k); 
        P_I = Pm; P_I(:,k) = [];
        Imk = norm(hmk'*P_I)^2 + norm(Pm, 'fro')^2/para.Pt; 
        R(k,m) = log2( 1 + abs(hmk'*pmk)^2/Imk );
    end
end
R_sum = sum(sum(R));


end

