function [a, r_ave] = beta_func(r, theta, N, d, f)
%The vector used for initialize the analog beamforming
%  [a, r_ave] = beamfocusing(r, theta, N, d, f)
%Inputs:
%   r: distance
%   theta: angle of direction
%   N: number of antennas
%   d: antenna spacing
%   f: carrier frequency
%Outputs:
%   a: modified array response vector
%   r_ave: average distance for all antenna arrays
%Date: 04/04/2024
%Author: Zhaolin Wang

c = 3e8;

n = 0:(N-1);
delta_n = (n-(N-1)/2) * d;
delta_n = delta_n';

% distance
r_n = sqrt(r^2 + delta_n.^2 - 2*r*delta_n*cos(theta)) - r;
r_ave = mean(r_n);
r_n = r_n - r_ave;
% beamfocusing vector
a = exp( -1i * 2 * pi * f/c * r_n );

end

