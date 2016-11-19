% function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)
% This function should perform the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           3X1
%           sigma_bar(t)        3X3
%           Q                   2X2
%           z_i(t)              2X1
%           M                   2XN
%           Lambda_m            1X1
% Outputs: 
%           c(t)                1X1
%           outlier             1X1
%           nu^i(t)             2XN
%           S^i(t)              2X2XN
%           H^i(t)              2X3XN
function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)
% FILL IN HERE
N = size(M, 2);
nu = zeros(2, N);
H = zeros(2, 3, N);
S = zeros(2, 2, N);
phi = zeros(1, N);
D = zeros(1, N)';

for j = 1 : N
    z_ihat = observation_model(mu_bar, M, j);
    H(:, :, j) = jacobian_observation_model(mu_bar, M, j, z_ihat, 1);
    S(:, :, j) = H(:, :, j) * sigma_bar * H(:, :, j)' + Q;
    nu(:, j) = z_i - z_ihat;
    nu(2, j) = mod(nu(2, j) + pi, 2 * pi) - pi; %Attention!
    D(j) = nu(:, j)' * inv(S(:, :, j)) * nu(:, j);
    phi(j) = det(2 * pi * S(:, :, j))^(-1/2) * exp(-(1/2) * D(j));
end
c = find(phi==max(phi));
outlier = logical(D(c) >= Lambda_m);
end