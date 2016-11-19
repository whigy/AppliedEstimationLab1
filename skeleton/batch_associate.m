% function [c,outlier, nu_bar, H_bar] = batch_associate(mu_bar,sigma_bar,z,M,Lambda_m,Q)
% This function should perform the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           3X1
%           sigma_bar(t)        3X3
%           Q                   2X2
%           z(t)                2Xn
%           M                   2XN
%           Lambda_m            1X1
% Outputs: 
%           c(t)                1Xn
%           outlier             1Xn
%           nu_bar(t)           2nX1
%           H_bar(t)            2nX3
function [c,outlier, nu_bar, H_bar] = batch_associate(mu_bar,sigma_bar,z,M,Lambda_m,Q)
% FILL IN HERE
n = size(z, 2);
c = zeros(1, n);
outlier = zeros(1, n);
nu_bar = zeros(2 * n, 1);
H_bar = zeros(2 * n, 3);

for i = 1 : n
    [cT,outlierT, nu_barT, S, H_barT] = associate(mu_bar, sigma_bar, z(:, i) , M, Lambda_m, Q);
    c(i) = cT;
    outlier(i) = outlierT;
    nu_barT2 = nu_barT(:, cT);
    nu_barT2(2, :) = mod(nu_barT2(2, :) + pi, 2 * pi) - pi; %Attention!
    ind = (i - 1) * 2 + 1;
    nu_bar(ind:ind+1, :) = nu_barT2;
    H_bar(ind:ind+1, :) = H_barT(:, :, cT);
end

end