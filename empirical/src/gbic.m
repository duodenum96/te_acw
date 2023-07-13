function GBIC = gbic(N, cov, q)
%% Calculate generalized Bayesian Information Criterion
% Input: 
    % N: Number of time points
    % cov: Covariance Matrix
    % q: Order parameters of different variables
% Output: 
    % GBIC: Generalized Bayesian Information Criterion
% Authored by Yasir Çatal a.k.a. Duodenum

GBIC = N * log(det(cov)) + sum(q) * log(N);
end