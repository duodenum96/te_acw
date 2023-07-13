function cov = getcov(x, y, a, b, c, d)
%% Build covariance matrix with estimated coefficients. 
% Everything are column matrices
% Input: 
%   x and y: Time series
%   a (y to y), b (x to y), c (x to x) and d (y to x): AR coefficients
% Output: 
%   cov: Covariance matrix (see documentation 'AR OLS solution.docx')
% Authored by Yasir Çatal a.k.a. Duodenum

x = x(:);
y = y(:);

q_xx = length(c);
q_xy = length(b);
q_yx = length(d);
q_yy = length(a);
N = length(x);

% Estimate yhat and xhat

q1 = max([q_yy, q_xy]);
y_hat = y(1:q1);

for i = (q1+1):N
    y_hat(i) = a' * flip(y_hat( (i-q_yy):(i-1) ) ) + b' * flip(x( (i-q_xy+1):(i) ) );
end

q2 = max([q_yx, q_xx]);
x_hat = x(1:q2);

for i = (q2+1):N
    x_hat(i) = c' * flip(x_hat( (i-q_xx):(i-1) ) ) + d' * flip(y( (i-q_yx+1):(i) ) );
end

% Estimate covariances
q = max([q1, q2]);

x_hat = x_hat(:);
y_hat = y_hat(:);

cov_yy = (1 / (N-q)) * sum((y((q+1):N) - y_hat((q+1):N)).^2);
cov_xx = (1 / (N-q)) * sum((x((q+1):N) - x_hat((q+1):N)).^2);
cov_xy = (1 / (N-q)) * sum( ( x((q+1):N) - x_hat((q+1):N) ) .* ( y((q+1):N) - y_hat((q+1):N) ));

cov = [cov_xx, cov_xy; cov_xy, cov_yy];
end