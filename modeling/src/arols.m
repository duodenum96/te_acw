function [a, b] = arols(x, y, q_x, q_y)
    %% OLS solution to ar equation with exogenous input
    % Allows flexibility of selecting different orders for variables x and y
    % Input:
    %    x: Exogenous input
    %    y: Output
    %    q_x: Order of x (i-0:i-q_x) (q_xy in bivariate version)
    %    q_y: Order of y (i-1:i-q_y) (q_yy in bivariate version)
    % Output: 
    %    a: AR coefficients of y
    %    b: AR coefficients of x
    %    cov: Covariance matrix of the solution
    % Authored by Yasir Çatal a.k.a. Duodenum
    
    assert(length(x) == length(y), ...
        'Inputs x and y must be of same length');
    
    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    
    t = max([q_x, q_y]) + 1;
    N = length(y);
    Y = y(t:N);
    Z = [];
    
    for i = 1:q_y
        Z = [Z, y( (t-i):(N-i) )];
    end
    
    for i = 1:q_x
        Z = [Z, x( (t-i+1):(N-i+1) )];
    end
    
    % phi = inv(Z' * Z) * Z' * Y; % Usual way
    phi = Z \ Y; % With backslash operator
    
    a = phi(1:q_y);
    b = phi(q_y+1:end);
    
    end
    