function [q_xx, q_xy, q_yx, q_yy] = greedy(x, y, qmax)
    %% Use the an algorithm described in Yang et al. 2013 (IEEE Transactions on Biomedical Engineering) to 
    % find the best model order parameters (minimizing gBIC)
    % Input: 
    %     x: Time series (exogenous input)
    %     y: Time series (output)
    %     qmax: Maximum model order (defaults to 20)
    % Output: 
        % q_ij: Optimal model order from i to j
    % Authored by Yasir Çatal a.k.a. Duodenum
    
    assert(length(x) == length(y), ...
        'Inputs x and y must be of same length');
    
    N = length(x);
    
    % Step 0: Search through qmax to find the best qmax (assume both variates have same q, 
    % update qmax to the min gBIC value)
    
    for i = 1:qmax
        [a, b] = arols(x, y, i, i);
        [c, d] = arols(y, x, i, i);
        cov = getcov(x, y, a, b, c, d);
        bic(i) = gbic(N, cov, [i, i, i, i]);
    end
    
    [~, indx] = min(bic);
    qmax = indx;
    
    % Step 1: Set parameters of eq. 2 (q_xx, q_yx) to qmax and find the best parameters for eq. 1 
    % (q_yy, q_xy)
    q_xx = qmax;
    q_yx = qmax;
    bic = [];
    
    for q_yy = 1:qmax
        for q_xy = 1:qmax
            [a, b] = arols(x, y, q_xy, q_yy);
            [c, d] = arols(y, x, q_yx, q_xx);
            cov = getcov(x, y, a, b, c, d);
            bic(q_yy, q_xy) = gbic(N, cov, [q_xx, q_xy, q_yx, q_yy]);
        end
    end
    
    % Step 2: Set parameters of eq. 1 (q_yy, q_xy) to optimal ones and do the grid search for equations in eq. 2
    
    [~, indx] = min(bic(:));
    [q_yy, q_xy] = ind2sub(size(bic), indx);
    
    bic = [];
    for q_xx = 1:qmax
        for q_yx = 1:qmax
            [a, b] = arols(x, y, q_xy, q_yy);
            [c, d] = arols(y, x, q_yx, q_xx);
            cov = getcov(x, y, a, b, c, d);
            bic(q_xx, q_yx) = gbic(N, cov, [q_xx, q_xy, q_yx, q_yy]);
    end
    
    [~, indx] = min(bic(:));
    [q_xx, q_yx] = ind2sub(size(bic), indx);
    
    end
    