function xt = pagetranspose(x)
    %% Replicate new matlab versions function pagetranspose to not mess with 
    % already written scripts. Only works for 3d arrays.
    
    xt = permute(x, [2 1 3]);
    end