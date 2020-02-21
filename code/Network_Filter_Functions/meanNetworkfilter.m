function mean_data = meanNetworkfilter(A, x)
%MEANNETWORKFILTER
%   Explain the function usage
    
    % Check data input
    
    % Is x a column vector?
    if(size(x,1) > 1 && size(x, 2) > 1)
        error('Metadata vector not correct dimensions')
    elseif(size(x,1) ==1 && size(x,2) >1)
        % transpose x into a column vector
        x = x';
    end
    % Is A a square matrix?
    if(size(A,1) ~= size(A,2))
        error('Adjacency Matrix must be square')
    end
    % Are A and X are same size?
    if(size(A,1) ~= size(x,1))
        error('Metadata and Network are different sizes')
    else
            num_nodes = size(A, 1);
    end
    
    % Calculate Mean Network Filter
    S = A + eye(num_nodes); % add 1 to diagonal using identity matrix
    K = 1./sum(S, 2); % calculate degrees
    K = K * ones(1, num_nodes);
    mean_data = S.*K*x; % calculate mean filter
    mean_data = mean_data';
end
