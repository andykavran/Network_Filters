function [r, numerator, denom] = metadassort(x, A)
%METADASSORT Calculates the node attribute assortativity of a graph
%   METADASSORT takes in node numerical metadata ROW vector and sparse
%   adjacency matrix A and returns the assortativity coefficient r
%
%   Syntax:
%    r = metadassort(x, A);
%
% Written by Andy Kavran July 2017
% Editted November 3, 2017. Vectorized to speed up calculation.

N = size(A,1);
if(length(x)~=N)
    error('The size of the adjacency matrix does not match the size of the metadata vector')
end
if(size(x,1)~=1)
    x=x';
end
K = sum(A); % calculate degree
twom = sum(K); % calculate twice the number of edges

numerator = full(sum(sum((A - K'*K/twom).*(x'*x))));
denom = full(sum(sum((diag(K) - K'*K/twom).*(x'*x))));
r = numerator/denom; % return r
end

