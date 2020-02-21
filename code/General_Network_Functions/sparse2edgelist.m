function [edge_list] = sparse2edgelist(A)
%SPARSE2EDGELIST converts a sparse adjacency matrix into an edge list
%   edge_list = sparse2edgelist(A);
%   
%   The adjacency matrix, A, should be an NxN symmetric sparse matrix. A
%   2xM edge list is returned. The row/column index of A are taken to be the
%   node IDs for the edge list.

% Version 1.0 (July 2017)
% Written by Andy Kavran
[row, col] = find(A);
duplicate = row > col;
edge_list = [row';col'];
edge_list(:,duplicate)=[];
end