function A = chungLuModel(W)
%CHUNGLUMODEL generates a simple random graph with a given expected degree distribution 
% 
% A function to efficiently create a Chung Lu Random Graph.
% Input W is a sorted 1xN vector of  decreasing weights (expected degrees)
% for the corresponding node. Output is an Mx2 edge list.
% This uses an O(N+M) algorithm described in 
% J.C Miller and A. Hagberg, Efficient Generation of Networks with Given
% Expected Degrees, Algorithms and Models for the Web Graph, 2011/5/27
%
% Syntax:
% A = chungLuModel(W);
%
% Example:
% W = randht(100, 'powerlaw', 3); % create powerlaw degree distribution
% W = sort(W, 'descend'); % sort in descending order
% A = chungLuModel(W); % get adjacency matrix
%
% Written by Andy Kavran, July 2017
s = sum(W);
temp = size(W);
if temp(1) > temp(2) % if column vector change to row vector
    W = W';
end
N = size(W,2); % Number of Nodes
A = sparse(N, N);
%u and v are node indices. For each u, determine which v will be neighbors
for u = 1:(N-1)
   v = u+1;
   p = min((W(u)*W(v))/s, 1); %don't want p greater than 1
  while (v<N+1 && p>0)
       if p ~= 1
           r = rand;
           v = v + floor(log(r)/log(1-p));
       end
       if v < N
           q = min((W(u)*W(v))/s, 1); %don't want q greater than 1
           r = rand;
           if r < q/p
               % add edge between u and v
               A(u,v) = 1;
               A(v,u) = 1;
           end
            p = q;
            v = v + 1;
       end
   end
end
