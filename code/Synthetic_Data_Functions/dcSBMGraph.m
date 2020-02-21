function [ edge_list, B, bins, g2 ] = dcSBMGraph(N, w, partition, G, g1 )
%DCSBMGRAPH Summary of this function goes here
%   Detailed explanation goes here

    %% Define Graph Parameters
    % generate power law distributed degrees
    invalid_graph = true;
    while(invalid_graph==true)
        theta = round(randht(N, 'powerlaw', 3));
        invalid_graph = max(theta)>20;
    end
    
    theta_norm = zeros(size(theta));
    % normalize theta
    for ii=1:G
        index = find(g1==ii);
        theta_norm(index) = theta(index)/sum(theta(index));
    end

    % sum of edges in each partition, kappa
    k = sum(w,2);
    % number of edges in graph
    m = sum(k)/2; %half the sum of the adjacency matrix

    % define omega random according to degree corrected SBM. wr means w random
    %
    wr = zeros(G);
    for ii=1:G
        for jj=1:G
            wr(ii,jj) =( k(ii) * k(jj) ) ./(2*m);
        end
    end
    wr = round(wr);

    % define L as lambda, the mixing coefficient of planted and random
    % partitions. The greater lambdais , the greater the planted partition 
    % model rather than random mixing.
    L = 0.85;
    % mix random graph with planted partition structure. linear combination
    combination = L * w + (1-L) * wr;
    for r = 1:(G)
        for s = (r):G
            omega(r,s) = poissrnd(combination(r,s));
        end
    end
    %% Generate network
    % uses a stub matching algorithm
    % draw edge from a vertex
    edgeA = [];
    edgeB = [];
    for jj = 1:G % for each community
        for kk = (jj):G %to other communities not already seen (symmetric)
            R=[]; 
            S=[];
            stubA = cumsum(theta_norm(partition{jj}));
            probs = rand(omega(jj,kk),1); % rand # for each edge
            for mm = 1:length(probs) % assign each edge to a vertex
                R(mm) = find(stubA>probs(mm),1, 'first');
            end
            edgeA = [edgeA, partition{jj}(R)];

            stubB = cumsum(theta_norm(partition{kk}));
            probs = rand(omega(jj,kk),1);% rand # for each edge
            for mm = 1:length(probs) % assign each edge to a vertex
                S(mm) = find(stubB>probs(mm), 1, 'first');
            end
            edgeB = [edgeB, partition{kk}(S)];
        end
        A = sparse(edgeA, edgeB,ones(length(edgeA),1),N,N); % match stubs
    end
    A = A + A'; % the matrix is symmetric
    for ii = 1:length(A) 
        A(ii,ii) = 0;% remove self edges
    end
    A(A>0) = 1; % remove multi edges
    g = graph(full(A)); % create graph object
    cc = conncomp(g); % find all components
    temp = mode(cc); % find the largest component
    keep = find(cc==temp); % keep the largest component
    g2 = g1(keep); % keep the largest component
    % count number of nodes remaining in each community
    [bins, ~] = histcounts(g2,G); 
    % keep largest connectd component in new adjacency matrix
    B = A(keep,keep);
    % define edge list from B
    edge_list = sparse2edgelist(B);
end

