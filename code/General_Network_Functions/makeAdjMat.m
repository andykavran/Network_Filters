function A = makeAdjMat(edgelist)
    
    if(size(edgelist,1)~=2)
        edgelist = edgelist';
    end
    nEdges=size(edgelist,2);
    nNodes=max(max(edgelist));
    A=sparse(nNodes, nNodes);
    for ii =1:(nEdges-1)
        s = edgelist(1,ii);
        t = edgelist(2,ii);
        A(s,t)=A(s,t)+1;
        A(t,s)=A(t,s)+1;
    end
end