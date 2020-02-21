function writeWG2(A, metadata, filename)
%writeWG2 
%   
    fileID = fopen(filename, 'w');
    degrees = full(sum(A));
    num_nodes = size(A,1);
    for ii = 1:num_nodes
        fprintf(fileID, '%d node, %f %d [ ', ii-1, metadata(ii), degrees(ii));
        neighbors = find(A(ii,:));
        for jj = 1:length(neighbors)
            fprintf(fileID, '(%d, 1) ', neighbors(jj)-1);
        end
        fprintf(fileID, ']\n');
    end
end

