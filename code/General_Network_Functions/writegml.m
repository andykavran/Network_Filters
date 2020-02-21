function writegml(edgelist, metadata, filename)
%WRITEGML creates a .gml file for a node annotated graph
%   WRITEGML takes in an edge list and a metadata vector to generate a .gml
%   file. Data in the metadata vector is written to the 'label' field to
%   the corresponding node in edge list. Nodes in edge list should be
%   numbered according to their position in the metadata vector.

num_nodes = length(metadata);
num_edges = size(edgelist,2);
fileID = fopen(filename, 'w');
%write nodes
fprintf(fileID, 'graph\n[\n');
for jj = 0:(num_nodes-1)
    fprintf(fileID, ' node\n [\n');
    fprintf(fileID, '  id %d\n  label "%d"\n', jj, metadata(jj+1));
    fprintf(fileID, ' ]\n');
end

for kk = 1:num_edges
    fprintf(fileID, ' edge\n [\n');
    fprintf(fileID, '  source %d\n  target %d\n', edgelist(1,kk)-1, edgelist(2,kk)-1);
    fprintf(fileID, ' ]\n');
end
fprintf(fileID,']\n');
fclose(fileID);
end

