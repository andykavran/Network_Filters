function network_metadata = network_filter(metadata, edge_list, filter, varargin)
% NETWORK_FILTER Modifies node metadata using smoothen, sharpen, or median filters.
%   network_metadata = network_filter(metadata, edge_list, filter);
%
% Metadata is a numeric 1xN row vector of metadata of the nodes in the
% network. Edge_list is an 2xM edge list. The index of the metadata
% should match the id of the nodes in the edge list  (i.e. metadata(i)
% describes the node of with ID i in edge_list). This only works
% for undirected networks and does not account for edge weight.
%
% Created by Andy Kavran July 2017
ii=1;
while ii <= length(varargin)
    argok = true;
    if (ischar(varargin{ii}) == true)
        switch varargin{ii}
            case 'scale'
                scale = varargin{ii+1};
            otherwise
                argok = false;
        end
    end
    if (argok == false)
        warning('[NETWORK_FILTER] Ignoring invalid argument %d.', ii); 
    end
    ii = ii+1;
end
neighbors = zeros(1,500); %initialize a big neighbor list. 500 is big enough for PPIN
num_nodes=length(metadata);
network_metadata = zeros(1,num_nodes);
% prepare search algorithm
[nrow, ncol]=size(edge_list);
if(nrow ~= 2 && ncol ~=2)
    error('Edge_List has incorrect dimensions. The size of one dimension must equal 2');
end
if(ncol<nrow)
    edge_list = edge_list'; %transpose to make wide and short edge_list
end
edge_listdouble=horzcat(edge_list,edge_list([2,1],:));
twonedges = size(edge_listdouble,2);
if(issorted(edge_listdouble(1,:))==false)
    [sort_edge_list, I] = sort(edge_listdouble(1,:));
    sort_edge_list(2,:) = edge_listdouble(2,I);
    edge_list2 = horzcat(sort_edge_list, [0;0]); %important for search alg
else
    edge_list2 = horzcat(edge_listdouble, [0;0]); %important for search alg
end
%create edge location directory
dict = zeros(1,num_nodes); % make space for dict
prev = 0; % store value of node in previous edge. intialized to zero
for ii = 1:twonedges
    if((prev ~= edge_list2(1,ii)) && (dict(edge_list2(1,ii)) == 0))
        prev = edge_list2(1,ii); % update prev
        dict(edge_list2(1,ii)) = ii; % store location 
    end
end
%%%%%%%%%%%%%% Main Algorithm %%%%%%%%%%%%%%
for ii=1:num_nodes
    self_meta = metadata(ii);
    % Search for Neighbors %
    position = dict(ii);
    count=1;
    if(position~=0)
        while(edge_list2(1,position) == ii)
            neighbors(count) = edge_list2(2,position);
            position = position + 1;
            if(edge_list2(1,position) == ii)
                count = count + 1;
            end
        end
    else
        %fprintf('Node %d looks like a singleton, cannot modify its metadata\n', ii)
        network_metadata(ii) = metadata(ii); %keep same value
        continue;
    end
    try
        neigh_metadata = metadata(neighbors(1:count));
    catch
        disp(neighbors(1:count))
        error('Houston, we have a problem with ii = %d and count = %d', ii, count)
    end
    switch(filter)
        case {'mean', 'Mean', 'smooth','Smooth','smoothen', 'Smoothen'}
            %average all neighbors
            network_metadata(ii) = mean([neigh_metadata, self_meta]);
        case {'sharp', 'Sharp','sharpen', 'Sharpen'}
            multiplier = scale;
            local_mean = mean([neigh_metadata, self_meta]);
            network_metadata(ii) = self_meta + multiplier .* (self_meta- local_mean);
        case {'median', 'Median'}
            network_metadata(ii) = median([neigh_metadata, self_meta]);
        otherwise
            error('Invalid filter input');
    end
end 
end