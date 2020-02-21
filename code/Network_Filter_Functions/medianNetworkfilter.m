function median_data = medianNetworkfilter(edge_list, x)
%MEDIANNETWORKFILTER
%   x is a 
    
    % Prepare adjacency list
    if (size(x,1) ~=1 && size(x,2) ==1)
        x = x';
    end
    neighbors = zeros(1,500); %initialize a big neighbor list. 500 is big enough for PPIN
    num_nodes=length(x);
    median_data = zeros(1,num_nodes);
    % prepare search algorithm
    [nrow, ncol]=size(edge_list);
    if(isempty(edge_list))
        warning('The edge list is empty. Returning uncorrected metadata')
        median_data = x;
        return;
    elseif(nrow ~= 2 && ncol ~=2)
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
        self_meta = x(ii);
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
            median_data(ii) = x(ii); %keep same value
            continue;
        end
        try
            neigh_metadata = x(neighbors(1:count));
        catch
            disp(neighbors(1:count))
            error('Houston, we have a problem with ii = %d and count = %d', ii, count)
        end
        median_data(ii) = median([neigh_metadata, self_meta]);
    end 
end
