%% Synthetic Graphs with community structure
% Test network filters on synthetic graphs with community structure
% Code takes about half an hour for 100 reps.

addpath(genpath('/Users/Andy/Google Drive/Grad.Research/Manuscripts/Actual Network Filter Paper/github_code/code'))
N = 500; %nodes in a network
G = 5; % number of partitions
multiplier=0.8; % sharp filter multiplier
reps = 100; % number of graphs to generate
%%
%instatiate vectors to hold errors
no_filter_error = zeros(1, reps);
patchwork_error = no_filter_error;
mean_filter_mae = no_filter_error;
sharp_filter_mae = no_filter_error;
median_filter_error = no_filter_error;
community_medain_permute_error = no_filter_error;

% assign 100 nodes to each community
g1 = [repelem(1,100), repelem(2,100), repelem(3, 100), repelem(4,100), repelem(5,100)];

for ii = 1:G
    partition{ii} = find(g1==ii);
end

% SBM structure matrix, omega planted. call it w.
w(1,1) = 250; 
w(2,2) = 250;
w(3,3) = 250;
w(4,4) = 250;
w(5,5) = 250;

% range of assortativity coefficients
upper_r = .7;
lower_r = .4;
% means of each community metadata distribution
mu = [110, 80, 60, 40, 20];
sigma = 5;
cd('Newman_Clauset_Community_Detection')
%%
tic
for loop = 1:reps % create graph
    %print out loop number 
    if(mod(loop,floor(0.1*reps))==0 && loop ~= 1)
        toc
        fprintf('Loop %d\n',loop)
        tic
        fclose('all');
    end
    % get dcSBM generated graph
    [edge_list, B, bins, g2] = dcSBMGraph(N, w, partition, G, g1 );
    assort_coeffs = (upper_r - lower_r) .* rand(1, G) + lower_r; % get assortativity coeffs
    for negatives = 0:5
        % generate metadata
        disassort = randi(length(assort_coeffs), negatives, 1); % pick which comms are disassort
        ac2 = assort_coeffs;
        ac2(disassort) = ac2(disassort) .* -1; % make these negative
        x_true = zeros(length(B), 1);
        x_cell = mat2cell(x_true, bins);
        x_cell_err = x_cell;
        bins_2 = cumsum([0, bins]);
        perm_save = [];
        for ii = 1:length(mu) % for each community
            %pd = makedist('Poisson', mu(ii));
            pd = makedist('normal', 'mu', mu(ii), 'sigma', sigma);
            % create ground truth metadata
            this_community = find(g2 == ii);
            num_comm_nodes = length(this_community);
            [x_cell{ii}, r_hat] = assortMCMC(B(this_community, this_community), ac2(ii), pd);
            x_cell{ii} = x_cell{ii}';
            
            % randomly permute a fraction of the data to create noise
            % first pick the nodes
            perm = randsample(num_comm_nodes, floor(.25 * num_comm_nodes));
            order = randsample(length(perm), length(perm)); % reorder them
            x_cell_err{ii} = x_cell{ii};
            temp = x_cell{ii}(perm); % hold the values to be permuted
            temp = temp(order); % permute them
            x_cell_err{ii}(perm) = temp; % reassign the values
            perm_save = [perm_save; bins_2(ii)+perm];% save which nodes were permuted
        end
        x_true = cell2mat(x_cell);
        x_error = cell2mat(x_cell_err); 
        % community detection using c code
        output = 1;
        norm_metadata = (x_error - min(x_error)+1) ./(max(x_error)+1 - min(x_error));
        metadata_out = horzcat((0:length(x_error)-1)', norm_metadata);
        dlmwrite('synthetic_metadata.txt', metadata_out, '\t');
        writeWG2(B, x_error, 'dcSBM_synthetic_data.wg2')
        restarts = 10;
        logLikelihood = zeros(1,restarts);
        for rand_restart = 1:restarts
            output = 1;
            while(output ~= 0)
                out_file = sprintf('msbm_output_%d.txt', rand_restart);
                err_file = sprintf('std_err_%d.txt', rand_restart);
                command = sprintf('./metadata dcSBM_synthetic_data.wg2 synthetic_metadata.txt >%s 2>%s', out_file, err_file);
                output = system(command);
                stream = fileread(err_file);
                LL = regexp(stream,'Log-likelihood = ([-+]?\d*\.?\d+)', 'tokens');
                if(isempty(LL)~=true)
                    logLikelihood(rand_restart) = str2double(char(LL{end}));
                else
                    logLikelihood(rand_restart) = -Inf;
                end
                
            end
        end
        [~, ind] = max(logLikelihood);
        maxLL = sprintf('msbm_output_%d.txt', ind);
        comm_file = dlmread(maxLL);
        % get communities
        group = zeros(1,length(comm_file));
        for ii = 1:length(comm_file)
            [~, group(ii)] = max(comm_file(ii,3:end));
        end
        % apply network filters on x_err
        
        patchwork_filter_data = zeros(size(x_error));
        
        for ii = 1:max(group)
            this_community = find(group==ii);
            within_adjacency = B(this_community, this_community);
            this_metadata = x_error(this_community);
            this_r = metadassort(this_metadata, within_adjacency);
            community_edges = sparse2edgelist(within_adjacency);
            if this_r >=0
                patchwork_filter_data(this_community) = meanNetworkfilter(within_adjacency, this_metadata);
            else
                patchwork_filter_data(this_community) = multiplier.*(this_metadata' - meanNetworkfilter(within_adjacency, this_metadata)) + mean(this_metadata);
            end
        end
        mean_filter_data = meanNetworkfilter(B, x_error);
        sharp_filter_data = multiplier.*(x_error' - mean_filter_data) + mean(x_error);
        median_filter_data = medianNetworkfilter(edge_list, x_error');

        %calculate error
        
        %the vectors sometimes change shape for some reason that I haven't
        %figured out yet.
        if(size(patchwork_filter_data)~=size(x_true))
            patchwork_filter_data = patchwork_filter_data';
        end
        if(size(mean_filter_data)~=size(x_true))
            mean_filter_data = mean_filter_data';
        end
        if(size(sharp_filter_data)~=size(x_true))
           sharp_filter_data = sharp_filter_data'; 
        end
        if(size(median_filter_data)~=size(x_true))
           median_filter_data = median_filter_data'; 
        end
        
        no_filter_error(negatives+1, loop) = maerr(x_true(perm_save), x_error(perm_save));
        patchwork_error(negatives+1, loop)= maerr(x_true(perm_save), patchwork_filter_data(perm_save));
        mean_filter_mae(negatives+1, loop)= maerr(x_true(perm_save), mean_filter_data(perm_save));
        sharp_filter_mae(negatives+1, loop) = maerr(x_true(perm_save), sharp_filter_data(perm_save));
        median_filter_error(negatives+1, loop)= maerr(x_true(perm_save), median_filter_data(perm_save));
    end
end
cd ..
save('../results/Community_Random_Graphs_Output/workspace.mat')
%% MAE Plot

mean_md_permute = flipud(mean((no_filter_error),2));
mean_community_permute = flipud(mean((patchwork_error),2));
mean_smooth_permute = flipud(mean((mean_filter_mae),2));
mean_sharp_permute = flipud(mean((sharp_filter_mae),2));
mean_median_permute = flipud(mean((median_filter_error),2));

boot_reps = 200;
md_permute_ci        = bootci(boot_reps, {@mean, no_filter_error'}, 'alpha', .01); 
community_permute_ci = bootci(boot_reps, {@mean, patchwork_error'}, 'alpha', .01); 
mean_permute_ci      = bootci(boot_reps, {@mean, mean_filter_mae'}, 'alpha', .01); 
median_permute_ci    = bootci(boot_reps, {@mean, median_filter_error'}, 'alpha', .01); 
sharp_permute_ci     = bootci(boot_reps, {@mean, sharp_filter_mae'}, 'alpha', .01); 

fig3 = figure;
fig3.Renderer = 'painters';
hold on
patch([fliplr(0:G), 0:G], [md_permute_ci(1,:),   fliplr(md_permute_ci(2,:))], 'black', 'FaceAlpha', .3);
patch([fliplr(0:G), 0:G], [mean_permute_ci(1,:), fliplr(mean_permute_ci(2,:))], 'red', 'FaceAlpha', .3);
patch([fliplr(0:G), 0:G], [median_permute_ci(1,:), fliplr(median_permute_ci(2,:))], [0,.3922,0], 'FaceAlpha', .3);
patch([fliplr(0:G), 0:G], [sharp_permute_ci(1,:), fliplr(sharp_permute_ci(2,:))], 'blue', 'FaceAlpha', .3);
patch([fliplr(0:G), 0:G], [community_permute_ci(1,:), fliplr(community_permute_ci(2,:))], [.6627, .0, .6627], 'FaceAlpha', .3);

marker_size = 10;
l1 = plot(0:G,mean_md_permute,'kx:', 'MarkerSize', marker_size); % no filter
l2 = plot(0:G,mean_smooth_permute, 'ro:','MarkerSize', marker_size); % mean filter
l3 = plot(0:G,mean_median_permute, '+:', 'MarkerSize', marker_size, 'Color', [0,.3922,0]); % median filter
l4 = plot(0:G,mean_sharp_permute, 'b*:', 'MarkerSize', marker_size); % sharp filter
l5 = plot(0:G,mean_community_permute, 's:', 'MarkerSize', marker_size, 'Color', [.6627, .0, .6627]); % patchwork

xlabel('Fraction of communities that have assortative values', 'FontSize', 16)
ylabel('Mean Absolute Error', 'FontSize', 16)
[~, legendIcons] = legend([l1, l2, l3, l4, l5], {'No Filter', 'Mean Filter',...
    'Median Filter', 'Sharp Filter','Patchwork Filter'},'Location', 'best');

ax = gca;
ax.FontSize = 14;
ax.YLim = [0, 27];
ax.XTick = [0, 1, 2, 3, 4, 5];
ax.XLim = [-0.1, 5.1];
ax.XTickLabel = {'0%', '20%', '40%', '60%', '80%', '100%'};

title(['Reps: ', int2str(reps)]);
savefig(fig3, '../results/Community_Random_Graphs_Output/permuted_nodes_mean_absolute_error.fig', 'compact')
plot2svg('../results/Community_Random_Graphs_Output/permuted_nodes_mean_absolute_error.svg', fig3, '', legendIcons)
%% Zoomed in MAE Plot
fig3 = figure;
fig3.Renderer = 'painters';
hold on
patch([fliplr(0:G), 0:G], [md_permute_ci(1,:),   fliplr(md_permute_ci(2,:))], 'black', 'FaceAlpha', .3);
patch([fliplr(0:G), 0:G], [mean_permute_ci(1,:), fliplr(mean_permute_ci(2,:))], 'red', 'FaceAlpha', .3);
patch([fliplr(0:G), 0:G], [median_permute_ci(1,:), fliplr(median_permute_ci(2,:))], [0,.3922,0], 'FaceAlpha', .3);
patch([fliplr(0:G), 0:G], [community_permute_ci(1,:), fliplr(community_permute_ci(2,:))],'cyan', 'FaceAlpha', .3);

l1 = plot(0:G,mean_md_permute,'kx:', 'MarkerSize', marker_size); % no filter
l2 = plot(0:G,mean_smooth_permute, 'ro:','MarkerSize', marker_size); % mean filter
l3 = plot(0:G,mean_median_permute, '+:', 'MarkerSize', marker_size, 'Color', [0,.3922,0]); % median filter
l5 = plot(0:G,mean_community_permute, 'cs:', 'MarkerSize', marker_size); % patchwork

xlabel('% of communities that have assortative values', 'FontSize', 14)
ylabel('Mean Absolute Error', 'FontSize', 14)
[~, legendIcons] = legend([l1, l2, l3, l5], {'No Filter', 'Mean Filter',...
    'Median Filter','Patchwork Filter'},'Location', 'best');
ax = gca;
ax.YLim = [3, 7];
ax.XTick = [0, 1, 2, 3, 4, 5];
ax.XTick = [0, 1, 2, 3, 4, 5];
ax.XLim = [0, 5];
ax.XTickLabel = {'0%', '20%', '40%', '60%', '80%', '100%'};
title(['Zoom in. Reps: ', int2str(reps)], 'FontSize', 14);
savefig(fig3, '../results/Community_Random_Graphs_Output/zoomed_permuted_nodes_mean_absolute_error.fig', 'compact')
plot2svg('../results/Community_Random_Graphs_Output/zoomed_permuted_nodes_mean_absolute_error.svg', fig3, '', legendIcons)
