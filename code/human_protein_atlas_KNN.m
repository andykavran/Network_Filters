%% Healthly --> Cancerous Protein Expression KNN Regression 
% 
addpath(genpath('.'))

%import and prep the data
edges = csvread('../data/edge_list.csv', 1, 0);
cancer = csvread('../data/cancer_data.csv', 1, 0)'; % note the transpose
tissue = csvread('../data/healthy_data.csv', 1, 0)'; % note the transpose
fileID = fopen('../data/cancerNames.txt', 'r');
cancer_names = textscan(fileID, '%s\n');
fclose(fileID);
cancer_names = cancer_names{:}; % names of cancers/healthy tissues
no_metadata_nodes_one_index = (cancer(1,end)+2):max(max(edges)+1);
communities = dlmread('../data/FoundComms1.tsv');
communities(no_metadata_nodes_one_index,:) = [];
comms = communities(:,2)+1;
%%

% edge list is zero indexed, so add 1 to everything
edges = edges +1;
A = makeAdjMat(edges);
proteins = size(A,1);
for kk = 1:proteins
    A(kk,kk)=0;
end
G = graph(A);

%remove nodes that don't have metadata. metadata files zero-indexed

G2 = rmnode(G, no_metadata_nodes_one_index);
%find connected component
bins = conncomp(G2);
nbins= max(bins);
N=histcounts(bins, nbins);
[~, I]=max(N);
rm_nodes=find(bins~=I);
keep_nodes=find(bins==I);
G2 = rmnode(G2,rm_nodes);
cancer2 = cancer(2:end,keep_nodes); % cancer metadata
tissue2 = tissue(2:end,keep_nodes); % healthy tissue metadata
comms = comms(keep_nodes);

%remove missing metadata and corresponding nodes
keep_nodes = find((isnan(tissue2(1,:)) == 0));
rm_nodes = find((isnan(tissue2(1,:)) == 1));
G2 = rmnode(G2,rm_nodes);
cancer2 = cancer2(:,keep_nodes);
tissue2 = tissue2(:,keep_nodes);
comms = comms(keep_nodes);

%find connected component again
bins = conncomp(G2);
nbins= max(bins);
N=histcounts(bins, nbins);
[~, I]=max(N);
rm_nodes=find(bins~=I);
keep_nodes=find(bins==I);
G2 = rmnode(G2,rm_nodes);
cancer2 = cancer2(:,keep_nodes);
tissue2 = tissue2(:,keep_nodes);
comms = comms(keep_nodes);

edge_list = G2.Edges.EndNodes;
[edge_list(:,1), I] = sort(edge_list(:,1));
edge_list(:,2) = edge_list(I,2); % PPIN edge_list
A = makeAdjMat(edge_list); % PPIN adjacency matrix

%% Run network filters
% patchwork mean-only filter
num_coms = max(comms);
for kk= 1:20 % for each healthy tissue and cancer
    % Run network filter 
    for ii=1:num_coms % cancer
        this_community = find(comms==ii);
        within_adjacency = A(this_community, this_community);
        cancer_metadata = cancer2(kk, this_community);
        cancer_patchwork_filter_mean_only(kk, this_community) = meanNetworkfilter(within_adjacency, cancer_metadata);
    end
    
    for ii = 1:num_coms % healthy
        this_community = find(comms==ii);
        within_adjacency = A(this_community, this_community);
        healthy_metadata = tissue2(kk, this_community);
        healthy_patchwork_filter_mean_only(kk, this_community) = meanNetworkfilter(within_adjacency,healthy_metadata); 
    end
end

%typical patchwork filter
num_coms = max(comms);
for kk= 1:20 % for each healthy tissue and cancer
    % Run network filter 
    for ii=1:num_coms % cancer
        this_community = find(comms==ii);
        within_adjacency = A(this_community, this_community);
        cancer_metadata = cancer2(kk, this_community);
        assort_coef = metadassort(cancer_metadata, within_adjacency);
        if assort_coef < 0
            cancer_patchwork_filter(kk, this_community) = 0.8*(cancer_metadata - meanNetworkfilter(within_adjacency, cancer_metadata)) + mean(cancer_metadata);
        else
            cancer_patchwork_filter(kk, this_community) = meanNetworkfilter(within_adjacency, cancer_metadata);
        end
    end
    
    for ii = 1:num_coms % healthy
        this_community = find(comms==ii);
        within_adjacency = A(this_community, this_community);
        healthy_metadata = tissue2(kk, this_community);
        assort_coef = metadassort(healthy_metadata, within_adjacency);
        if assort_coef < 0
            healthy_patchwork_filter(kk, this_community) = 0.8*(healthy_metadata - meanNetworkfilter(within_adjacency, healthy_metadata)) + mean(healthy_metadata);
        else
            healthy_patchwork_filter(kk, this_community) = meanNetworkfilter(within_adjacency, healthy_metadata);
        end 
    end
end

% mean and median filters
for jj=1:20 % for each healthy tissue and cancer 
    cancer_mean(jj,:) = meanNetworkfilter(A, cancer2(jj,:));
    tissue_mean(jj,:) = meanNetworkfilter(A, tissue2(jj,:));
    cancer_median(jj,:) = medianNetworkfilter(edge_list, cancer2(jj,:));
    tissue_median(jj,:) = medianNetworkfilter(edge_list, tissue2(jj,:));
    
    cancer_sharp(jj,:) = 0.8*(cancer2(jj,:) - cancer_mean(jj,:)) + mean(cancer2(jj,:));
    tissue_sharp(jj,:) = 0.8*(tissue2(jj,:) - tissue_mean(jj,:)) + mean(tissue2(jj,:));
end

%% KNN Regression
K_min = 1;  % smallest choice of K neighbors
K_max=6; % largest choice of K neighbors
dims = 4; % how many PC's to use

% Delta Vectors--No Filter
difference = cancer2 - tissue2;
mae_no_filter = healthy_cancer_tissue_regression(tissue2, difference, K_min, K_max, dims);

% Delta Vectors--Mean Filter
difference = cancer_mean - tissue_mean;
mae_mean = healthy_cancer_tissue_regression(tissue_mean, difference, K_min, K_max, dims);

% Delta Vectors--Median Filter
difference = cancer_median - tissue_median;
mae_median = healthy_cancer_tissue_regression(tissue_median, difference, K_min, K_max, dims);

% Delta Vectors--Sharp Filter
difference = cancer_sharp - tissue_sharp;
mae_sharp = healthy_cancer_tissue_regression(tissue_sharp, difference, K_min, K_max, dims);

% Delta Vectors--Patchwork-- Mean Filter only
difference = cancer_patchwork_filter_mean_only - healthy_patchwork_filter_mean_only;
mae_patchwork_mean_only = healthy_cancer_tissue_regression(healthy_patchwork_filter_mean_only, difference, K_min, K_max, dims);

% Delta Vectors--Patchwork Filter -- Smooth and Sharp
difference = cancer_patchwork_filter - healthy_patchwork_filter;
mae_patchwork = healthy_cancer_tissue_regression(healthy_patchwork_filter, difference, K_min, K_max, dims);


%% Plot the results
%bootstrap confidence intervals
boot_reps = 1000; % number of samples to bootstrap
alpha_level = 0.01; % 99% confidence
no_filter_ci = bootci(boot_reps, {@mean, mae_no_filter}, 'alpha', alpha_level); 
mean_ci = bootci(boot_reps, {@mean, mae_mean}, 'alpha', alpha_level); 
median_ci = bootci(boot_reps, {@mean, mae_median}, 'alpha', alpha_level); 
patchwork_mean_ci = bootci(boot_reps, {@mean, mae_patchwork_mean_only}, 'alpha', alpha_level); 
sharp_ci = bootci(boot_reps, {@mean, mae_sharp}, 'alpha', alpha_level); 
patchwork_ci = bootci(boot_reps, {@mean, mae_patchwork}, 'alpha', alpha_level); 

fig4 = figure;
ax = gca;
xlabel('# of K Nearest Neighbors')
ylabel('MAE')
ax.FontSize = 14;
hold on
patch([K_min:K_max, fliplr(K_min:K_max)], [no_filter_ci(1,:),   fliplr(no_filter_ci(2,:))], 'black', 'FaceAlpha', .3);
patch([K_min:K_max, fliplr(K_min:K_max)], [mean_ci(1,:), fliplr(mean_ci(2,:))], 'red', 'FaceAlpha', .3);
patch([K_min:K_max, fliplr(K_min:K_max)], [median_ci(1,:), fliplr(median_ci(2,:))], [0, 0.3922, 0], 'FaceAlpha', .3);
patch([K_min:K_max, fliplr(K_min:K_max)], [patchwork_mean_ci(1,:), fliplr(patchwork_mean_ci(2,:))],'cyan', 'FaceAlpha', .3);
patch([K_min:K_max, fliplr(K_min:K_max)], [patchwork_ci(1,:), fliplr(patchwork_ci(2,:))], 'magenta', 'FaceAlpha', .3);
patch([K_min:K_max, fliplr(K_min:K_max)], [sharp_ci(1,:), fliplr(sharp_ci(2,:))], 'blue', 'FaceAlpha', .3);

aa= plot(K_min:K_max, mean(mae_no_filter(:,1:K_max)), 'kx:');
bb= plot(K_min:K_max, mean(mae_mean(:,1:K_max)), 'ro:');
cc= plot(K_min:K_max, mean(mae_median(:,1:K_max)), '+:', 'Color', [0, 0.3922, 0]);
dd= plot(K_min:K_max, mean(mae_patchwork_mean_only(:,1:K_max)), 'c*:');
ff= plot(K_min:K_max, mean(mae_patchwork(:,1:K_max)), 'm*:');
ee= plot(K_min:K_max, mean(mae_sharp(:,1:K_max)), 'b*:');
[l, legendIcons] = legend([aa, bb, cc, dd, ff, ee], {'No Filter', 'Mean Filter', 'Median Filter', 'Patchwork Filter (mean-only)', 'Patchwork Filter',...
    'Sharp Filter'}, 'Location', 'best');
l.FontSize = 14;
ax.YLim(2) = .65;
savefig(fig4, '../results/Protein_Atlas_Output/KNN_Regression.fig', 'compact')
plot2svg('../results/Protein_Atlas_Output/KNN_Regression.svg', fig4, '', legendIcons)
%% Histogram of assortativity coefficients
for kk = 1:20
    for ii=1:num_coms % cancer
        this_community = find(comms==ii);
        within_adjacency = A(this_community, this_community);
        cancer_metadata = cancer2(kk, this_community);
        assort_coefs_cancer(ii,kk) = metadassort(cancer_metadata, within_adjacency);
    end
end

for kk = 1:20
    for ii=1:num_coms % cancer
        this_community = find(comms==ii);
        within_adjacency = A(this_community, this_community);
        healthy_metadata = cancer2(kk, this_community);
        assort_coefs_healthy(ii,kk) = metadassort(healthy_metadata, within_adjacency);
    end
end
fig5 = figure;
ax = gca;
xlabel('Assortativity Coefficients')
ylabel('Count')
hold on
histogram([assort_coefs_cancer(:), assort_coefs_healthy(:)])
savefig(fig5, '../results/Protein_Atlas_Output/HPA_assortativity_coeffs.fig', 'compact')