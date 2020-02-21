%% Noisy Chung-Lu Synthetic Data
%
addpath(genpath('.'))
rep = 5000; % how many graphs should be made
r = -0.8:.1:0.8; % range of assortativity coefficients

% arrays to hold error
no_filter_error = zeros(length(r), length(rep));
sharp_error = no_filter_error;
smooth_error = no_filter_error;
median_error = no_filter_error;
no_filter_r = no_filter_error;
smooth_r = no_filter_error;
sharp_r = no_filter_error;
median_r = no_filter_error;
no_filter_permute_error = no_filter_error;
smooth_permute_error = no_filter_error;
sharp_permute_error = no_filter_error;
median_permute_error = no_filter_error;

tic
for ii = 1:rep
    if (mod(ii,floor(max(rep)*.1))==0)
        disp(ii)
    end
    % create graph
    invalid_graph = true;
    while(invalid_graph==true)
        W = randht(100, 'powerlaw', 3, 'xmin',1); % create degree sequence
        W = sort(W, 'descend');
        invalid_graph = max(W)>17;
    end
    A = chungLuModel(W); % create random graph
    g = graph(full(A));
    cc = conncomp(g); % remove unconnected nodes for this task
    remove = find(cc~=mode(cc));
    g = rmnode(g, remove);
    A = adjacency(g);
    loop_count=1;
    for jj = r
        % create metadata
        pd = makedist('Normal', 'mu', 100, 'sigma', 10);
        [x, rhat, counter] = assortMCMC(A,jj,pd);
        % now we want to permute some values of x to create noise
        perm = randsample(length(x), floor(.25*length(x)));
        order = randsample(length(perm), length(perm));
        temp = x(perm);
        temp = temp(order);
        x_err = x;
        x_err(perm) = temp; % data vector with noise

        % use network filters
        edge_list = sparse2edgelist(A);
        mean_filtered  = meanNetworkfilter(A, x_err);
        sharp_filtered = 0.8.*(x_err - mean_filtered) + mean(x_err);
        median_filterd = medianNetworkfilter(edge_list, x_err);

        % calculate mean absolute errors of real and filtered data
        no_filter_error(loop_count,ii) = maerr(x, x_err);
        smooth_error(loop_count,ii) = maerr(x, mean_filtered);
        sharp_error(loop_count,ii)  = maerr(x, sharp_filtered);
        median_error(loop_count,ii) = maerr(x, median_filterd);

        % calculate assortativity of real and filtered data
        no_filter_r(loop_count,ii) = metadassort(x_err, A);
        smooth_r(loop_count,ii) = metadassort(mean_filtered, A);
        sharp_r(loop_count, ii) = metadassort(sharp_filtered, A);
        median_r(loop_count, ii) = metadassort(median_filterd, A);

        % calculate mean absolute errors of only permuted values
        no_filter_permute_error(loop_count,ii) = maerr(x(perm), x_err(perm));
        smooth_permute_error(loop_count,ii) = maerr(x(perm), mean_filtered(perm));
        sharp_permute_error(loop_count,ii) = maerr(x(perm), sharp_filtered(perm));
        median_permute_error(loop_count,ii) = maerr(x(perm), median_filterd(perm));

        loop_count = loop_count + 1;
    end
end
toc
save('../results/Chung_Lu_Output/workspace.mat')

%% Plot MAE of Permuted Elements
total_permute = mean((no_filter_permute_error),2);
smooth_permute = mean((smooth_permute_error),2);
sharp_permute = mean((sharp_permute_error),2);
median_permute = mean((median_permute_error),2);

boot_reps = 200;

B = permute(no_filter_permute_error, [2,1]);
total_permute_ci = bootci(boot_reps, {@mean, B}, 'alpha', .01); 
B = permute(smooth_permute_error, [2,1,3]);
smooth_permute_ci = bootci(boot_reps, {@mean, B}, 'alpha', .01); 
B = permute(sharp_permute_error, [2,1,3]);
sharp_permute_ci = bootci(boot_reps, {@mean, B}, 'alpha', .01); 
B = permute(median_permute_error, [2,1,3]);
median_permute_ci = bootci(boot_reps, {@mean, B}, 'alpha', .01); 


fig1 = figure(1);
marker_size = 10;
set(fig1, 'rend', 'painters', 'Position', [.25, .25, 480, 480])
hold on
ax = gca; 
ax.YLim = [5, 14];
ax.XLim = [-0.82, 0.82];
patch([r, fliplr(r)], [total_permute_ci(1,:), fliplr(total_permute_ci(2,:))], 'black', 'FaceAlpha', .3);
patch([r, fliplr(r)], [smooth_permute_ci(1,:), fliplr(smooth_permute_ci(2,:))], 'red', 'FaceAlpha', .3);
patch([r, fliplr(r)], [median_permute_ci(1,:), fliplr(median_permute_ci(2,:))], [0, 0.3922,0], 'FaceAlpha', .3);
patch([r, fliplr(r)], [sharp_permute_ci(1,:), fliplr(sharp_permute_ci(2,:))], 'blue', 'FaceAlpha', .3);
l1 = plot(r, total_permute,'kx:', 'MarkerSize', marker_size);
l2 = plot(r, smooth_permute,'ro:', 'MarkerSize', marker_size);
l3 = plot(r, median_permute,'+:', 'Color', [0, 0.3922,0], 'MarkerSize', marker_size);
l4 = plot(r, sharp_permute,'b*:', 'MarkerSize', marker_size);

xlabel('Assortativity Coefficient', 'FontSize', 15)
ylabel('Mean Absolute Error', 'FontSize', 15)
ax.XTick = [-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8];
[l, legendIcons] = legend([l1, l2, l3, l4],{'No Filter', 'Mean Filter', 'Median Filter', 'Sharp Filter'}, 'location', 'best');
title(['Reps: ', int2str(rep)])
savefig(fig1, '../results/Chung_Lu_Output/permuted_nodes_error.fig', 'compact')
plot2svg('../results/Chung_Lu_Output/permuted_nodes_error.svg', fig1, '', legendIcons)