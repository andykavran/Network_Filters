function [meanabserr, y_values] = healthy_cancer_tissue_regression(X_data, Y_data, K_min, K_max, dims)
% [meanabserr, y_values] = healthy_cancer_tissue_regression(X_data, Y_data, K_min, K_max, dims)
% main function for kNN healthy tissue - cancer regression task
% leave one out cross validation

    k_runs = K_max - K_min + 1;
    meanabserr = zeros(size(X_data,1), k_runs); % holds rmse values
    y_values = zeros(size(X_data,1), K_max, dims); % hold predicted data
    for ii=1:size(X_data, 1) %leave one observation out each interation
        index = 1:(size(X_data,1)); %create a list that will index the matrix
        index(ii)=[]; % leave index ii out for training
        raw_x = X_data(index,:);%training predictor points
        raw_y = Y_data(index,:);%training response points
        [mapx, x_features, ~] = pca(raw_x, 'NumComponents', dims, 'Centered', false);
        test_x_features =  X_data(ii,:)*mapx;
        [mapy, y_features, ~] = pca(raw_y, 'NumComponents', dims, 'Centered', false);
        count = 1;
        for K = K_min:K_max % use different number of neighbors
            [neighbors, weights] = getKNeighbors(x_features, test_x_features, K);
            fitted_test_y = kNNRegression(y_features(neighbors,:), weights);
            original_space_y = fitted_test_y * mapy';
            meanabserr(ii,count) = maerr(original_space_y, Y_data(ii,:));
            y_values(ii,count,:) = fitted_test_y;
            count = count + 1;
        end
    end
end