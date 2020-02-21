function [neighbors, weights] = getKNeighbors(train, test,k)
% A function that give the indicies and weights of the K nearest neighbors
% of the test data.
% 
% The training and test should be matricies of n rows of observations of m
% columns of features. K should be an integer specifiying how many
% neighbors should be returned.
%
% Created by Andy Kavran March 2017
    [test_row, test_col] = size(test);
    [train_row, train_col] = size(train);
    if(test_col~=train_col)
        error('The training and test set do not have the same number of features')
    end
    if(sum(isnan(train)>0))
        error('This function cannot handle NaNs in the training data')
    end
    if(sum(isnan(test)>0))
        error('This function cannot handle NaNs in the test data')
    end
    neighbors = zeros(test_row,k);
    weights = zeros(test_row,k);
    for ii = 1:test_row
      dist = zeros(1,train_row); 
      V = ones(train_row,1) * test;
      D = train - V;
      for jj = 1:train_row
          dist(jj) = ( D(jj,:) * D(jj,:)' ) .^ (1/2) ;
      end
      [dist_sort, I] = sort(dist);
      sum_dist = sum(dist_sort(1:k));
      similarities(ii,:)=(sum_dist./dist_sort(1:k));
      sum_metric=sum(similarities(ii,:));
      weights(ii,:)=similarities./sum_metric;
      neighbors(ii,:)=I(1:k);
    end
end