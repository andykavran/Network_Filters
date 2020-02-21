function [x, r_hat, counter] = assortMCMC(A, r, varargin)
% metadata = assortMCMC(edge_list, assort_coef, varargin)
% 
% A function to create metadata on a network defined by the sparse matrix A
% with the assortativity coefficient given by r. A user defined
% distribution object created with the matlab function makedist() can be
% passed as an optional argument. Default metadata distribution is
% Normal(0,1).
%
% Written by Andy Kavran July 2017
% Edited November 3, 2017 to vectorize recalc assortativity code.
if (abs(r)>1)
    error('[assortMCMC] Invalid imput for Assortativity Coefficient r. Must be in range [-1, 1]')
end
if (isempty(varargin)==false)
    ii = 1;
    while ii <= length(varargin)
        if isa(varargin{ii}, 'prob.ProbabilityDistribution')
           pd = varargin{ii};
        elseif strcmp(varagin{ii}, 'verbose')
            verbose = varargin{ii+1};
        else
           warning('[assortMCMC] Ignoring invalid argument %d', ii) 
        end
        ii = ii+1;
    end
else
pd = makedist('Normal', 'mu', 0, 'sigma', 1);
end
counter = 0;
N= size(A,1);
x = random(pd,1,N); %metadata
% calculate assorativity
K = sum(A);
twom = sum(K);
[r_hat, numerator, denom] = metadassort(x,A);
converge_factor = .009;
while(abs(r_hat - r) > converge_factor)
    counter = counter+1;
    x_star = x;
    i = randi(N);
    x_star(i) = random(pd);
    [r_hat_star, numerator_star, denom_star] = recalc_assortativity(i, x_star);
    if (abs(r-r_hat) > abs(r-r_hat_star))
        x(i) = x_star(i);
        r_hat = r_hat_star;
        numerator = numerator_star;
        denom = denom_star;
    end
    if counter >2000
        warning('R hat did not reach R in the Max allowed steps. Final R value is %0.3f', r_hat)
        return
    end

end

    function [r, numerator2, denom2] = recalc_assortativity(node, xx)
        numerator_subtract = (A(node,:)' - (K(node)*K(:))./twom)'* x(node)*x(:);
        denom_subtract = (K(node)*(node==1:N)' - (K(node)*K(:))/twom)' *x(node) *x(:);
        
        numer_add = (A(node,:)' - (K(node)*K(:))./twom)'* xx(node)*xx(:);
        denom_add = (K(node)*(node==1:N)' - (K(node)*K(:))/twom)' *xx(node) *xx(:);
        
        numerator_subtract = (2*numerator_subtract) - ((A(node,node) - (K(node)*K(node)/twom)) *x(node) *x(node));
        denom_subtract = (2*denom_subtract) - ((K(node)*(node==node) - (K(node)*K(node)/twom)) *x(node) *x(node));
        numer_add = (2*numer_add) - ((A(node,node) - (K(node)*K(node)/twom)) *xx(node) *xx(node));
        denom_add = (2*denom_add) - ((K(node)*(node==node) - (K(node)*K(node)/twom)) *xx(node) *xx(node));
        numerator2 = numerator;
        denom2 = denom;
        numerator2 = numerator2 - numerator_subtract + numer_add;
        denom2 = denom2 - denom_subtract + denom_add;
        r = numerator2/denom2;
    end
end




