MU1 = [1 2]; % input matrix of means mu
SIGMA1 = [2 0; ...
		  0 1]; % input array of covariances

MU2 = [-3 -5];
SIGMA2 = [1 0; ...
	      0 1];

% generate data from a mixture of two bivariate Gaussian distributions
% using the mvnrnd function.
X = [mvnrnd(MU1, SIGMA1, 1000); ...
	 mvnrnd(MU2, SIGMA2, 1000); ...
	 mvnrnd(MU1, SIGMA1, 5000)];
size(X)
scatter(X(:,1), X(:,2), 10, '.');
hold on;

%To fit a Gaussian mixture distribution model to data.
%Use function [fitgmdist]
options = statset('Display', 'final');
obj = fitgmdist(X, 2, 'options', options);
% with 2 clusters, the log-likelihood value is smalles.
h = ezcontour(@(x,y)pdf(obj, [x,y]), [-8 6], [-8 6]);