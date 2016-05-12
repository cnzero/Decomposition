% clear
clear,close all,clc

% cluster-1
MU1 = [1 2]; % input matrix of means mu
SIGMA1 = [3   0.2; ...
		  0.2 2]; % input array of covariances
% cluster-2
MU2 = [-1 -2];
SIGMA2 = [2 0; ...
	      0 1];

% generate data from a mixture of two bivariate Gaussian distributions
% using the mvnrnd function.
X = [mvnrnd(MU1, SIGMA1, 10); ...
	 mvnrnd(MU2, SIGMA2, 10)];
scatter(X(:,1), X(:,2), 10, '.');
hold on;

% %To fit a Gaussian mixture distribution model to data.
% %Use function [fitgmdist]
options = statset('Display', 'final');
obj = fitgmdist(X, 2, 'options', options)
% with 2 clusters, the log-likelihood value is smalles.
h = ezcontour(@(x,y)pdf(obj, [x,y]), [-8 6], [-8 6]);


% % classification and visualization
idx = cluster(obj, X);
cluster1 = (idx==1);
cluster2 = (idx==2);

% figure
% scatter(X(cluster1, 1), X(cluster1, 2), 10, 'r+');
% hold on;
% scatter(X(cluster2, 1), X(cluster2, 2), 10, 'bo');
% hold off;
% legend('Cluster 1', 'Cluster 2', 'Location', 'NW')

% % posterior probability and visualization
% figure
% P = posterior(obj, X);
% size(P)
% scatter(X(cluster1, 1), X(cluster1, 2), 10, P(cluster1, 1), '+');
% hold on;
% scatter(X(cluster2, 1), X(cluster2, 2), 10, P(cluster2, 1), '+')
% hold off;
% legend('Cluster 1', 'Cluster 2', 'Location', 'NW');
% clrmap = jet(80); colormap(clrmap(9:72, :));
% ylabel(colorbar, 'component 1 posterior probability')

% % soft clustering and visualization
% figure
% P = posterior(obj, X);
% [~, order] = sort(P(:,1));
% plot(1:size(X,1), P(order, 1), 'r-', 1:size(X,1), P(order, 2), 'b-');

% soft clustering and visualization and diagonal coveriance.
% figure
% gm2 = fitgmdist(X, 2, 'CovType', 'diagonal', 'SharedCov', true);
% % with 2 clusters, the log-likelihood value is smalles.
% h = ezcontour(@(x,y)pdf(gm2, [x,y]), [-8 6], [-8 6]);

% figure
% P2 = posterior(gm2, X);
% [~, order] = sort(P2(:,1));
% plot(1:size(X,1), P2(order, 1), 'r-', 1:size(X,1), P2(order, 2), 'b-');


% % Assign New Data to Clusters
% Y = [mvnrnd(MU1, SIGMA1, 50); ...
% 	 mvnrnd(MU2, SIGMA2, 25)];
% idx = cluster(gm2, Y);




% ---------------Discussion
% ----------case-1
% X = [X1, X1(:,1)+X2(:,2)];
% ----problem:
% The columns of X are linearly dependent.
% This can cause ill-conditioned covariance estimates.
% ----solution:
% Fit a Gaussian mixture model again, but use the regularization.
% GMModel = fitgmdist(X, 2, 'Regularize', 0.1)

% -----------case-2
% Select the number of GMM components Using PCA

% -----------case-3
% Determine the Best GMM using AIC--AIC fit statistic
% Akaike information criterion, which provides a means of model selection.
% A trade-off between the goodness of fit of the model and the complexity.

% -----------case-4
% Set Initial Values When Fitting Gaussian Mixture Models

% -----------case-5
% Cluster with Gaussian Mixtures
% functions: cluster and posterior are very useful.
% completed.

% -----------case-6
% Soft Clustering Using Gaussian Mixtures Distributions

% -----------case-7
% Assign New Data to Clusters