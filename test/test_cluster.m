clear,close all,clc
rng('default');
X = [gallery('uniformdata', [10 3], 12); ...
	 gallery('uniformdata', [10 3], 13)+1.2; ...
	 gallery('uniformdata', [10 3], 14)+2.5];
% 30x3

T = clusterdata(X, 'maxclust', 4);
% cluster labels 30x1

% For deep computation, 
% you can use Matlab logic comprehension with T.

% Can you give out a objective function 
% --under which we can decide the best number of clusters.

% Fisher algorithm may be a good one. 


scatter3(X(:,1), ...
		 X(:,2), ...
		 X(:,3), ...
		 100, T, 'filled')