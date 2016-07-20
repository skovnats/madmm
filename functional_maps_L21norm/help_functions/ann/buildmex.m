mex -output mexann -I./ ANN.cpp bd_fix_rad_search.cpp bd_tree.cpp kd_fix_rad_search.cpp kd_split.cpp perf.cpp bd_pr_search.cpp brute.cpp kd_pr_search.cpp kd_tree.cpp ann_mex_new.cpp bd_search.cpp kd_dump.cpp kd_search.cpp kd_util.cpp

% Find nearest neighbors of m query points in the collection of n points
d = 2;   % dimensions
n = 1e2; % number of points in the dataset
m = 1e1; % number of query points
k = 3;   % number of nearest neighbors
eps = 1; % tolerance (eps=0 - exact search)

X = randn(d,n);
Y = randn(d,m);

tree = ann('init', X);
[idx, dist] = ann('search', tree, Y, k, 'eps', eps);
% idx = kxm matrix with indices to the n points in the dataset denoting the
%       nearest neigbors of each query point
%       idx(k,m) - index of the k-th nearest neighbor of m-th point
% dist = kxm matrix of corresponding Euclidean distances
ann('deinit', tree);

ann('close');


dist0 = sqrt( bsxfun(@minus, Y(1,:),X(1,:)').^2 + bsxfun(@minus, Y(2,:),X(2,:)').^2 );
[dist0,idx0] = sort(dist0,1);
dist0 = dist0(1:k,:);
idx0 = idx0(1:k,:);
