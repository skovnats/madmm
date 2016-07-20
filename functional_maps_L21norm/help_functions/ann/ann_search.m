function [indices,dists] = ann_search( data, testset, n )

%- data - each column is a point
tree = ann('init', data );
[indices, dists] = ann('search', tree, testset, n, 'eps', 1.1);
ann('deinit', tree);
clear ann;