function [shape] = loadfaust_i(j,k)
% Loads shape i with k bases
%% params
if ismac
path='data/faust/diam=200';
pathdesc='data/faust/diam=200/descs/shot_9_7.16_3'; %% 'descs\shot_9_7.16_3'
patheigs='data/faust/diam=200/eigendec';
pathLaplacian='data/faust/diam=200/lbos';
pathMeshes='data/faust/diam=200/meshes';    
else
path='D:\Research\deynard\data\faust\diam=200';
pathdesc='D:\Research\deynard\data\faust\diam=200\descs\shot_9_7.16_3'; %% 'descs\shot_9_7.16_3'
patheigs='D:\Research\deynard\data\faust\diam=200\eigendec';
pathLaplacian='D:\Research\deynard\data\faust\diam=200\lbos';
pathMeshes='D:\Research\deynard\data\faust\diam=200\meshes';
end


%%
name=sprintf('tr_reg_%.3d.mat',j);
%% mesh
load(fullfile(pathMeshes,name));
%% Laplacian
load(fullfile(pathLaplacian,name));
% shape.W=W; 
shape.A=A;
%% Eigs
load(fullfile(patheigs,name));
shape.Phi=Phi(:,1:k);
shape.DD=Lambda(1:k);
%% Desc/Shot
load(fullfile(pathdesc,name))
shape.F=desc;
shape.fcoef=shape.F.' * shape.A * shape.Phi;

%% Load symmetry indeces
% load('D:\Research\deynard\data\faust\FAUST_registrations_idxs_sym.mat');
% shape.idxs=idxs;

%% sampling
shape.idx=fps(1,150-1,shape);
end