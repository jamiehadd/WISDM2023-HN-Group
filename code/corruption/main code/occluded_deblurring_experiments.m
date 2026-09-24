%% Code for QTRK and mQTRK Circular Deblurring Experiments using MRI Data %%

clear; clc; close all;
addpath(genpath(pwd));
addpath("../circular deblurring/")
addpath('../tproduct toolbox 2.0 (transform)/')
%warning('off','all')

% Import data
load mri;                               % loads mri video as order-4 tensor
X = squeeze(double(D));                 % removes dimension 1 mode
X = mat2gray(X(:,:,1:12));              % only first 12 frames   

% Hyperparametes
num_its = 10000; % number of iterations
num_corrupt = 10; % number of corruption
%<<<<<<< Updated upstream
q = 1; % quantile value
%=======
%q = 0.9999; % quantile value
%>>>>>>> Stashed changes
k = 1; % number of corrupted rows

% Corruptions magnitude distribution
mean_corrupt = 1;
deviation_corrupt = 0;

% Gaussian Filter
h = fspecial('gaussian',[5,5],2);

% Run Experiments
occluded_deblurring_plots(X, h, num_its, num_corrupt, q, k, mean_corrupt, deviation_corrupt)

