%% Code to Generate Synthetic Data for Experiments %%

clear; clc;
addpath(genpath(pwd));
addpath('../../matlab/tproduct toolbox 2.0 (transform)/')

% Define all parameters below to generate syntehtic data

l = 5; p = 4; n = 10; m = 25;
tdims = [l,p,n,m];

num_corrupt_array = [25,75,100];
k_array = [5,10,20];

q_array = [1, 1 - num_corrupt_array/(m*p*n)]; % OR under-estimate OR over-estimate array

num_trials = 150;
num_its = 2000;

% Specify model
alg = "mQTRK"; % or alg = "mQTRK";

% String or Array
corr_option = "large";

% Examples:
if corr_option == "large"
    cor_size = [100,20];
elseif corr_option == "small"
    cor_size = [10,5];
end


ALG_plots(alg, tdims, num_corrupt_array, k_array, q_array, num_trials, num_its, cor_size, corr_option)

