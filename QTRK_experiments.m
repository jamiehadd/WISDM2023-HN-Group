%% Code to Generate Synthetic Data for Experiments %%

clear; clc;
addpath(genpath(pwd));

% Define all parameters below to generate syntehtic data

l = 5; p = 4; n = 10; m = 12;
tdims = [l,p,n,m];

num_corrupt_array = [10,20]; 
k_array = [1,3,6];

q_array = 1-num_corrupt_array/(m*p*n); % OR under-estimate OR over-estimate array

num_trials = 50;
num_its = 2000;

% String or Array
corr_option = "small"; 

% Examples:
if corr_option == "large"
    cor_size = [100,20]; 
elseif corr_option == "small"
    cor_size = [10,5];
else 
    cor_size = corr_option;
end


QTRK_trials_plots(tdims, num_corrupt_array, k_array, q_array, num_trials, num_its, cor_size)
