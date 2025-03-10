%% Code to Generate Synthetic Data for Experiments %%

clear; clc; close all;
addpath(genpath(pwd));
addpath('../../matlab/tproduct toolbox 2.0 (transform)/')
%warning('off','all')

% Define all parameters below to generate syntehtic data (see hyperparameters.txt)

l = 5; p = 4; n = 10; m = 25;
tdims = [l,p,n,m];

beta_array = [0.025,0.075,0.1];
betarow_array = [0.4,0.6];

q_array = [1, 1 - beta_array]; % OR under-estimate OR over-estimate array

num_trials = 150;
num_its = 2000;

% String or Array
corr_option = "large";

% Examples:
if corr_option == "large"
    cor_size = [100,20];
elseif corr_option == "small"
    cor_size = [10,5];
end


mQTRK_QTRK_plots(tdims, beta_array, betarow_array, q_array, num_trials, num_its, cor_size, corr_option)

