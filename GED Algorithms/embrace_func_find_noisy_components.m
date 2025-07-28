function [removed_evals, no_removed_eigmax, no_removed_mad, no_removed_knee, no_removed_all] = embrace_func_find_noisy_components(evals,eig_threshold)
% This function uses four methods to find the srtifsctusl components:
%   1- Based on the distribution of eigenvalue_max that is identified from clean data
%   2- Based on MAD (the rmoutlier function in MATLAB)
%   3- Based on Knee/Elbow detection
%   4- Based on the Effective number (eq 15, https://doi.org/10.1016/j.jneumeth.2014.01.024)

% based on eigenvalue_max
removed_evals_1 = evals >= eig_threshold;
% removed_evals_1 = true(length(evals),1);

% based on MAD
mad_threshold = 2 * 1.4826 * median(abs(evals-median(evals)));
% removed_evals_2 = evals >= mad_threshold;
removed_evals_2 = true(length(evals),1);

% based on knee/elbow
% removed_evals_3      = false(length(evals),1);
% m                    = embrace_func_knee_point(log10(evals));
% removed_evals_3(1:m) = true;
removed_evals_3 = true(length(evals),1);

% based on the effective number

% final result
removed_evals = removed_evals_1 & removed_evals_2 & removed_evals_3;

no_removed_eigmax = sum(removed_evals_1);
no_removed_mad    = sum(removed_evals_2);
no_removed_knee   = sum(removed_evals_3);
no_removed_all    = sum(removed_evals);

end