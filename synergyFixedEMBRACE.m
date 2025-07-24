function [W, H, VAF] = synergyFixedEMBRACE(dataIn, nSyn)

[W, H, VAF] = NN_mat_fact_sparse(dataIn, nSyn, 500);