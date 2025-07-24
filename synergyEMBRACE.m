function [W, H, VAF] = synergyEMBRACE(envIn)

nMuscles = size(envIn,2);
VAF = zeros(1,nMuscles);

for i=1:nMuscles
    [W_tmp{i},H_tmp{i},VAF(i)] = NN_mat_fact_sparse(envIn,i,500);
end

nSyn = min(find(VAF>=0.9));

W = W{nSyn};
H = H{nSyn};