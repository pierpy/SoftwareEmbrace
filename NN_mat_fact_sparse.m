function [W,H,VAF] = NN_mat_fact_sparse(M,k,maxiter)

%[W,H,VAF] = NN_mat_fact(M,k,maxiter)
%
%applies nonnegative matrix factorization to the matrix M, decomposing it
%into two matrixes W and H based on k synergies. The algorithm uses the
%multiplicative update rules defined by Lee&Seung (1999,2001). The
%algorithm ends if the error is not changing more than 0.001 along 30
%iterations.
%
%INPUT:
%
%M = input matrix to be decomposed, containing the muscle envelopes
%k = number of synergies to be extracted from the data M
%maxiter = maximum number of iterations to be applied before the algorith
%automatically stops
%
%OUTPUT:
%
%W = matrix of the synergy vectors
%H = matrix of the synergy activation coefficients
%VAF = VAF

if size(M,1)>size(M,2)
    
    M = M';
    
end

%initialization of W and H matrixes

W = 0.05*rand(size(M,1),k)+0.05;

m = randperm(size(M,1));

for i = 1:k

    W(m(i),i) = W(m(i),i) + 0.7 + 0.1*rand(1);

end



% idx=randi([1 k],1,size(M,1));
% 
% W=0.05*rand(size(M,1),k);
% 
% for i=1:length(idx)
%     W(i,idx(i))=W(i,idx(i))+0.7+0.1*rand(1);
% end

Winit = W;

H = rand(k,size(M,2));

W_story(:,:,1) = Winit;

count = 1;

while count < maxiter
    
    
    num_H = W'*M;
    den_H = W'*W*H;
    H = H.*num_H./den_H;
    
    
    num_W = M*H';
    den_W = W*H*H';
    W = W.*num_W./den_W;
    
    W_story(:,:,count+1) = W;
    
    err_tmp(count) = sqrt(sum(sum((M-W*H).^2)));
    
    if count>30 && abs(err_tmp(count) - err_tmp(count-30)) < 0.00001
        
        break
        
    end
    

    
    count = count+1;
    
end

VAF = 1 - sum(sum((M-W*H).^2))./sum(sum((M).^2));