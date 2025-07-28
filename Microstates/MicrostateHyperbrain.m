function [CT, CoOcc] = MicrostateHyperbrain(seq1, seq2)

% Calculate an index of microstate co-occurrence
% Input: Sequences of microstates: seq1 and seq2 
% Output: CT (cooccurrence index for each microstate-angle association
%         CoOcc: mean of CT

if length(sequence1) ~= length(sequence2)
    error('Sequences must have the same length')
end
Nstates = 7;
phi = linspace(0,2*pi,Nstates+1);

Pe = perms(1:7);
ind = find(Pe(:,1)==1);
Pe = Pe(ind,:);

seqT1 = zeros(size(Pe,1),length(seq1));
seqT1(1,:) = seq1;
seqT2 = zeros(size(Pe,1),length(seq2));
seqT2(1,:) = seq2;

cont = 2;
for iter = size(Pe,1)-1:-1:1
    for k=1:7
        seqT1(cont,find(seq1==k)) = Pe(iter,k);
        seqT2(cont,find(seq2==k)) = Pe(iter,k);
    end
    cont = cont+1;
end

CT = zeros(1,size(seqT1,1));
for iter = 1:size(seqT1,1)
    % iter
    seqR1 = seqT1(iter,:);
    seqR2 = seqT2(iter,:);

    for k=1:length(seqR1)
        diff_phi = phi(seqR1(k))-phi(seqR2(k));
        pos(k) = exp(sqrt(-1)*diff_phi );
    end
    CT(iter) = abs(mean(pos));

    clear seqR1 seqR2
end
CoOcc = mean(sqrt(CT)); %mean(CT)

