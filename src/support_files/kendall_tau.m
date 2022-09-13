function tau = kendall_tau(scoreX,scoreY)
%Calculates the Kendall rank correlation coefficient between two score
%vectors. Assumes NO TIES.
%
%INPUTS 
% scoreX,scoreY:    1D vectors of scores according to different metrics
%
%OUTPUTS
% -1<=tau<=1     1: perfect agreement
%               -1: perfect disagreement (reverse order)

%Examples:
%Reverse order=> tau = -1
%sX=sort(rand(10,1));
%sY=sort(sX,'descend');
%tau = kendall_tau(sX,sY);
%
%Random scores=> tau = close to zero
%tau = kendall_tau(rand(100,1),rand(100,1));

%Y.Kamer
%20200614
%Frankfurt

    numS        = numel(scoreX);
    allPairs    = nchoosek(1:numS,2); %generate pairs i<j
    numAllP = size(allPairs,1); %total pairs
    matX    = scoreX(allPairs); %i,j according to scoreX
    matY    = scoreY(allPairs); %i,j according to scoreY
    
    %number of concurrent pairs
    numConc = sum((matX(:,1)>matX(:,2) & matY(:,1)>matY(:,2))...
                | (matX(:,1)<matX(:,2) & matY(:,1)<matY(:,2)));
    
    tau = (2*numConc-numAllP)/numAllP;
end
        