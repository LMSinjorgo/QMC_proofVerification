% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
clear;
s = 11;
timer = tic;
fprintf("Verifying case s = 11. Expected runtime: 1 minute. \n")

% load the graphs for which we need to verify c(G,k) <= floor(s/2);
numGraphs = 26360;
fileName = "humanReadable_s" + num2str(s) + ".txt";
[edgeIndicator] = nautyToMATLAB(fileName,s,numGraphs);

% compute the cover numbers as upper bound on c(G,k)
% we have precomputed the cover numbers, so we will only recompute the
% cover numbers for the graphs for which we know tau(G) <= 5
load("hasLowCoverNum_s11.mat");
[coverNumsUB] = computeCoverNumbers(edgeIndicator(:,hasLowCoverNum_s11));

% remove the graphs for which tau(G) <= floor(s/2)
edgeIndicator(:,hasLowCoverNum_s11(coverNumsUB <= floor(s/2))) = [];

% now we compute c(G,s) via computing lambda_max(H_G)
[cGk_val] = computeEigenValueBounds(edgeIndicator,s);

% check if the lemma is true
if any(cGk_val > floor(s/2))
    error("False lemma!")
else
    fprintf("Case s = " + num2str(s) + " verified correctly, in " + num2str(toc(timer)) + " seconds.\n");
end
