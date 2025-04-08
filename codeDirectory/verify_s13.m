% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
clear;
s = 13;
fprintf("Verifying case s = 13. Expected runtime: ~16 hours. \n")
timer = tic;    

% load the graphs for which we need to verify c(G,k) <= floor(s/2)
numGraphs = 9035088;
fileName = "humanReadable_s" + num2str(s) + ".txt";
[edgeIndicator] = nautyToMATLAB(fileName,s,numGraphs);

% verify which graphs satisfy |N(S)| <= |S| (with S a stable set, that we restrict |S| = 2). We have precomputed which
% graphs satisfy this condition, so we only recheck for those graphs
load("satisfiesStableSetCond_s13.mat")
[satisfiesStableSetCondition] = stableSetCondition(edgeIndicator(:,satisfiesStableSetCond_s13));

% remove the graphs that satisfy |N(S)| <= |S|
edgeIndicator(:,satisfiesStableSetCond_s13(satisfiesStableSetCondition)) = [];

% compute the cover numbers
% we have precomputed for which graphs tau(G) <= floor(s/2), so that we
% only compute tau(G) for those graphs
load("hasLowCoverNum_s13.mat")
[coverNumsUB] = computeCoverNumbers(edgeIndicator(:,hasLowCoverNum_s13));

% remove graphs that satisfy tau(G) <= floor(s/2)
edgeIndicator(:,hasLowCoverNum_s13(coverNumsUB <= floor(s/2))) = [];

% now we compute c(G,s) for the remaining graphs via computing lambda_max(H_G)
[cGk_val] = computeEigenValueBounds(edgeIndicator,s);

% check if the lemma is true
if any(cGk_val > floor(s/2))
    error("False lemma!")
else
    fprintf("Case s = " + num2str(s) + " verified correctly, in " + num2str(toc(timer)/3600) + " hours.\n");
end
