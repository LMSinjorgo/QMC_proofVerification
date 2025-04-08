% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
clear;
timer = tic;
fprintf("Verifying case s = 9. Expected runtime: 20 seconds. \n")
s = 9;
[Y,h,edgeList,cons] = parseQMC_swap(s);

% load the graphs (in matlab) for which we need to verify c(G,k) <= 4
numGraphs = 219;
fileName = "humanReadable_s" + num2str(s) + ".txt";
[edgeIndicator] = nautyToMATLAB(fileName,s,numGraphs);

% compute the cover numbers as upper bound on c(G,k)
[coverNumsUB] = computeCoverNumbers(edgeIndicator);

% remove the graphs that have tau(G) <= floor(s/2);
edgeIndicator(:,coverNumsUB <= floor(s/2)) = [];
numGraphs = size(edgeIndicator,2);

% set YALMIP parameters
ops = sdpsettings;
ops.verbose = 0;
ops.solver = 'mosek';
ops.dualize = 0;

% construct objective functions
symbolicObjVec = sdpvar(numGraphs,1);
objVal = zeros(1,numGraphs);
for k = 1:numGraphs
    symbolicObjVec(k) = sum(h(edgeIndicator(:,k)),1);
end

% solve using yalmip
yalmipDiagnostics = optimize(cons,-symbolicObjVec,ops);

totalSolverTime = sum(yalmipDiagnostics.solvertime);
yalmipTime = yalmipDiagnostics.yalmiptime;
totalTime = totalSolverTime+yalmipTime;

for k = 1:numGraphs
    objVal(k) = value(symbolicObjVec(k),k);
end

% check if the lemma is true
if any(objVal > floor(s/2))
    error("False lemma!")
else
    fprintf("Case s = " + num2str(s) + "  verified correctly, in " + num2str(toc(timer)) + " seconds.\n");
end
