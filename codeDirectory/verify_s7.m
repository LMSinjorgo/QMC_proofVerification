% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
clear;
s = 7;
timer = tic;
[Y,h,edgeList,cons] = parseQMC_swap(s);

% load the 6 graphs (in matlab) for which we need to verify c(G,k) <= 3
numGraphs = 6;
fileName = "humanReadable_s" + num2str(s) + ".txt";
[edgeIndicator] = nautyToMATLAB(fileName,s,numGraphs);

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

if any(objVal > floor(s/2))
    error("False lemma!")
else
    fprintf("Case s = " + num2str(s) + "  verified correctly, in " + num2str(toc(timer)) + " seconds.\n");
end
