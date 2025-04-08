% S. Gribling, L. Sinjorgo and R. Sotirov (April 2025)
clear;
s = 5;
timer = tic;
[Y,h,fullEdgeList,cons] = parseQMC_swap(s);

% for s = 5, we only need to check the cycle graph C_5, 
% consisting of the following edges
cycleEdges = [1,2;2,3;3,4;4,5;1,5];

% we need to match these edges with the corresponding h variable
[~, h_idx] = ismember(cycleEdges, fullEdgeList, 'rows');
obj = sum(h(h_idx));

% set YALMIP parameters
ops = sdpsettings;
ops.verbose = 0;
ops.solver = 'mosek';
ops.dualize = 1;

% solve the SDP using YALMIP (default in YALMIP is minimization)
yalmipDiagnostics = optimize(cons,-obj,ops);

optValue = value(obj);

% check c(C_5,2) <= 2
if optValue > floor(s/2)
    error("False lemma!")
else
    fprintf("Case s = " + num2str(s) + "  verified correctly, in " + num2str(toc(timer)) + " seconds.\n");
end