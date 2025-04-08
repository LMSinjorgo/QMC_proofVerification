function [edgeIndicator] = nautyToMATLAB(fileName,s,numGraphs)
% convert the human readable output from nauty into a matlab matrix.
% Each column of the returned matrix corresponds to the edge-incidence
% matrix of a single graph.

fid = fopen(fileName);
graphIdx = 1;

lowerTriuIdx = find(tril(ones(s),-1));
mappingMat = zeros(s);

maxNumEdges = numel(lowerTriuIdx);

mappingMat(lowerTriuIdx) = 1:maxNumEdges;
mappingMat = mappingMat+mappingMat';

edgeIndicator = false(maxNumEdges,numGraphs);
lineCounter = 1;
sparseIdx = zeros(round(0.5*maxNumEdges*numGraphs),2,"uint32");
sparseIdxCounter = 1;
while true
    % we are only interested in the line number when the line number is
    % multiple of 4. These lines contain the edges, and the other lines are
    % not interesting
    fgets(fid);
    fgets(fid);
    fgets(fid);
    newLine = fgets(fid);
    
    % read the edges in the form {v1,v2}
    vertexValues = sscanf(newLine,'%u');
    numEdges = numel(vertexValues)/2;
    edges = reshape(vertexValues,2,numEdges)';
    linIdxEdges = edges*[1;s]-s;
    linIdxEdges = mappingMat(linIdxEdges);
    sparseIdx(sparseIdxCounter:sparseIdxCounter+numEdges-1,:) = [linIdxEdges,repmat(graphIdx,numEdges,1)];
    sparseIdxCounter= sparseIdxCounter+size(linIdxEdges,1);

    if sparseIdxCounter >= size(sparseIdx,1)-2*maxNumEdges
        sparseIdx = [sparseIdx; zeros(round(0.5*maxNumEdges*(numGraphs-graphIdx)),2,"uint32")];
    end

    graphIdx = graphIdx + 1;
    if graphIdx > numGraphs
        break;
    end
end
fclose('all');

firstZero = find(sparseIdx(:,1) == 0,1,'first');
if ~isempty(firstZero)
    sparseIdx = sparseIdx(1:firstZero-1,:);
end
edgeIndicator = sparse(sparseIdx(:,1),sparseIdx(:,2),true,maxNumEdges,numGraphs);


end

