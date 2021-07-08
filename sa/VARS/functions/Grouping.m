function [ Groups, relGrp, relGrp_all ] = Grouping( IVARS, H, ST, numDim, numGrp, num_samples, rankST, rankIVARS, ...
                                                   bnchmrkrankST, bnchmrkrankIVARS, plotFlg, name )
% Author: Saman Razavi in November 2019

%%cluster_analysis (ST);
sens_MatrixIvars50 = squeeze(IVARS((H==0.5), :,:));
sens_MatrixST = ST';

if nargin == 12
[numGrp_IVARS50,IVARS50_grp,ClustersIVARS50] = factorGrouping ( squeeze(IVARS((H==0.5), :,:)) , numGrp, plotFlg, name);
[numGrp_ST,ST_grp,ClustersST] = factorGrouping ( ST' , numGrp, plotFlg, name);
else
[numGrp_IVARS50,IVARS50_grp,ClustersIVARS50] = factorGrouping ( squeeze(IVARS((H==0.5), :,:)) , numGrp, plotFlg, name);
[numGrp_ST,ST_grp,ClustersST] = factorGrouping ( ST' , numGrp, plotFlg, name);
end

Groups = [ST_grp, IVARS50_grp]';

for g = 1:numGrp_ST
    cluster_ST{1,g} = find(ST_grp==g)';
    cluster_rankST{1,g} = bnchmrkrankST(cluster_ST{1,g});
    cluster_rankST{1,g} = sort(cluster_rankST{1,g});
end
for g = 1:numGrp_IVARS50
    cluster_IVARS50{1,g} = find(IVARS50_grp==g)';
    cluster_rankIVARS50{1,g} = bnchmrkrankIVARS((H==0.5),cluster_IVARS50{1,g});
    cluster_rankIVARS50{1,g} = sort(cluster_rankIVARS50{1,g});
end

reliST = zeros(num_samples,numDim); reliIVARS50 = zeros (num_samples,numDim);
for D = 1 : numDim
    Match = cellfun(@(x) find(x==D), cluster_ST, 'UniformOutput', 0);
    rank_rangeST = ~cellfun(@isempty,Match);
    rankST_benchmark_grp = cluster_rankST{1,rank_rangeST};
    
    Match = cellfun(@(x) find(x==D), cluster_IVARS50, 'UniformOutput', 0);
    rank_rangeIVARS50 = ~cellfun(@isempty,Match);
    rankIVARS50_benchmark_grp = cluster_rankIVARS50{1,rank_rangeIVARS50};
    
    for iter = 1 : num_samples
        reliST(iter, D) = length ( find( rankST(iter, D) == rankST_benchmark_grp ) ) / num_samples;
        reliIVARS50(iter, D) = length ( find( rankIVARS((H==0.5),D, iter) == rankIVARS50_benchmark_grp ) ) / num_samples;
    end
end
relGrp(1, :) = sum(reliST);
relGrp(2, :) = sum(reliIVARS50);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for kk = 1:numDim-1
    numgrpST = max(ClustersST(:,kk));
    for g = 1:numgrpST
        clusterST{1,g} = find(ClustersST(:,kk)==g)';
        clusterrankST{1,g} = bnchmrkrankST(clusterST{1,g});
        clusterrankST{1,g} = sort(clusterrankST{1,g});
    end
    numgrpivars = max(ClustersIVARS50(:,kk));
    for g = 1:numgrpivars
        clusterIVARS50{1,g} = find(ClustersIVARS50(:,kk)==g)';
        clusterrankIVARS50{1,g} = bnchmrkrankIVARS((H==0.5),clusterIVARS50{1,g});
        clusterrankIVARS50{1,g} = sort(clusterrankIVARS50{1,g});
    end
    
    reliST1 = zeros(num_samples,numDim); reliIVARS501 = zeros (num_samples,numDim);
    for D = 1 : numDim
        Match = cellfun(@(x) find(x==D), clusterST, 'UniformOutput', 0);
        rankrangeST = ~cellfun(@isempty,Match);
        rankSTbenchmarkgrp = clusterrankST{1,rankrangeST};
        
        Match = cellfun(@(x) find(x==D), clusterIVARS50, 'UniformOutput', 0);
        rankrangeIVARS50 = ~cellfun(@isempty,Match);
        rankIVARS50benchmark_grp = clusterrankIVARS50{1,rankrangeIVARS50};
        
        for iter = 1 : num_samples
            reliST1(iter, D) = length ( find( rankST(iter, D) == rankSTbenchmarkgrp ) ) / num_samples;
            reliIVARS501(iter, D) = length ( find( rankIVARS((H==0.5),D, iter) == rankIVARS50benchmark_grp ) ) / num_samples;
        end
    end
    relGrp_all{1}(kk,:) = sum(reliST1);
    relGrp_all{2}(kk,:) = sum(reliIVARS501);
end

end

%%
function [optm_numGrp, rank_grp, Clusters,leafOrder,S] = factorGrouping ( sensIndx , numGrp, plotFlg, name )
[n,m] = size(sensIndx);
% Transformation
%Replacing zero elements. This is due to numerical reason
M=sort(sensIndx);
uu=sum(~any(sensIndx,2))+1;
replace = M(uu,:);
sensIndx( ~any(sensIndx,2), : ) = repmat(replace,uu-1,1);
%
R=[];
for i=1:m
    R=[R;sensIndx(:,i)];
end

[TRANSDAT, LAMBDA] = BoxCox(R); %BOX-COX transformation
if LAMBDA <= 0.0099
    TRANSDAT = log(R);
end
[idx,jdx]=find(isinf(TRANSDAT));
TRANSDAT(idx,jdx) = log(min(R(R>0)));
S = reshape(TRANSDAT,[n,m]);

% Agglomerative hierarchical cluster
Z = linkage(S, 'ward', 'euclidean');
%Z(:,3) = (Z(:,3)-min(Z(:,3)))./(max(Z(:,3))-min(Z(:,3)));

% Optimal group number
Clusters = cluster(Z,'maxclust',2:n);
if numGrp
    rank_grp = cluster(Z,'maxclust',numGrp);
    optm_numGrp = numGrp;
    nn=1; id = size(Z,1);
    while  nn ~= optm_numGrp
        cutoff = Z(id,3);
        rank_grp = cluster(Z,'cutoff',cutoff,'criterion','distance');
        nn = max(rank_grp);
        id = id-1;
    end
    clrThrshl = 0.5*(Z(id+1,3) + Z(id+2,3));
else
    [cutoff , clrThrshl] = elbow_method(Z);
    % [cutoff , clrThrshl] = maxGap_method(Z);
    rank_grp = cluster(Z,'cutoff',cutoff,'criterion','distance');
    optm_numGrp = max(rank_grp);
end

% Generating dendogram
if plotFlg == 1
    fig3 = figure(3);
    fig3.WindowState = 'fullscreen';
    D = pdist(S);
    leafOrder = optimalleaforder(Z,D);

    if nargin==4
        labels = name';
        for i = 1:numel(labels)
            labels{i} = strcat('$$', labels{i},'$$');
        end
         labels{end} = '$$\Omega$$';
    else
        for kk=1:size(X,2)
            labels{kk} = strcat('\it x_','{',num2str(kk),'}');
        end
    end
    
    H = dendrogram(Z,0,'Reorder',leafOrder,'Orientation','top','Labels',labels,'ColorThreshold',clrThrshl);
    set(H,'LineWidth',2.5)
    xlabel('Factors', 'interpreter','latex','FontSize',24) % x-axis label
    ylabel('Dissimilarity Metric', 'interpreter','latex','FontSize',24) % x-axis label

    %Changing the colours
    lineColours = cell2mat(get(H,'Color'));
    colourList = unique(lineColours, 'rows');
    myColours = [ 0.9961    0.4961         0
        1         0.2       0.6
        0.1293    0.8060    0.8679
        0.7674    0.2331    0.0641
        0.4067    0.3567    0.6615
        0.0724    0.9594    0.0988
        0.5533    0.5938    0.4948
        0.9439    0.1561    0.7330
        0.8500    0.3250    0.0980
        1 0 0
        1 1 0
        1 0 1
        0 1 0
        0.1922    0.7216    0.7647
        0.3008    0.6836    0.2891
        0.4940    0.1840    0.5560
        0.6350    0.0780    0.1840
        0.8906    0.1016    0.1094
        0.85       0.21       0.52
        0.6       0.4       0.67
        0.48       0.61       0.7
        0.23       0.85       0.87
        0         1         0.9
        0.5137    0.2510    0.6588
        0.5       0.5       0.5
        1         0.2       0.6
        0         0.4       0.4
        0.4       1         1
        0.1958    0.8952    0.8133
        0.6559    0.9467    0.5456
        0.8300    0.7180    0.8609
        0.6015    0.7941    0.1034
        0.5082    0.3411    0.1797
        0.0803    0.7204    0.1911
        0.8334    0.4518    0.9423
        0.9037    0.2637    0.8813
        0.0995    0.6714    0.3421
        0.7394    0.6726    0.2611
        0.4270    0.6758    0.5042
        0.7625    0.5729    0.9960
        0.7468    0.5388    0.5591
        0.4596    0.0827    0.6215];
    
    %Replace each colour (colour by colour). Start from 2 because the first colour are the "unclustered" black lines
    for colour = 2:size(colourList,1)
        %Find which lines match this colour
        idx = ismember(lineColours, colourList(colour,:), 'rows');
        %Replace the colour for those lines
        lineColours(idx, :) = repmat(myColours(colour-1,:),sum(idx),1);
        legidx(colour-1) = find(idx, 1, 'first');
    end
    %Apply the new colours to the chart's line objects (line by line)
    for line = 1:size(H,1)
        set(H(line), 'Color', lineColours(line,:));
    end

    colourList = unique(lineColours, 'rows');
    leg = {'1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th'};
    lgd = legend(H(legidx),leg(1:size(colourList,1)), 'interpreter', 'latex', 'Location', 'best');
    title(lgd,'Influence of the clusters')

end

end

%%
% Finding optimal group number

function [cutoff , clrThrshl] = elbow_method(Z)

Q1 = [1,Z(1,3)]; Q2 = [size(Z,1), Z(end,3)];
for i = 1:(size(Z,1)-2)
    P = [i+1,Z(i+1,3)];
    d(i) = abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1); % for row vectors.
end

[~,id] = max(d);
cutoff = Z(id+1,3);
clrThrshl = 0.5*(Z(id+1,3) + Z(id+2,3));
end

% function [cutoff , clrThrshl] = maxGap_method(Z)
% 
% for i=1:(size(Z,1)-1)
%     gap(i)=abs(Z(i,3)-Z(i+1,3));
% end
% 
% [~,lsor] = findpeaks(gap,'SortStr','descend');
% %findpeaks(y)
% if isempty(lsor)
%     M = max(Z);
%     cutoff = M(3)-eps;
% else
%     cutoff = Z(lsor(1,1)+1,3);
%     clrThrshl = 0.5*(Z(lsor(1,1)+1,3) + Z(lsor(1,1)+2,3));
% end
% end
% EOF

