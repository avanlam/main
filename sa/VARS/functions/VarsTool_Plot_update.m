function VarsTool_Plot_update(VARS_out,fig1,delta, name, btsrpFlg)
% (do bootstarp) param name
% =============================================================== %
% Description:
% This function is a part of VARS TOOL software. It plots the Variograms,
% piecharts for (a)Morris (Morris-ABE),(b)Sobol (VARS-TO), and (c)IVARS50 metrics,
% and the respective confidence intervals (CIs) estimated based on the bootstrapping
%
% =============================================================== %
% Inputs:
%       VARS_out ==========>>> VARS output file
%         fig   ==========>>> primary figure with default property value
%
%  Outputs:
%         figure
% =============================================================== %
% Note:
% The size of each slice in piecharts represents the proportion of the sensitivity metric
% value of the respective parameter to the sum of all sensitivity metric values
% over all parameters.
% =============================================================== %
% References:
% (1) S. Razavi, & H.V. Gupta (2016):
%     A new framework for comprehensive, robust, and efficient global sensitivity analysis:
%     1. Theory. Water Resources Research, 52,423-439.
% (2) S. Razavi & H.V. Gupta (2016):
%     A new framework for comprehensive, robust, and efficient global sensitivity analysis:
%     2. Application, Water Resources Research, 52, 440–455.
%----------------------------------------------------------------
% Programmed by Razi Sheikholeslami & Saman Razavi (2017)
% Global Institute for Water Security, School of Environment
% and Sustainability, University of Saskatchewan
% E-mail: razi.sheikholeslami@usask.ca;
%         saman.razavi@usask.ca
% ----------------------------------------------------------------

fig1.WindowState = 'fullscreen';

% assigning the sensitvity metrics
gamma = VARS_out.Gamma{1, end}(1:end,:);
X = VARS_out.IVARS{1, end};
H = VARS_out.IVARSid{end};
X(X == 0) = eps;
TO = VARS_out.ST{1,end};
TO(TO == 0) = eps;
morris = VARS_out.MAEE{1,end}(1,:);
morris(morris == 0) = eps;
MAEE=morris/sum(morris);
Vars_TO=TO/sum(TO);

p3 = H==0.5;
%p1 = h==0.1; p2 = h==0.3;

for JJ=1:size(X,2)
    %IVARS30(JJ) = X(p2,JJ)/sum(X(p2,:));
    IVARS50(JJ) = X(p3,JJ)/sum(X(p3,:));
end

% assigning factor rankings
[~, order] = sort(morris);
temp = [ order; 1: length(morris) ]';
temp2 = sortrows(temp, 1)';
rankMEE = temp2(2, :);
rank50 = VARS_out.rnkIVARS{1,end}(p3,:);
rankTO = VARS_out.rnkST{1,end};

if nargin==5
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
%%
% plotting directional Variograms
hh = linspace(delta,1-delta,size(gamma,1));
%q(1) = subplot(4,4,[1:3 5:7]);
q(1) = subplot(2,1,1);
for kk = 1:size(gamma,2)
    [mrk,lS,FaceClrs,clrs] = GetLineStyleForPlot(kk);
    plot([0,hh],[0; gamma(:,kk)],...
        'color',clrs,...
        'LineWidth',1.25,...
        'LineStyle',lS,...
        'Marker',mrk,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',FaceClrs,...
        'MarkerSize',6);
    hold on
    xlim([0 0.5])
end
hold off
ax=gca;
ax.XTick = hh;
xlabel('h (perturbation scale)',...
    'interpreter','latex','FontSize',18) % x-axis label
ylabel('Variogram $$\it \gamma(h)$$','interpreter','latex','FontSize',18) % y-axis label
set(ax,'fontsize',18,'FontName','CMU Serif','LineWidth',1.25); %'YScale','log'
legend(labels(1:size(X,2)),'Location','best','Orientation','vertical','LineWidth',1,...
    'interpreter','latex','FontSize',18,'NumColumns',2);
%%
% plotting bar charts
q(2) = subplot(2,1,2);
GSA_matrix = [Vars_TO' IVARS50' MAEE'];
bar(GSA_matrix)
ax11 = gca;
bar1 = bar(GSA_matrix,'Parent',ax11);
set(bar1(3),'FaceColor',[227,26,28]./255,...
    'EdgeColor',[0.0784313753247261 0.168627455830574 0.549019634723663]);
set(bar1(2),'FaceColor',[51,160,44]./255);
set(bar1(1),'FaceColor',[31,120,180]./255,...
    'EdgeColor',[0.0705882385373116 0.211764708161354 0.141176477074623]);
xlabel('Factors',...
    'interpreter','latex','FontSize',18) % x-axis label
ylabel('Ratio of Factor Sensitivity','interpreter','latex','FontSize',18) % y-axis label
xticks(1:numel(labels)); set(gca,'xticklabel',labels)
set(ax11,'fontsize',18,'LineWidth',1.25,'ygrid','off','GridLineStyle','--','GridColor',[0 0 1])

%%
if btsrpFlg==1
    % ploting the estimated confidence intervals based on the bootstrapping for
    % IVARS50 & VARS_TO
    
    IVARSb_lower_bound = VARS_out.IVARSlb{1,end}(p3,:);
    IVARSb_lower_bound(IVARSb_lower_bound == 0) = eps;
    IVARSb_lower_bound = IVARSb_lower_bound/sum(X(p3,:));
    
    IVARSb_upper_bound = VARS_out.IVARSub{1,end}(p3,:);
    IVARSb_upper_bound(IVARSb_upper_bound == 0) = eps;
    IVARSb_upper_bound = IVARSb_upper_bound/sum(X(p3,:));
    
    hold on
    Data = VARS_out.IVARS{1, end}(p3,:);
    Data(Data == 0) = eps;
    x_ax = 1:size(Data,2);
    for nn = 1:size(Data,2)
        plot([x_ax(nn) x_ax(nn)],[IVARSb_lower_bound(nn) IVARSb_upper_bound(nn)],...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerSize',6,...
            'Marker','+',...
            'LineWidth',1,...
            'LineStyle','--',...
            'Color',[0 0 0]);
        hold on
    end
    leg1 = legend('VARS-TO (Sobol)', 'IVARS-50', 'VARS-ABE (Morris)', 'Normalized confidence interval', 'location', 'best');
    set(leg1,...
        'LineWidth',1,...
        'FontWeight','normal',...
        'FontSize',18,...
        'interpreter','latex');
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [mrk,lS,FaceClrs,clrs] = GetLineStyleForPlot(i)
        % This function creates the style for your plot
        % Colors iterate over all colors except for white one
        markers = {'d','o','+','x','h','p','s','^','v','>','<','*','.','none'};
        lineStyles = {'-.', '--', '-', ':'};
        colors = [0 0.4470    0.7410
            0.8500    0.3250    0.0980
            0.4023    0.6602    0.8086
            0.1922    0.7216    0.7647
            1.0000    0.4980    0.0549
            0.2148    0.4922    0.7188
            0.3008    0.6836    0.2891
            0.5938    0.3047    0.6367
            0.9961    0.4961         0
            0.9290    0.6940    0.1250
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.9336    0.5391    0.3828
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840
            0.8906    0.1016    0.1094
            1         0         0.4
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
            1 0 0
            1 1 0
            1 0 1
            0 1 0
            0 0 0
            0.1958    0.8952    0.8133
            0.6559    0.9467    0.5456
            0.8300    0.7180    0.8609
            0.6015    0.7941    0.1034
            0.5082    0.3411    0.1797
            0.0803    0.7204    0.1911
            0.8334    0.4518    0.9423
            0.4769    0.9214    0.2941
            0.1293    0.8060    0.8679
            0.7674    0.2331    0.0641
            0.4067    0.3567    0.6615
            0.4766    0.8316    0.6768
            0.4332    0.7702    0.4648
            0.0724    0.9594    0.0988
            0.5533    0.5938    0.4948
            0.9439    0.1561    0.7330
            0.4580    0.2751    0.8735
            0.9037    0.2637    0.8813
            0.5348    0.8339    0.9787
            0.6090    0.8005    0.7543
            0.4245    0.1317    0.6681
            0.2054    0.1198    0.8910
            0.6722    0.4197    0.2705
            0.0995    0.6714    0.3421
            0.7394    0.6726    0.2611
            0.4270    0.6758    0.5042
            0.7625    0.5729    0.9960
            0.4854    0.9417    0.0432
            0.7468    0.5388    0.5591
            0.4596    0.0827    0.6215];
        mrk = markers{mod(i,numel(markers))+1};
        lS = lineStyles{mod(i,numel(lineStyles))+1};
        clrs = colors(mod(i,numel(colors))+1,:);
        Fcolors = {'y', 'm', 'b', 'r', 'c', 'g', 'k'};
        FaceClrs = Fcolors{mod(i,numel(Fcolors))+1};
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end