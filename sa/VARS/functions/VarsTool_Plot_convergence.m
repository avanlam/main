function VarsTool_Plot_convergence(VARS_out,Freq,iter,fig2, name, btsrpFlg)

% =============================================================== %
% Description:
% This function is a part of VARS TOOL software. It plots the convergence
% of the selected sensitivity metric, the factor rankings, and reliability of
% estimated rankings base on selected metric for the specified reporting
% frequency.
%
% =============================================================== %
% Inputs:
%     VARS_out ==========>>> VARS output file
%       Freq   ==========>>> Reporting frequency : report after each x stars are completed
%       iter   ==========>>> Iteration number
%
%  Outputs:
%         figure (1x3 subplots)
% =============================================================== %
% Note:
% The default sensitivity metric is "IVARS50".
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

fig2.WindowState = 'fullscreen';

h = VARS_out.IVARSid{end};
idx = h==0.5; %default value. Any value from (0-1). For example H==0.3 corresponds to 'IVARS30'
counter = iter/Freq;

% assigning the calculated sensitivity metric, factor ranking, and estimated reliability
for i=1:counter
    metric(i,:) = VARS_out.IVARS{1, i*Freq}(idx,:);
    metric_rank(i,:) = VARS_out.rnkIVARS{1,i*Freq}(idx,:);
    if btsrpFlg == 1
    metric_rel(i,:) = 100*VARS_out.relIVARS{1,i*Freq}(idx,:);
    end
end

%creating label for factors
if nargin==6
    labels = name';
    for i = 1:numel(labels)
        labels{i} = strcat('$$', labels{i},'$$');
    end
    labels{end} = '$$\Omega$$';
else
    for kk=1:size(metric,2)
        labels{kk} = strcat('x_','{',num2str(kk),'}');
    end
end

% finding the perfect position for figure
plotheight = 20; plotwidth=16; subplotsx = 3;
subplotsy = 1; leftedge=1.35; rightedge = 0.4;
topedge = 1; bottomedge = 1.5; spacex = 0.65; spacey = 0.2;
sub_pos = subplots_pos(plotwidth,plotheight,leftedge,rightedge,...
    bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

f(1)=subplot(1,3,1);

% This subplot shows convergence of the selected sensitivity metric
for kk = 1:size(metric,2)
    [mrk,lS,FaceClrs,clrs] = GetLineStyleForPlot(kk);
    plot(metric(:,kk),...
        'color',clrs,...
        'LineWidth',1.25,...
        'LineStyle',lS,...
        'Marker',mrk,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',FaceClrs,...
        'MarkerSize',6);
    xlabel('Star Number',...
        'interpreter','latex','FontSize',18) % x-axis label
    ylabel('IVARS50','interpreter','latex','FontSize',18) % y-axis label
    xlim([1 counter+1])
    hold on
end
ax = gca;
ax.XTick = (1:counter);
ax.XTickLabel = (Freq:Freq:iter);
set(ax,'fontsize',18,'fontname','Calibri','LineWidth',1.25,'YScale','log');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperSize', [plotwidth plotheight]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
% ax.Position = sub_pos{1,1};

f(2)=subplot(1,3,2);

% This subplot shows convergence of the factor rankings based on selected sensitivity metric
for kk = 1:size(metric_rank,2)
    [mrk,lS,FaceClrs,clrs] = GetLineStyleForPlot(kk);
    plot(metric_rank(:,kk),...
        'color',clrs,...
        'LineWidth',1.25,...
        'LineStyle',lS,...
        'Marker',mrk,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',FaceClrs,...
        'MarkerSize',6);
    xlabel('Star Number',...
        'interpreter','latex','FontSize',18) % x-axis label
    ylabel('Ranking based on IVARS50',...
        'interpreter','latex','FontSize',18) % y-axis label
    xlim([1 counter+1])
    ylim([1 size(metric_rank,2)])
    hold on
end
ax = gca;
ax.XTick = (1:counter);
ax.YTick = (1:size(metric_rank,2)+1);
ax.XTickLabel = (Freq:Freq:iter);
set(ax,'FontName','CMU Serif','FontSize',18,'LineWidth',1.25);
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperSize', [plotwidth plotheight]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
% ax.Position = sub_pos{2,1};

if btsrpFlg == 1
f(3)=subplot(1,3,3);
% This subplot shows convergence of the estimated robustness of the factor ranking
for kk = 1:size(metric_rel,2)
    [mrk,lS,FaceClrs,clrs] = GetLineStyleForPlot(kk);
    plot(metric_rel(:,kk),...
        'color',clrs,...
        'LineWidth',1.25,...
        'LineStyle',lS,...
        'Marker',mrk,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',FaceClrs,...
        'MarkerSize',6);
    xlabel('Star Number',...
        'interpreter','latex','FontSize',18) % x-axis label
    ylabel('Robustness of Rankings based on IVARS50 ($\%$)',...
        'interpreter','latex','FontSize',18) % y-axis label
    xlim([1 counter+1])
    ylim([1 100])
    hold on
end
leg = legend(labels(1:size(metric,2)),'Location','best','Orientation','vertical','LineWidth',1,...
    'interpreter','latex','FontSize',18,'NumColumns',2);
axl = gca;
axl.XTick = (1:counter);
axl.XTickLabel = (Freq:Freq:iter);
% axl.Position = sub_pos{3,1};
set(axl,'fontsize',18,'fontname','CMU Serif','LineWidth',1.25);
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperSize', [plotwidth plotheight]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
end

% for ki = 1:3
%     if ki==3
%         copyobj([leg,axl],fig2);
%     else
%         copyobj(f(ki),fig2);
%     end
% end


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [ positions ] = subplots_pos(plotwidth,plotheight,leftmargin,rightmargin,...
            bottommargin,topmargin,nbx,nby,spacex,spacey)
        % This function finds the perfect positions for subplots
        
        subxsize=(plotwidth-leftmargin-rightmargin-spacex*(nbx-1.0))/nbx;
        subysize=(plotheight-topmargin-bottommargin-spacey*(nby-1.0))/nby;
        
        for mn=1:nbx
            for j=1:nby
                
                xfirst=leftmargin+(mn-1.0)*(subxsize+spacex);
                yfirst=bottommargin+(j-1.0)*(subysize+spacey);
                
                positions{mn,j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
                
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end