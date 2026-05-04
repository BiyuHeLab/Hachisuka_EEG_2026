function pval_pasbycorr = C2F_plotsuperandbasic_byPAS(corrbyPAS)
% Function used on main_behavior_script.m
numSubjects=30;
condColors = {[0.9290, 0.6940, 0.1250],[0.6350, 0.0780, 0.1840],[0.4660 0.6740 0.1880]};

figure; set(gcf,'Color','w');
hold on

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
groupLabels = {'PAS 0','PAS 1','PAS 2','PAS 3'};

for r = 1:3 
    subplot(1,3,r); hold on
    for superorbasic = 1:2 
        meanvals = squeeze(nanmean(corrbyPAS{superorbasic}(:,r, :).*100));
        semvals = squeeze(nanstd(corrbyPAS{superorbasic}(:,r, :).*100))./sqrt(numSubjects);

        for g = 1:4  
            xPos = repmat(g+2,1,numSubjects);
            % boxchart(xPos, corrbyPAS{superorbasic}(:,r, g).*100,'BoxWidth',0.2,'BoxFaceColor',condColors{r}, ...
            %     'JitterOutliers','off','MarkerStyle','none');  
            individDataPts = corrbyPAS{superorbasic}(:,r, g)'.*100;
            if superorbasic == 1
                s = swarmchart(xPos,individDataPts,30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', condColors{r}); %,'*','LineWidth',2,'Color','k');
            else
                s = swarmchart(xPos,individDataPts,30,'black','o', 'MarkerEdgeColor',condColors{r},'MarkerFaceColor', 'none'); %,'*','LineWidth',2,'Color','k');
            end
            s.XJitter = 'rand';
            s.XJitterWidth = 0.1;
            s.MarkerFaceAlpha = .7;

            pval_pasbycorr{superorbasic}(g,r) = signrank(individDataPts,50);
        end
    end
    for superorbasic = 1:2 
        meanvals = squeeze(nanmean(corrbyPAS{superorbasic}(:,r, :).*100));
        semvals = squeeze(nanstd(corrbyPAS{superorbasic}(:,r, :).*100))./sqrt(numSubjects);
        if superorbasic == 1
            errorbar(3:6,meanvals,semvals,'Color','k','LineWidth',1);
            plot(3:6,meanvals,'Color',condColors{r},'LineWidth',1);
        else
            errorbar(3:6,meanvals,semvals,'Color','k','LineWidth',1,'LineStyle','--');
            plot(3:6,meanvals,'Color',condColors{r},'LineWidth',1,'LineStyle','--');
        end
    end
    xlim([2 7]);                    
    xticks(2:7);                    
    xticklabels({'', '0', '1', '2', '3', ''});
    ylabel('% Correct');
    xlabel('PAS');
    set(gca, 'YGrid', 'on');
    yline(50,'--');
    set(gca,'FontSize',14);
    hold off
end
end
