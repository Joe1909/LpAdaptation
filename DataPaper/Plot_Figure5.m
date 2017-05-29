%plot 20D
dim = 20; 
Con = 3;

cols =[   0.8941    0.1020    0.1098
    0.2157    0.4941    0.7216
    0.3020    0.6863    0.2902
    0.5961    0.3059    0.6392
    1.0000    0.4980         0
    1.0000    1.0000    0.2000
    0.6510    0.3373    0.1569
    0.9686    0.5059    0.7490
    0.6000    0.6000    0.6000
    0.3216    0.3216    0.3216
         0         0         0];

styles_set2 = {'*-','d-','o-','s-'};
fontSize2 = 14; 
fontSizeAll = 18;


filepath2 = [pwd,'/'];
pnormVec = [0.5,1,2,101];
%% 

load('xVsTime_rep_dim20_numSteps40'); % run get_DataPlot_Lp to get this file; 

%%%%%%%%%%%%%%%%%%% normalized volume
figure;
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.

h=nan(4,1);
 for pf=1:4

    subplot(2,2,pf);
    set(gca,'FontSize',fontSizeAll);

    Vol_true = Vol_lp(dim,1,pnormVec(pf));
    
    for pt=1:4
        Vol_est = Vol_est_cell{pf,pt};
        Vol_est_normalized = Vol_est./Vol_true;
        stdVol_est = nanstd(Vol_est_normalized,[],3);
        meanVol_est = nanmean(Vol_est_normalized,3);

        xvalues = xvalues_cell{pf,pt};
        xvalues_m = nanmean(xvalues,2);
        xvalues_std = nanstd(xvalues,[],2);


        pp=2;
        h(pt) = plot(xvalues_m(1:2:end),meanVol_est(1:2:end,pp),styles_set2{pt},'Color',cols(pt,:),'LineWidth',2,'MarkerSize',10);
        hold on
        errorbar(xvalues_m(1:2:end),meanVol_est(1:2:end,pp),stdVol_est(1:2:end,pp),styles_set2{pt},'Color',cols(pt,:),'LineWidth',1,'MarkerSize',10);
        hold on

        xlim([0 max(xvalues(:))]);
   end
    if pnormVec(pf) == 101
        title(['feasible region: $L_\infty$-ball'],'FontSize',fontSizeAll,'Interpreter','latex');
    elseif pnormVec(pf) == 0.5
        title(['feasible region: $L_{',num2str(pnormVec(pf)),'}$ -ball'],'FontSize',fontSizeAll,'Interpreter','latex');
    else
        title(['feasible region: $L_{',num2str(pnormVec(pf)),'}$-ball'],'FontSize',fontSizeAll,'Interpreter','latex');
    end
    Vol_t = 1;
    fline=plot([0 max(xvalues(:))],[Vol_t Vol_t],'--k','LineWidth',2);

    if pf >=3
        xlabel('Evaluations','Interpreter','latex');
    else
        xlabel('');
    end
    if pf ==1 || pf == 3
        ylabel('Relative Volume','Interpreter','latex');
    else
        ylabel('');
    end
    set(gca,'FontSize',fontSizeAll);

    ypos = 1.6; 
    verticalLines = verticalLines_cell{pf,pt};
    verticalLines2 = mean(verticalLines,2);
    len_hitP = length(verticalLines2);
    for l=1:len_hitP
        plot([verticalLines2(l) verticalLines2(l)],[0 ypos],'Color',[.7 .7 .7]);
        hold on
        if l==1
            text(verticalLines2(l)/2,1.55,['hittP ',num2str(hitP_desired(l))],'BackgroundColor',[240/255 240/255 240/255],'HorizontalAlignment','center','FontSize',fontSizeAll);
        else
            text((verticalLines2(l)+verticalLines2(l-1))/2,1.55,['',num2str(hitP_desired(l))],'BackgroundColor',[240/255 240/255 240/255],'HorizontalAlignment','center','FontSize',fontSizeAll);
        end
    end
            
    ylim([0 ypos]);

 end

legend([h;fline],{'$L_{0.5}$ -proposal',...
            '$L_{1}$-proposal', '$L_{2}$-proposal','$L_{\infty}$-proposal','true'},'Orientation','horizontal','Position',[0.5 0.51 0 0],'FontSize',fontSizeAll,'Interpreter','latex');

    
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
set(gcf, 'Color', 'w');

