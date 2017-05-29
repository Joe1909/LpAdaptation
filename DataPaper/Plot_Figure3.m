%close all
clear
load('out_Storn10')

%% colours
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

%% 2D example from Storn1999     
xConstraints = [0,5,5,10,10,5,5,0,0];
yConstraints = [-0.1,-0.1,-10,-10,10,10,0.1,0.1,-0.1];    
dim=2;

iterSave = length(out.rVec);
pn = out.opts.pn;
cnt_old = 1;

figure('Position', [0, 0,1200,1000]); 
out.cntVec(1)=1;

iterVec = [1,61,131,iterSave];

plot(xConstraints,yConstraints,'k','linewi',3);
hold on
fill(xConstraints,yConstraints,[0.9,0.9,0.9])
xlim([-12 14]);

hold on;
axis equal
for ii = 1:3
    iter = iterVec(ii);
    cc = ii;    

    eval_i = out.cntVec(iter);
    cnt = find(out.cntAcc <= eval_i,1,'last');
    cnt_old = cnt; 
    mu_tmp = out.muVec(iter,:);
    Q_tmp = out.QCell{iter};
    r_tmp = out.rVec(iter);

    if ii>1
      arrow(mu_old,mu_tmp,'linew',1.5,'length',17,'tipangle',11);%'BaseAngle',10,
    end
    mu_old = mu_tmp; 
    ball=lpBall2surface(7000*dim,2,pn,r_tmp,mu_tmp',Q_tmp,0);
    plot(ball(:,1),ball(:,2),'.','col',cols(ii,:),'markerSize',6);
    set(gca,'fontsize',18);
    
    m = plot(mu_tmp(1),mu_tmp(2),'x','Color',cols(cc,:),'MarkerSize',17,'linew',3);%m = plot(mu_tmp(1),mu_tmp(2),'bx','MarkerSize',15,'linew',3);
    m = plot(mu_tmp(1),mu_tmp(2),'x','Color','k','MarkerSize',17,'linew',1);%m = plot(mu_tmp(1),mu_tmp(2),'bx','MarkerSize',15,'linew',3); 
  
end
hold on;
numIterAv = 100; 

outLast.QCell = out.QCell(1:iterVec(end));
outLast.rVec = out.rVec(1:iterVec(end));
outLast.P_empVecAll = out.P_empVecAll(1:iterVec(end));


outCov=averageCov(out,numIterAv);
C_p = outCov.C_previous; 
% normalize C_p, such that det=1
C_normp = C_p./det(C_p)^(1/dim);
[eigvec,eigval]=eig(C_p);

 muVecLast = out.muVec(1:iterVec(end),:);
 mu = mean(muVecLast(end-numIterAv:end,:));

h = error_ellipse('mu',mu,'C',C_p,'style','m');

m = plot(mu(1),mu(2),'x','Color',cols(ii+1,:),'MarkerSize',17,'linew',3);
m = plot(mu(1),mu(2),'x','Color','k','MarkerSize',17,'linew',1);
arrow(mu_old,mu,'linew',1.5,'length',17,'tipangle',11);

xlabel('x_1');
ylabel('x_2');

 %%%%hitting probability
     axes('position', [0.26 0.67 0.23 0.2]);

    for ii = 1:4
        iter = iterVec(ii);
        eval_i =  out.cntVec(iter);
        cc=ii;
        plot([eval_i eval_i],[0 1],'Color',cols(cc,:),'linewidth',3);
        hold on
        
    end
    plot([out.cntVec(iterVec(end)-numIterAv),out.cntVec(iterVec(end)-numIterAv)],[0 1],'Color',cols(ii,:),'linewidth',3);

    fill([out.cntVec(iterVec(end)-numIterAv),out.cntVec(iterVec(end)),out.cntVec(iterVec(end)),out.cntVec(iterVec(end)-numIterAv)],[1,1,0,0],cols(ii,:),'facealpha',0.3);
     
    len = length(out.P_empVecAll);
    smooth_hitP = nan(len,1);
    Vol_est = nan(len,1);
    
    window = 50;
    
    for k=1:len
        eval = out.cntVec(k);
        tmp = max(k-window,1);
        low = out.cntVec(tmp);
        evalSum = eval - low + 1;  

        numFeas = sum((out.cntAcc <= eval) & (out.cntAcc >= low));
        smooth_hitP(k) = numFeas/evalSum;

        rm = mean(out.rVec(max(1,k-window):k));
        hittProbm = mean(smooth_hitP(max(1,k-window):k));
        
        Vol_est(k) = hittProbm * Vol_lp(dim,rm,pn);

    end

    
    plot(out.cntVec(1:iterVec(end)),smooth_hitP(1:iterVec(end)),'k','linewidth',2);%'color',[120,120,120]./256,'linewidth',1);
    hold on
    eval_i =  out.cntVec(1);
    cc = 1; 
    plot([eval_i eval_i],[0 1],'Color',cols(cc,:),'linewidth',3);
    
    plot([0 out.cntVec(iterVec(end))],[0.35 0.35],'k--','linew',3);
    
    xlim([-10 out.cntVec(iterVec(end))+10]);
    ylim([0 1]);
    xlabel('Evaluations','fontSize',18);
    ylabel('Hitting Probability','fontSize',18);
    set(gca,'fontsize',18);
    
    
    %volume
   axes('position', [0.26 0.2 0.23 0.2]);%axes('position', [0.25 0.2 0.26 0.2]);
%     grid on
    xlim([-20 out.cntVec(iterVec(end))+20]);
    max_Vol = max(out.volVec.*out.P_empVecAll);
    
    
    for ii = 1:4
        iter = iterVec(ii);
        eval_i =  out.cntVec(iter);
        cc=ii;
        plot([eval_i eval_i],[0 max_Vol],'Color',cols(cc,:),'linewidth',3);
        hold on
        
    end
    plot([out.cntVec(iterVec(end)-numIterAv),out.cntVec(iterVec(end)-numIterAv)],[0 max_Vol],'Color',cols(ii,:),'linewidth',3);
    fill([out.cntVec(iterVec(end)-numIterAv),out.cntVec(iterVec(end)),out.cntVec(iterVec(end)),out.cntVec(iterVec(end)-numIterAv)],[max_Vol,max_Vol,0,0],cols(ii,:),'facealpha',0.3);
   
    vol = out.volVec.*smooth_hitP;
    
    plot(out.cntVec(1:iterVec(end)),vol(1:iterVec(end)),'k','linewidth',2);
    plot([0 out.cntVec(iterVec(end))],[101 101],'k--','linewidth',3);
    
 
    ylim([0 max_Vol]);
    xlim([-10 out.cntVec(iterVec(end))+10]);
    xlabel('Evaluations','fontSize',18);
    ylabel('Volume Estimation','fontSize',18);
    set(gca,'fontsize',18);        
    set(gcf, 'Color', 'w');

