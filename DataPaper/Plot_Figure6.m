
fontSizeAll = 22;

trueScale = 1; 

params_names = {'$v_1$';'$v_{12}$';'$v_{32}$';'$v_{132}$';'$v_{532}$';'$v_{13}$';'$v_{33}$';'$v_{133}$';'$v_{533}$'};%params_names = {'v1';'v12';'v32';'v132';'v532';'v13';'v33';'v133';'v533'};

%one solution found by storn after 22344 evaluations of transfer function
v1	= 0.000427; 
v12 = 12.775360; 
v32 = 2.853888;
v132 = 15.606793; 
v532 = 12.293596; 
v13 = 6.117639; 
v33 = 0.690770; 
v133 = 6.829826; 
v533 = 0.574528;

paramsDE = [v1;v12;v32;v132;v532;v13;v33;v133;v533];
paramsDE1 = log(paramsDE);

load('out_SCfilter.mat');        
numRep=length(outCell);
[~,dim] = size(outCell{1}.muVec);

xAccAll = nan(10*90000,dim);

xMuLast = nan(10,dim,numRep);        
xStart = nan(numRep,dim);
        
xMuRunCell = cell(numRep,1);

every = 10;
idx_old = 1; 
for rep=1:numRep
    out = outCell{rep};
    [len,~] = size(out.xAcc(1:every:end,:));
    idx = idx_old+len-1;
    xAccAll(idx_old:idx,:) = out.xAcc(1:every:end,:);
    idx_old = idx+1; 
    
    xStart(rep,:) = out.xAcc(1,:);
    xMuRunCell{rep} = out.muVec;
    
end
      
if trueScale==1
    xAccAll = exp(1).^xAccAll;
    xStart = exp(1).^xStart;
    paramsDE1 = exp(1).^paramsDE1;
    paramsDE = exp(1).^paramsDE;

end
xData = xAccAll(1:idx,:);
xData = unique(xData,'rows');
[~,N]=size(xData);

if size(xStart,2)~= N
    error('xMu and xData need to have same dimension');
end

numXmu = size(xStart,1);

vname = [1;12;32;132;532;13;33;133;533];
figure('Position', [100, 100, 900, 900]);
maxlevels =0;
idx_i = [2,4,6,8];
idx_j=[3,5,7,9];

for k=1:4
    i = idx_i(k);
    j = idx_j(k);

    x=xData(:,i);y=xData(:,j);
    n=49;
    xi = linspace(min(x(:)),max(x(:)),n);
    yi = linspace(min(y(:)),max(y(:)),n);
    xr = interp1(xi,1:numel(xi),x,'nearest')';
    yr = interp1(yi,1:numel(yi),y,'nearest')';

    Z = accumarray([xr' yr'],1,[n n]);
    maxlevels = max(maxlevels,max(max(Z)));

end

levels=[1:15:maxlevels]';

for k=1:4
    i = idx_i(k);
    j = idx_j(k);


    % Plot the data for the current pair of dimensions
    subplot(2,2,k);
    
    %http://blogs.mathworks.com/videos/2010/01/22/advanced-making-a-2d-or-3d-histogram-to-visualize-data-density/
    x=xData(:,i);y=xData(:,j);
    n=49;

    xi = linspace(min(x(:)),max(x(:)),n);
    yi = linspace(min(y(:)),max(y(:)),n);
    xr = interp1(xi,1:numel(xi),x,'nearest')';
    yr = interp1(yi,1:numel(yi),y,'nearest')';

    Z = accumarray([xr' yr'],1,[n n]);
    contour(xi,yi,Z',levels,'LineWidth',1);
    grid on
    hold on; 
    %trajectory of run 5 (run 5 and 9 have biggest volume)
    for rep=5%[5,9]%1:numRep
        xMuRun1 = xMuRunCell{rep};
        if trueScale==1
            xMuRun1 = exp(1).^xMuRun1; 
        end
        plot(xMuRun1([1:20:100,200:200:end],i),xMuRun1([1:20:100,200:200:end],j),'-*k','linew',1,'Markersize',7);
    end

    plot(paramsDE1(i),paramsDE1(j),'h','color',[0.3,0.3,0.3],'markersize',15,'linew',3)

    xlim([min(x(:)) max(x(:))]);
    ylim([min(y(:)) max(y(:))]);
    %axis square
    set(gca,'fontsi',fontSizeAll);

    xlabel([params_names{i}],'fontsi',fontSizeAll,'interpreter','latex');%xlabel(['D',num2str(i)],'fontsi',18);
    ylabel([params_names{j}],'fontsi',fontSizeAll,'interpreter','latex');%ylabel(['D',num2str(j)],'fontsi',18);
    xlabh = get(gca,'XLabel');
    if k==1 
        set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0]);
    else
        set(xlabh,'Position',get(xlabh,'Position') - [0 .05 0]);
    end


end
  