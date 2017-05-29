
numSteps = 40; 
dim=20; 

clear outCell;
clear out;
strInFileName = ['_numSteps',num2str(numSteps)];
filepath = [pwd,'/'];

pnormVec = [0.5,1,2,101];
numPf = 4; 
numPt = 4; 

numRep = 10;

% save volume
Vol_est_cell = cell(numPf,numPt);
Vol_true = nan(numPf);

xvalues_cell = cell(numPf,numPt);
verticalLines_cell = cell(numPf,numPt);

rounds = 1; %to display progress
for pf=1:4
    
    allFiles = dir([filepath,['out_Lp_stretched20D_','pf',num2str(pf),'*.mat']]);
    load([filepath,allFiles(1).name]);%just load to get feasible volume
    
    out = outCell{1};
    targetC = out.opts.oracleInopts{4}^2*out.opts.oracleInopts{2}.^2;
    targetC_n = targetC./det(targetC)^(1/dim);
    r2_n = nthroot(det(targetC),dim*2);
    Vol_true(pf) = Vol_lp(dim,r2_n,pnormVec(pf));

    for pt=1:4
        file = dir([filepath,['out_Lp_stretched20D_','pf',num2str(pf),'_pt',num2str(pt),'.mat']]);
        disp(num2str(rounds)); 
        rounds = rounds+1;

        load([filepath,file(1).name]);
        out = outCell{1};
        
        len_hitP = length(out.opts.para_hitP_adapt.PVec);
        verticalLines = nan(len_hitP,numRep);

       
        [numRep,~]=size(outCell);


        maxEvalVec = nan(numRep,1);
        for rep=1:numRep
            maxEvalVec(rep) = outCell{rep}.cntVec(end);
        end
        maxE = max(maxEvalVec);

        numEnd = numSteps;% min(numSteps,ls);

        Vol_est = nan(numEnd,2,numRep);
        xvalues = nan(numEnd,numRep);

        for rep=1:numRep
                out = outCell{rep};
                
                %specify points, where approxVol etc should be evaluated
                idxVec_funEval = round(linspace(1,maxE,numSteps));
                cnt = out.cntVec;
                cnt(1)=1;
                idx_cnt=nan(numSteps,1);
                for k=1:numSteps
                    tmp= find(cnt>=idxVec_funEval(k),1,'first');
                    if ~isempty(tmp)
                        idx_cnt(k) = tmp;
                    end
                end
                tmp = find(isnan(idx_cnt),1,'first');
                if isempty(tmp)
                    l_idx_cnt = length(idx_cnt);
                else
                    l_idx_cnt = tmp-1;
                end
                int_ls = idx_cnt;

                hitP_desired(1:len_hitP) = out.opts.para_hitP_adapt.PVec;

                for v=1:l_idx_cnt%numEnd    
                    idx = int_ls(v);
                    Q = out.QCell{idx};
                    r = out.rVec(idx);

                    C = r^2*Q*Q';
                    if det(C) == 0
                        C_n = C./1e-4^(1/dim);
                    else
                        C_n = C./det(C)^(1/dim);
                    end
                    hittProb = out.P_empVecWindow(idx);
                    p = out.opts.pn;
                    Vol_est(v,1,rep) = hittProb * abs(det(Q))* Vol_lp(dim,r,p);

                    %take mean of part of rValues and hitting probabilities
                    %mean over at least 2*msize samples
                    msize = 5;
                    rm = mean(out.rVec(max(1,idx-msize):idx));
                    hittProbm = mean(out.P_empVecWindow(max(1,idx-msize):idx));
                    Vol_est(v,2,rep) = hittProbm * abs(det(Q)) * Vol_lp(dim,rm,p);
                end

                    xval=nan(numSteps,1);
                    xval(1:l_idx_cnt) = out.cntVec(idx_cnt(1:l_idx_cnt));
                    if idx_cnt(1)==1
                        xval(1)=1;
                    end
                    xvalues(:,rep) = xval;
                    len_h = length(out.adaptation.iterVec);
                    verticalLines(1:len_h,rep) = out.adaptation.iterVec;
        end
        
        xvalues_cell{pf,pt} = xvalues;
        verticalLines_cell{pf,pt}=verticalLines;
        Vol_est_cell{pf,pt} = Vol_est;

    end
end
clear outCell
save(['xVsTime_rep_dim',num2str(dim),strInFileName,'.mat'], 'Vol_est_cell','Vol_true','xvalues_cell','verticalLines_cell','hitP_desired');%,...
