function out = aurocpanel(dummy,evL,thresh)
%% Function to compute the confusion table and the AUROC
% Loops have been avoided by vectorizing the input. Bear in mind the size
% of the grid. Up to 100000 should be fine. From then on a loop might be
% more efficient. Trial and error...
%
% Pablo Rovira Kaltwasser
% Financial Stability Division
% Central Bank of Luxembourg
%
%% step 1: parameters and vectorize
NT     = size(thresh,1);            % determine the number of thresholds
realiz = dummy;                     % binary variables indicating the crises or precrises periods (depending on the H0)
pred   = evL;                       % prediction (indicator variable)
nans   = ~isnan(pred);              % determine the position of the available data
nobs   = size(realiz,1);            % number of observations
PT     = repmat(thresh',nobs,1);    % vectorize
Ptemp  = pred(nans);
Rtemp  = realiz(nans);
P      = repmat(Ptemp,1,NT);        % vectorize
R      = repmat(Rtemp,1,NT);        % vectorize

%% step 2: confusion matrix
Itp   = R.*double(P>=PT);
Itn   = (1-R).*double(P<PT); 
Ifp   = (1-R).*double(P>=PT);
Ifn   = R.*double(P<PT);

Prev  = sum(Rtemp)/length(Rtemp);                   % Prevalence
TP    = sum(Itp)';                                  % True +
TN    = sum(Itn)';                                  % True -
FP    = sum(Ifp)';                                  % False +
FN    = sum(Ifn)';                                  % False -
TPR   = TP./(TP+FN);                                % True + rate (sensitivity)
FNR   = 1-TPR;                                      % False - rate
FPR   = FP./(TN+FP);                                % False + rate
TNR   = 1-FPR;                                      % True - rate (specificity)
PPV   = TP./(TP+FP);                                % Positive Predictive Value
FDR   = FP./(TP+FP);                                % False Discovery Rate
FOR   = FN./(TN+FN);                                % False Omission Rate
NPV   = TN./(TN+FN);                                % Negative Prediction Value
NS    = FPR./TPR;                                   % Noise to signal ratio
ACC   = (TP+TN)./(TP+TN+FP+FN);                     % Accuracy
CondP = (TPR*Prev)./(TPR*Prev+(1-TNR)*(1-Prev));    % Conditional Probablitiy

%% step 3: compute the area under the ROC curve (trapezoidal approximation)
Ylinks  = TPR(2:end,1);    Yrechts = TPR(1:end-1,1);
Xlinks  = FPR(2:end,1);    Xrechts = FPR(1:end-1,1);
auc     = sum((Xrechts-Xlinks).*(Yrechts+Ylinks)/2);        % area under the curve (trapezoidal approx) 
out.AUROC = auc;

nnc    = sum(Rtemp==0);  % number of no crisis periods
nc     = sum(Rtemp);     % number of crisis periods
auc2   = auc^2;
q1     = auc/(2-auc);
q2     = (2*(auc^2))/(1+auc);
varauc = (auc*(1-auc)+(nc-1)*(q1-auc2)+(nnc-1)*(q2-auc2))/(nc*nnc);
stdauc = sqrt(varauc);
lbauc  = auc - 1.96*stdauc;
ubauc  = auc + 1.96*stdauc;

dist  = auc - 0.5;
zstat = dist/stdauc;        % z-score
df    = nnc+nc;             % degrees of freedom
pval  = tpdf(zstat,df);     % p-value

out.stdauroc = stdauc;
out.lbauroc  = lbauc;
out.ubouroc  = ubauc;
out.pcalroc   = pval;

%% step 3: compute the area under the PR curve (trapezoidal approximation)
Ylinks  = PPV(2:end,1);    Yrechts = PPV(1:end-1,1);
Xlinks  = TPR(2:end,1);    Xrechts = TPR(1:end-1,1);
auc     = sum((Xrechts-Xlinks).*(Yrechts+Ylinks)/2);        % area under the curve (trapezoidal approx) 
out.AUPR = auc;

nnc    = sum(Rtemp==0);  % number of no crisis periods
nc     = sum(Rtemp);     % number of crisis periods
auc2   = auc^2;
q1     = auc/(2-auc);
q2     = (2*(auc^2))/(1+auc);
varauc = (auc*(1-auc)+(nc-1)*(q1-auc2)+(nnc-1)*(q2-auc2))/(nc*nnc);
stdauc = sqrt(varauc);
lbauc  = auc - 1.96*stdauc;
ubauc  = auc + 1.96*stdauc;

dist  = auc - 0.5;
zstat = dist/stdauc;        % z-score
df    = nnc+nc;             % degrees of freedom
pval  = tpdf(zstat,df);     % p-value

out.stdaupr = stdauc;
out.lbaupr  = lbauc;
out.uboupr  = ubauc;
out.pcalpr   = pval;

%% step 4: prepare the output
out.prev = Prev;
out.tp   = TP;      out.tn  = TN;   out.fp   = FP;      out.fn    = FN;
out.tpr  = TPR;     out.fnr = FNR;  out.fpr  = FPR;     out.tnr   = TNR;
out.ppv  = PPV;     out.fdr = FDR;  out.for  = FOR;     out.npv   = NPV;
out.ns   = NS;      out.acc = ACC;  out.condp = CondP;
out.ITP  = Itp;     out.ITN = Itn;  out.IFP  = Ifp;     out.IFN   = Ifn;
out.thresh = thresh;

%% step % graphical output
% figure()
%     p = plot(FPR,TPR,'b',0:1/NT:1,0:1/NT:1,'r--');
%     set(p(1),'LineWidth',1.5,'LineWidth',2);
%     xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
%     ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');
%     legend(['ROC curve ' , num2str(auc)]);
%     set(gca,'FontSize',10,'fontWeight','bold')



