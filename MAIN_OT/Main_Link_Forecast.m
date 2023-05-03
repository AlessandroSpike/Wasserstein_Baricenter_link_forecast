clc; clearvars; close all
%% addpath
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\UTILITIES')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\DATI_PROCESSATI')
%% load dati
load DatiPerOT_QUARTILE_REVISION_TRADE.mat
%% parameter selection
MaxLag = 5; % scegli settore
NumObs =  numel(CountContainer{1});
% lunghezza serie storica
Lunghezza = size(CountContainer,1);
P_forecast_container = zeros(NumObs,Lunghezza,MaxLag);
TestLag_container= zeros(Lunghezza,MaxLag);
Lambda_container= cell(Lunghezza,MaxLag);
%% loop over different lags
for v = MaxLag +1 :Lunghezza
    %% costo
    C = Co; % CoGDP
    %% dependent variable
    Countreshape_y = reshape(CountContainer{v},NumObs,1);
    Countreshape_y(Countreshape_y==0)=10^-4;
    y = Countreshape_y./sum(Countreshape_y);
    %disp(['Variab Dip: ' ,num2str(anno(v))])
    for Lag = 1:MaxLag
        %% independent variables
        X = zeros(NumObs,Lag);
        for vv = 1:Lag
            posiz = v-vv;
            %disp(['Variab Indip: ' ,num2str(anno(posiz))])
            Countreshape = reshape(CountContainer{posiz},NumObs,1);
            Countreshape(Countreshape==0)=10^-4;
            % posiz 1 + recent 
            X(:,vv) = Countreshape./sum(Countreshape);    
        end
  
        %% parametri
        % num of  x histograms
        S = size(X,2);
        % num iter
        L = 50;
        % regularize param
        gamma = 8;
        % graadient descend alpha
        alphav = 1;
        % termination tolerance
        tol = 1e-5; 
        tol2= 1e-5;
        % maximum number of allowed iterations
        maxiter = 50;
        % minimum allowed perturbation
        dxmin = 1e-5;
        % initialize param
        lambdav = repmat(1/S,S,1);
        theta = Inf;
        %% Kernel approx of cost
        K = exp(-C/gamma);        
        KT = K';
        %% main GD
        if Lag==1
             lambdan=lambdav;
             [P,wn]=MainOT(X,y,size(C,1),K,KT,S,L,lambdav);
             disp('Solution Found')
        else
            [P,wv]=MainOT(X,y,size(C,1),K,KT,S,L,lambdav);
            lambdan = lambdav-alphav*wv;
            lambdan=exp(lambdan)/sum(exp(lambdan));
            [P,wn]=MainOT(X,y,size(C,1),K,KT,S,L,lambdan);
            dx = norm(lambdan-lambdav);
            gnorm = norm(wn-wv);
            disp(['Norma Gradiente: ',num2str(gnorm),' - Norma D-Lambda: ',num2str(dx)]);             
            if  gnorm<tol || dx<dxmin
                disp('Solution Found')
            else
                for t = 1:maxiter
                    alphan = min(sqrt(1+theta)*alphav,norm(lambdan-lambdav)/(2*norm(wn-wv)));
                    lambdav = lambdan;
                    wv = wn;
                    lambdan = lambdav-alphan*wv;
                    lambdan=exp(lambdan)/sum(exp(lambdan));
                    theta = alphan/alphav;
                    [P,wn]=MainOT(X,y,size(C,1),K,KT,S,L,lambdan);
                    gnorm = norm(wn-wv);
                    gnorm2 = norm(wn);
                    dx = norm(lambdan-lambdav);
                    disp(['Norma Gradiente: ',num2str(gnorm),' - Norma D-Lambda: ',num2str(dx)]);
                    if gnorm<tol || dx<dxmin
                        disp('Solution Found')
                        break;
                    end
                    if gnorm2<tol2 
                        disp('Solution Found')
                        break;
                    end
                end
            end
        end
        %% model selection
        WD = WassersteinDist(y,P,K,KT,C,maxiter);
        TestLag_container(v,Lag)=WD;
        Lambda_container{v,Lag}=lambdan;
         %% forecast
        X_update = [y,X(:,1:end-1)];
        P_forecast=MainOT_forecast(X_update,size(C,1),K,KT,S,L,lambdan);
        P_forecast_container(:,v,Lag)=P_forecast;
    end
end
%% save
save('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_OT\RisultatiOT_QUARTILE_REVISION_TRADE.mat','P_forecast_container','TestLag_container','Lambda_container')