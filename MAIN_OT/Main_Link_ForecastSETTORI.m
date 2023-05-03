clc; clearvars; close all
%% addpath
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\DATI_PROCESSATI')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\UTILITIES')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_ALTERNATIVES')
%% load dati
load DatiPerOT_SETTORI_PESATI_REVISION.mat
load ModelSelectionPESATI_REVISION.mat
%% parameter selection
MaxLag = 5; % scegli settore
NumObs =  numel(CountContainer{1});
% lunghezza serie storica
Lunghezza = size(CountContainer,1);
P_forecast_container = zeros(NumObs,Lunghezza,MaxLag);

%% loop over different lags
nik = 1;
for v = MaxLag +1 :Lunghezza
    %% costo
    C = Co; % CoGDP
    %% dependent variable
    Countreshape_y = reshape(CountContainer{v,8},NumObs,1);
    Countreshape_y(Countreshape_y==0)=10^-20;
    y = Countreshape_y./sum(Countreshape_y);
    
    for Lag = 1:MaxLag
        if Lag == posmin(nik)   
            %% independent variables
            X = zeros(NumObs,Lag+1);
            for vv = 1:Lag
                posiz = v-vv;
                %disp(['Variab Indip: ' ,num2str(anno(posiz))])
                Countreshape = reshape(CountContainer{posiz,8},NumObs,1);
                Countreshape(Countreshape==0)=10^-20;
                % posiz 1 + recent 
                X(:,vv) = Countreshape./sum(Countreshape);    
            end
            CountreshapeR = reshape(CountContainer{v-1,4},NumObs,1);
            CountreshapeR(CountreshapeR==0)=10^-20;
            X(:,end) = CountreshapeR./sum(CountreshapeR);  

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
            tol = 1e-5; % per discreti -7 per pesati -5
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
            K(K<10^-20)=10^-20;
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

            %% forecast
            X_update = [y,[X(:,1:end-2),X(:,1:end)]];
            P_forecast=MainOT_forecast(X_update,size(C,1),K,KT,S,L,lambdan);
            P_forecast_container(:,v,Lag)=P_forecast;
        end  
    end
    nik = nik+1;
end
%% save
save('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_OT\RisultatiOT_SETTORI_ITC_VS_ENERGY_LAG_PESATI_REVISION.mat','P_forecast_container')