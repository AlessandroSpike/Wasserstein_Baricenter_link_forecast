clc; clearvars; close all
%% addpath
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\DATI_PROCESSATI')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\UTILITIES')
%% load dati
load DatiPerOT_PESATI_REVISION_TRADE.mat
%% parameter selection
MaxLag = 5; % scegli settore
NumObs = size(reti,1)*size(reti,2);
% lunghezza serie storica
Lunghezza = size(reti,3);
%% contenitori
P_forecast_containerNMF_S = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerNMF2_S = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerK_S =zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerK2_S = zeros(NumObs,Lunghezza,MaxLag);

P_forecast_containerRWR_S = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerRWR2_S = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerAA_S = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerAA2_S = zeros(NumObs,Lunghezza,MaxLag);

%% loop over different lags
for v = MaxLag +1 :Lunghezza
    for Lag = 1:MaxLag
        %% reti
        X = zeros(size(reti,1),size(reti,2),Lag);
        for vv = 1:Lag
            posiz = v-vv;           
            % posiz 1 + recent 
            X(:,:,vv) = reti(:,:,posiz+1);          
        end
        X2 = X;
        %% aggrego
        if Lag>1           
            X = mean(X,3);          
            theta = 0.2;
            aus = zeros(Lag,1);
            for h = 1:Lag
                aus(h)= (1-theta)^(Lag-h);
                X2(:,:,h)=X2(:,:,h)*aus(h);
            end
            X2 = mean(X2,3);          
        end
        %% forecast  
        % NMF
        [W,H,D] = nnmf(X,1,'Algorithm','mult');
        forecasta = W*H;

        [W,H,D] = nnmf(X,2,'Algorithm','mult');
        forecastb = W*H;

        [W,H,D] = nnmf(X,3,'Algorithm','mult');
        forecastc = W*H;

        [W,H,D] = nnmf(X,4,'Algorithm','mult');
        forecastd = W*H;

        [W,H,D] = nnmf(X,5,'Algorithm','mult');
        forecaste = W*H;
        
         netN = forecasta+forecastb+...
                forecastc+forecastd+...
                forecaste;
            netN =netN/5;

        % KATZ
        beta = 0.1;
        S = (eye(size(X))-beta*X)^(-1)-eye(size(X));
        s = find(S<0);
        while isempty(s)==0
            beta = beta/5;
            S = (eye(size(X))-beta*X)^(-1)-eye(size(X));
            s = find(S<0);
        end

        % RWR
        net_prob=X./repmat(sum(X,2),1,size(X,1));
        net_prob(isnan(net_prob))=0;            
        forecast_rwa=inv(eye(size(net_prob))-0.1*net_prob)-diag(ones(size(net_prob,1),1));
        forecast_rwb=inv(eye(size(net_prob))-0.2*net_prob)-diag(ones(size(net_prob,1),1));
        forecast_rwc=inv(eye(size(net_prob))-0.6*net_prob)-diag(ones(size(net_prob,1),1));
        forecast_rwd=inv(eye(size(net_prob))-0.8*net_prob)-diag(ones(size(net_prob,1),1));
        forecast_rwe=inv(eye(size(net_prob))-0.9*net_prob)-diag(ones(size(net_prob,1),1));
        netN_rw = forecast_rwa+forecast_rwb+...
            forecast_rwc+forecast_rwd+...
            forecast_rwe;
        net_rwN =netN_rw/5;
        forecast_rw=net_rwN; 

        %AA
        netbin=X;
        netbin(netbin>0)=1;
        netaa=max(netbin,netbin');
        n = size(netaa,1);
        deg = sum(netaa,2);
        forecast_aa_aus = zeros(n);
        for i = 1:n
            for j = i+1:n
                commonNeighbors = (netaa(i,:) & netaa(j,:));
                forecast_aa_aus(i,j) = sum(1./log10(deg(commonNeighbors)));
                forecast_aa_aus(j,i) = forecast_aa_aus(i,j);
            end
        end
        forecast_aa=forecast_aa_aus;
 
        % CONTNITORI
        P_Forecast_NMF = reshape(netN,NumObs,1);
        P_forecast_containerNMF_S(:,v,Lag)=P_Forecast_NMF/sum(P_Forecast_NMF);
        P_Forecast_K = reshape(S,NumObs,1);
        P_forecast_containerK_S(:,v,Lag)=P_Forecast_K/sum(P_Forecast_K);        
        P_Forecast_RWR = reshape(forecast_rw,NumObs,1);
        P_forecast_containerRWR_S(:,v,Lag)=P_Forecast_RWR/sum(P_Forecast_RWR);
        P_Forecast_AA = reshape(forecast_aa,NumObs,1);
        P_forecast_containerAA_S(:,v,Lag)=P_Forecast_AA/sum(P_Forecast_AA);

        %NMF
        [W2,H2,D2] = nnmf(X2,1,'Algorithm','mult');
        forecast2a = W2*H2;

        [W,H,D2] = nnmf(X2,2,'Algorithm','mult');
        forecast2b = W*H;

        [W,H,D2] = nnmf(X2,3,'Algorithm','mult');
        forecast2c = W*H;

        [W,H,D2] = nnmf(X2,4,'Algorithm','mult');
        forecast2d = W*H;

        [W,H,D] = nnmf(X2,5,'Algorithm','mult');
        forecast2e = W*H;
    
          net2N = forecast2a+forecast2b+...
                forecast2c+forecast2d+...
                forecast2e;
            net2N =net2N/5;
        % KATZ
        beta2 = 0.1;
        S2 = (eye(size(X2))-beta*X2)^(-1)-eye(size(X2));
        s2 = find(S2<0);
        while isempty(s2)==0
            beta2 = beta2/5;
            S2 = (eye(size(X2))-beta*X2)^(-1)-eye(size(X2));
            s2 = find(S2<0);
        end

        % RWR
        net_prob=X2./repmat(sum(X2,2),1,size(X,1));
        net_prob(isnan(net_prob))=0;            
        forecast_rwa=inv(eye(size(net_prob))-0.1*net_prob)-diag(ones(size(net_prob,1),1));
        forecast_rwb=inv(eye(size(net_prob))-0.2*net_prob)-diag(ones(size(net_prob,1),1));
        forecast_rwc=inv(eye(size(net_prob))-0.6*net_prob)-diag(ones(size(net_prob,1),1));
        forecast_rwd=inv(eye(size(net_prob))-0.8*net_prob)-diag(ones(size(net_prob,1),1));
        forecast_rwe=inv(eye(size(net_prob))-0.9*net_prob)-diag(ones(size(net_prob,1),1));
        netN_rw = forecast_rwa+forecast_rwb+...
            forecast_rwc+forecast_rwd+...
            forecast_rwe;
        net_rwN =netN_rw/5;
        forecast_rw2=net_rwN; 

        %AA
        netbin=X2;
        netbin(netbin>0)=1;
        netaa=max(netbin,netbin');
        n = size(netaa,1);
        deg = sum(netaa,2);
        forecast_aa_aus = zeros(n);
        for i = 1:n
            for j = i+1:n
                commonNeighbors = (netaa(i,:) & netaa(j,:));
                forecast_aa_aus(i,j) = sum(1./log10(deg(commonNeighbors)));
                forecast_aa_aus(j,i) = forecast_aa_aus(i,j);
            end
        end
        forecast_aa2=forecast_aa_aus;
        
        % CONTENITORI
        P_Forecast_NMF2 = reshape(net2N,NumObs,1);
        P_forecast_containerNMF2_S(:,v,Lag)=P_Forecast_NMF2/sum(P_Forecast_NMF2);
        P_Forecast_K2 = reshape(S2,NumObs,1);
        P_forecast_containerK2_S(:,v,Lag)=P_Forecast_K2/sum(P_Forecast_K2);
        P_Forecast_RWR2 = reshape(forecast_rw2,NumObs,1);
        P_forecast_containerRWR2_S(:,v,Lag)=P_Forecast_RWR2/sum(P_Forecast_RWR2);
        P_Forecast_AA2 = reshape(forecast_aa2,NumObs,1);
        P_forecast_containerAA2_S(:,v,Lag)=P_Forecast_AA2/sum(P_Forecast_AA2);
        
    end
end
%% save
save('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_ALTERNATIVES\RisultatiNMF_PESATI_REVISION_TRADE.mat',...
    'P_forecast_containerNMF2_S','P_forecast_containerNMF_S','P_forecast_containerK_S','P_forecast_containerK2_S',...
    'P_forecast_containerRWR2_S','P_forecast_containerRWR_S','P_forecast_containerAA_S','P_forecast_containerAA2_S')