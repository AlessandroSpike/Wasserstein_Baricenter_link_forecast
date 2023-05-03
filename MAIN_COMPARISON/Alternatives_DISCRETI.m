clc; clearvars; close all
%% addpath
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\DATI_PROCESSATI')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\UTILITIES')
%% load dati
load DatiPerOT_PRCTILE_REVISION_TRADE.mat
%% parameter selection
MaxLag = 5; % scegli settore
NumObs =  numel(CountContainer{1});
% lunghezza serie storica
Lunghezza = size(CountContainer,1);
%% contenitori
P_forecast_containerNMF = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerNMF2 = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerK = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerK2 = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerRWR = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerRWR2 = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerAA = zeros(NumObs,Lunghezza,MaxLag);
P_forecast_containerAA2 = zeros(NumObs,Lunghezza,MaxLag);
%% creo reti
reti = zeros(length(NodiTutti),length(NodiTutti),length(LinkContainer));
nodiId = (NodiDtutti);
for t = 1:length(LinkContainer)
    link =  LinkContainer{t};
    for k = 1:size(link,1)
        cap = link(k,1);
        s = link(k,2);
        d = link(k,3);
        s2 = find(nodiId==s);
        d2 = find(nodiId==d);
        reti(s2,d2,t)=cap;
    end
end
%% loop over different lags
for v = MaxLag +1 :Lunghezza
    v
    for Lag = 1:MaxLag
        %% reti
        X = zeros(length(NodiTutti),length(NodiTutti),Lag);
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
        [a,b,val] = find(X>0);
        SS = NodiDtutti(a);
        DD = NodiDtutti(b);
        X_l = [X(val),SS,DD];
        [count,centri,pos] = histonet(X_l,Cap_bin2,NodiDtutti);
        forecast = zeros(size(count));
        forecastK = zeros(size(count));
        forecast_rw = zeros(size(count));
        forecast_aa = zeros(size(count));
        for r = 1:size(forecast,1)
            %NMF
            net = reshape(count(r,:,:),size(forecast,2),size(forecast,3));
            [W,H,D] = nnmf(net,1,'Algorithm','mult');
            forecasta = W*H;
    
            [W,H,D] = nnmf(net,2,'Algorithm','mult');
            forecastb = W*H;
    
            [W,H,D] = nnmf(net,3,'Algorithm','mult');
            forecastc = W*H;
    
            [W,H,D] = nnmf(net,4,'Algorithm','mult');
            forecastd = W*H;
    
            [W,H,D] = nnmf(net,5,'Algorithm','mult');
            forecaste = W*H;
            
             netN = forecasta+forecastb+...
                forecastc+forecastd+...
                forecaste;
            netN =netN/5;
            forecast(r,:,:)=netN;
            %KATZ
            beta = 0.1;
            S = (eye(size(net))-beta*net)^(-1)-eye(size(net));
            s = find(S<0);
            while isempty(s)==0
                beta = beta/5;
                S = (eye(size(net))-beta*net)^(-1)-eye(size(net));
                s = find(S<0);
            end
            forecastK(r,:,:)=S;
            % RW
            net_prob=net./repmat(sum(net,2),1,size(net,1));
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
            forecast_rw(r,:,:)=net_rwN;

            %AA
            netbin=net;
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
            forecast_aa(r,:,:)=forecast_aa_aus;
        end
        P_Forecast_NMF = reshape(forecast,NumObs,1);
        P_Forecast_NMF = P_Forecast_NMF./sum(P_Forecast_NMF);
        P_forecast_containerNMF(:,v,Lag)=P_Forecast_NMF;

        P_Forecast_RWR = reshape(forecast_rw,NumObs,1);
        P_Forecast_RWR = P_Forecast_RWR./sum(P_Forecast_RWR);
        P_forecast_containerRWR(:,v,Lag)=P_Forecast_RWR;

        P_Forecast_AA = reshape(forecast_aa,NumObs,1);
        P_Forecast_AA = P_Forecast_AA./sum(P_Forecast_AA);
        P_forecast_containerAA(:,v,Lag)=P_Forecast_AA;

        P_Forecast_K = reshape(forecastK,NumObs,1);
        P_Forecast_K = P_Forecast_K./sum(P_Forecast_K);
        P_forecast_containerK(:,v,Lag)=P_Forecast_K;


        [a2,b2,val2] = find(X2>0);
        SS2 = NodiDtutti(a2);
        DD2 = NodiDtutti(b2);
        X_l2 = [X2(val2),SS2,DD2];
        [count2,centri2,pos] = histonet(X_l2,Cap_bin2,NodiDtutti);
        forecast2 = zeros(size(count2));
        forecastK2 = zeros(size(count2));
        forecast_rw2 = zeros(size(count2));
        forecast_aa2 = zeros(size(count2));
        for r = 1:size(forecast2,1)
            net = reshape(count2(r,:,:),size(forecast2,2),size(forecast2,3));
            [W2,H2,D2] = nnmf(net,1,'Algorithm','mult');
            forecast2a = W2*H2;
    
            [W,H,D2] = nnmf(net,2,'Algorithm','mult');
            forecast2b = W*H;
    
            [W,H,D2] = nnmf(net,3,'Algorithm','mult');
            forecast2c = W*H;
    
            [W,H,D2] = nnmf(net,4,'Algorithm','mult');
            forecast2d = W*H;
    
            [W,H,D] = nnmf(net,5,'Algorithm','mult');
            forecast2e = W*H;
        
            netN2 = forecast2a+forecast2b+...
                forecast2c+forecast2d+...
                forecast2e;
            netN2 =netN2/5;
             forecast2(r,:,:)=netN2;

            beta2 = 0.1;
            S2 = (eye(size(net))-beta2*net)^(-1)-eye(size(net));
            s2 = find(S2<0);
            while isempty(s2)==0
                beta2 = beta2/5;
                S2 = (eye(size(net))-beta2*net)^(-1)-eye(size(net));
                s2 = find(S2<0);
            end
            forecastK2(r,:,:)=S2;

            % RW
            net_prob=net./repmat(sum(net,2),1,size(net,1));
            net_prob(isnan(net_prob))=0;            
            forecast_rw2a=inv(eye(size(net_prob))-0.1*net_prob);
            forecast_rw2b=inv(eye(size(net_prob))-0.2*net_prob);
            forecast_rw2c=inv(eye(size(net_prob))-0.6*net_prob);
            forecast_rw2d=inv(eye(size(net_prob))-0.8*net_prob);
            forecast_rw2e=inv(eye(size(net_prob))-0.9*net_prob);
            netN_rw2 = forecast_rw2a+forecast_rw2b+...
                forecast_rw2c+forecast_rw2d+...
                forecast_rw2e;
            net_rwN2 =netN_rw2/5;
            forecast_rw2(r,:,:)=net_rwN2;

            %AA
             netbin=net;
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
            forecast_aa2(r,:,:)=forecast_aa_aus;

        end

        % conteiner
        P_Forecast_NMF2 = reshape(forecast2,NumObs,1);
        P_Forecast_NMF2 = P_Forecast_NMF2./sum(P_Forecast_NMF2);
        P_forecast_containerNMF2(:,v,Lag)=P_Forecast_NMF2;

        P_Forecast_RWR2 = reshape(forecast_rw2,NumObs,1);
        P_Forecast_RWR2 = P_Forecast_RWR2./sum(P_Forecast_RWR2);
        P_forecast_containerRWR2(:,v,Lag)=P_Forecast_RWR2;

        P_Forecast_AA2 = reshape(forecast_aa2,NumObs,1);
        P_Forecast_AA2 = P_Forecast_AA2./sum(P_Forecast_AA2);
        P_forecast_containerAA2(:,v,Lag)=P_Forecast_AA2;

        P_Forecast_K2 = reshape(forecastK2,NumObs,1);
        P_Forecast_K2 = P_Forecast_K2./sum(P_Forecast_K2);
        P_forecast_containerK2(:,v,Lag)=P_Forecast_K2;
     
    end
end
%% save
save('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_ALTERNATIVES\RisultatiNMF_PRCTILE_REVISION_TRADE.mat','P_forecast_containerNMF2','P_forecast_containerNMF',...
    'P_forecast_containerK','P_forecast_containerK2',...
    'P_forecast_containerRWR','P_forecast_containerRWR2',...
    'P_forecast_containerAA','P_forecast_containerAA2')