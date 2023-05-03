clc; clearvars; close all
%% addpath
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_OT')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\RISULTATI_ALTERNATIVES')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\UTILITIES')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\DATI_PROCESSATI')
%% load results
load RisultatiOT_PRCTILE_REVISION_GDP_TRADE
P_forecast_containerGDP = P_forecast_container;
load RisultatiOT_PRCTILE_REVISION_TRADE
load RisultatiNMF_PRCTILE_REVISION_TRADE
load('DatiPerOT_PRCTILE_REVISION_TRADE','CountContainer')
anni = 2003:2021;
%% colori
rgb = vals2colormap(size(TestLag_container,2)+1:size(TestLag_container,1), 'turbo');
%% show optimal lag
leg = num2str(anni(size(TestLag_container,2)+1:size(TestLag_container,1))');
leg1 = ["Value-";"Optimum-"];
leg1 = repmat(leg1,length(leg),1);
legx = string(repelem(leg,2,1));
leg2 = strcat(leg1,legx);
leg = string(leg);
posmin=zeros(length(size(TestLag_container,2)+1:size(TestLag_container,1)),1);
figure
x = 1;
for k = size(TestLag_container,2)+1:size(TestLag_container,1)
    plot(TestLag_container(k,:),'color',rgb(x,:),'linewidth',2)
    hold on
    [a,b]=min(TestLag_container(k,:));
    posmin(x)=b;
    plot(b, a, 'o','color',rgb(x,:), 'MarkerSize', 10);
    x = x+1;
end
axis tight
grid on
ylabel('Wasserstein Distance')
xlabel('Lag')
xticks(1:size(TestLag_container,2))
xticklabels(1:size(TestLag_container,2))
title('Model Selection')
set(gca,'fontsize',10,'Fontweight','bold')
lgd = legend(leg2,'location','eastoutside');
title(lgd,'Years')
%% plot forecasted network
Forecast = P_forecast_container(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
ForecastGDP = P_forecast_containerGDP(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastk = P_forecast_containerK(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastk2 = P_forecast_containerK2(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastnmf = P_forecast_containerNMF(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastnmf2 = P_forecast_containerNMF2(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);

Forecastrwr = P_forecast_containerRWR(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastrwr2 = P_forecast_containerRWR2(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastaa = P_forecast_containerAA(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);
Forecastaa2 = P_forecast_containerAA2(:,size(TestLag_container,2)+1:size(TestLag_container,1),:);


ForecastOT = zeros(size(Forecast,1),size(Forecast,2));
ForecastOT_GDP = zeros(size(Forecast,1),size(Forecast,2));
ForecastK = zeros(size(Forecast,1),size(Forecast,2));
ForecastK2 = zeros(size(Forecast,1),size(Forecast,2));
ForecastNMF = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastNMF2 = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastRWR = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastRWR2 = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastAA = zeros(size(Forecastnmf,1),size(Forecast,2));
ForecastAA2 = zeros(size(Forecastnmf,1),size(Forecast,2));

Reali = zeros(size(Forecastnmf,1),size(Forecast,2));

for k = 1:size(Forecast,2)
    ForecastOT(:,k)=Forecast(:,k,posmin(k));
    ForecastOT_GDP(:,k)=ForecastGDP(:,k,posmin(k));
    ForecastNMF(:,k)=Forecastnmf(:,k,posmin(k));
    ForecastNMF2(:,k)=Forecastnmf2(:,k,posmin(k));
    ForecastK(:,k)=Forecastk(:,k,posmin(k));
    ForecastK2(:,k)=Forecastk2(:,k,posmin(k));
    ForecastRWR(:,k)=Forecastrwr(:,k,posmin(k));
    ForecastRWR2(:,k)=Forecastrwr2(:,k,posmin(k));
    ForecastAA(:,k)=Forecastaa(:,k,posmin(k));
    ForecastAA2(:,k)=Forecastaa2(:,k,posmin(k));
end

%% non parametric
figure
tiledlayout(2,5, 'TileSpacing', 'compact')

nexttile
ausrocRis=zeros(size(ForecastOT,2)-1,1);
ausprRis=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastOT(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRis(k) = out.AUROC;
    ausprRis(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Wasserstein')


nexttile
ausrocRisNMF=zeros(size(ForecastOT,2)-1,1);
ausprRisNMF=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastNMF(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRisNMF(k) = out.AUROC;
    ausprRisNMF(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Aver. NMF')

nexttile
ausrocRisNMF2=zeros(size(ForecastOT,2)-1,1);
ausprRisNMF2=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastNMF2(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRisNMF2(k) = out.AUROC;
    ausprRisNMF2(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Exp. NMF')

nexttile
ausrocRis_GDP=zeros(size(ForecastOT,2)-1,1);
ausprRis_GDP=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastOT_GDP(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRis_GDP(k) = out.AUROC;
    ausprRis_GDP(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Wasserstein GDP')


nexttile
ausrocRisK=zeros(size(ForecastOT,2)-1,1);
ausprRisK=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastK(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRisK(k) = out.AUROC;
    ausprRisK(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Aver. Katz')

nexttile
ausrocRisK2=zeros(size(ForecastOT,2)-1,1);
ausprRisK2=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastK2(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRisK2(k) = out.AUROC;
    ausprRisK2(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Exp. Katz')


nexttile
ausrocRisRWR=zeros(size(ForecastOT,2)-1,1);
ausprRisRWR=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastRWR(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRisRWR(k) = out.AUROC;
    ausprRisRWR(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Aver. RWR')

nexttile
ausrocRisRWR2=zeros(size(ForecastOT,2)-1,1);
ausprRisRWR2=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastRWR2(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRisRWR2(k) = out.AUROC;
    ausprRisRWR2(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Exp. RWR')

nexttile
ausrocRisAA=zeros(size(ForecastOT,2)-1,1);
ausprRisAA=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastAA(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRisAA(k) = out.AUROC;
    ausprRisAA(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Aver. AA')

nexttile
ausrocRisAA2=zeros(size(ForecastOT,2)-1,1);
ausprRisAA2=zeros(size(ForecastOT,2)-1,1);
for k = 1:size(ForecastOT,2)-1
    fore = ForecastAA2(:,k);
    thresh = linspace(min(fore),max(fore),1000);
    realX = reshape(CountContainer{5+1+k},numel(CountContainer{5+1+k}),1);
    out = aurocpanel(realX,fore,thresh');
    ausrocRisAA2(k) = out.AUROC;
    ausprRisAA2(k) = out.AUPR;
    plot(out.fpr,out.tpr,...
        'color',rgb(k,:),'LineWidth',1.5);
    hold on  
end
plot(0:1/length(thresh):1,0:1/length(thresh):1,'k--')
lgd = legend([leg(2:end);"45 degree line"]);
title(lgd,'Years')
lgd.Layout.Tile = 'east';
axis square
grid on
xlabel('False Positive Rate','FontSize',10,'fontWeight','bold');
ylabel('True Positive Rate','FontSize',10,'fontWeight','bold');    
set(gca,'FontSize',10,'fontWeight','bold')
title('Exp. AA')

%%

figure
tiledlayout(1,2, 'TileSpacing', 'compact')

nexttile
plot(ausrocRis,'ko--','linewidth',2)
hold on
plot(ausrocRis_GDP,'kd--','linewidth',2)
plot(ausrocRisNMF,'bo--','linewidth',2)
plot(ausrocRisNMF2,'bd--','linewidth',2)
plot(ausrocRisK,'ro--','linewidth',2)
plot(ausrocRisK2,'rd--','linewidth',2)
plot(ausrocRisRWR,'go--','linewidth',2)
plot(ausrocRisRWR2,'gd--','linewidth',2)
plot(ausrocRisAA,'mo--','linewidth',2)
plot(ausrocRisAA2,'md--','linewidth',2)

title('Auroc Values')
xticks(1:length(ausrocRis))
xticklabels(leg(2:end))
set(gca,'FontSize',10,'fontWeight','bold')
axis tight
grid on
nexttile
plot(ausprRis,'ko--','linewidth',2)
hold on
plot(ausprRis_GDP,'kd--','linewidth',2)
plot(ausprRisNMF,'bo--','linewidth',2)
plot(ausprRisNMF2,'bd--','linewidth',2)
plot(ausprRisK,'ro--','linewidth',2)
plot(ausprRisK2,'rd--','linewidth',2)
plot(ausprRisRWR,'go--','linewidth',2)
plot(ausprRisRWR2,'gd--','linewidth',2)
plot(ausprRisAA,'mo--','linewidth',2)
plot(ausprRisAA2,'md--','linewidth',2)
title('Aupr Values')
xticks(1:length(ausrocRis))
xticklabels(leg(2:end))
set(gca,'FontSize',10,'fontWeight','bold')
axis tight
grid on

lgd = legend('Wasserstein','Wasserstein GDP','Aver. NMF',...
    'Exp. NMF','Aver. Katz','Exp. Katz','Aver. RWR','Exp. RWR',...
    'Aver. AA','Exp. AA');
title(lgd,'Models')
lgd.Layout.Tile = 'east';