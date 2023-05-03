clc; clearvars; close all
%% addpath
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\DATI_PROCESSATI')
addpath('C:\Users\Alessandro\Desktop\Lavori2023\ROYAL\UTILITIES')
%% load dati
load('DatiPerOT_QUARTILE_REVISION_TRADE.mat','NodiDtutti','NodiTutti','LinkContainer','Ltot','CountContainer','centri')
anni = 2003:2021;
%% creo reti
reti = zeros(length(NodiTutti),length(NodiTutti),length(LinkContainer));
nodiId = NodiDtutti;
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
%% plot disaggregato
periodo = 1;
Count1 = CountContainer{periodo};
figure
subplot(1,2,1)
imagesc(reti(:,:,periodo))
axis square
title(['Weighted Adj Matrix: ',num2str(anni(periodo))])
c = colorbar('southoutside');
colormap turbo
c.Label.String = 'Weights';
set(gca,'fontsize',10,'Fontweight','bold')
xticks(1:length(NodiTutti))
xticklabels(NodiTutti)
yticks(1:length(NodiTutti))
yticklabels(NodiTutti)

subplot(1,2,2)
netL = zeros(size(Count1,2),size(Count1,3),size(Count1,1));
for k = 1:size(Count1,1)
    netL(:,:,k)=reshape(Count1(k,:,:),size(Count1,2),size(Count1,3));
end
for k = 1:size(Count1,1)
    H(k) = slice(netL,[],[], k); %slice() requires at least 2x2x2
    set(H(k),'EdgeColor',[220/255	220/255	220/255]) %required so image isn't just an edge
    set(H(k),'FaceAlpha',0.8)
    hold on
    colormap(1- gray)
end
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('Countries')
ylabel('Countries')
zlabel('Capital invested (QUARTILES)')
zticks(1:length(centri{1}))
zticklabels(centri{1})
title('Multilayer Adj')
set(gca,'fontsize',10,'Fontweight','bold')
%% in-out deg e strength
Ind = zeros(size(reti,1),size(reti,3));
Outd = zeros(size(reti,1),size(reti,3));
Ins = zeros(size(reti,1),size(reti,3));
Outs = zeros(size(reti,1),size(reti,3));
TotFlow = zeros(size(reti,3),1);
NumLink = zeros(size(reti,3),1);

for t = 1:size(reti,3)
    r = reti(:,:,t);
    rb = r;
    rb(rb>0)=1;
    Ind(:,t) = sum(rb)';
    Ins(:,t) = sum(r)';
    Outd(:,t) = sum(rb,2);
    Outs(:,t) = sum(r,2);
    TotFlow(t) = sum(sum(r));
    NumLink(t) = sum(sum(rb));
end
%% plot mappe

Periodo = [1 6 12 18];
figure
for z = 1:length(Periodo)
    periodo = Periodo(z)
    link1 = LinkContainer{periodo};
    cap = link1(:,1);
    s = link1(:,2);
    d = link1(:,3);
    [LiS,LocS] = ismember(s,NodiDtutti);
    [LiD,LocD] = ismember(d,NodiDtutti);
    PosD = Ltot(LocD,:);
    PosS = Ltot(LocS,:);
    
    LatLong = Ltot;
    sizeV = Ins(:,periodo);
    colV = Outs(:,periodo);
    tolgo = intersect(find(sizeV==0),find(colV==0));
    if isempty(tolgo)==0
        LatLong(tolgo,:)=[];
        sizeV(tolgo)=[];
        colV(tolgo)=[];
    end
    rgb = vals2colormap(colV, 'turbo');
    
    tolgo2 = find(cap<prctile(cap,60));
    PosD(tolgo2,:)=[];
    PosS(tolgo2,:)=[];
    link1(tolgo2,:)=[];
    subplot(2,2,z)
    for k = 1:size(link1,1)
        geoplot([PosD(k,1),PosS(k,1)],[PosD(k,2),PosS(k,2)],'k:','LineWidth',0.1)
        hold on
    end
    geoscatter(LatLong(:,1),LatLong(:,2),800*(sizeV/sum(sizeV)),rgb,'filled')
    title(['WTN: ',num2str(anni(periodo))])
    geobasemap grayland
end
%% figure
figure
yyaxis right
plot(NumLink,'ko-','LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
ylabel('Number of Links')
yyaxis left
bar(TotFlow,'FaceColor',[.2 .6 .5],'EdgeAlpha',0.1,'FaceAlpha',.7)
ylabel('Total Flow')
axis tight
grid on
legend('Flow','Link')
xticks(1:size(Ind,2))
xticklabels(2003:2021)
title('Network statistics')
set(gca,'fontsize',10,'Fontweight','bold')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

figure
subplot(1,2,1)
b= bar3h(Ind',0.5);
colormap hot           % Choose a colormap
for k = 1:length(b)
    ydata = b(k).YData;               % Use YData property to create color gradient
    b(k).CData = ydata;               % Set CData property to Ydata
    b(k).FaceColor = "interp";        % Set the FaceColor to 'interp' to enable the gradient 
end
view(-23.950731707317058,29.373883490218486)
c=colorbar('southoutside');
c.Label.String = 'In-degree';
xticks(1:size(Ind,1))
xticklabels(NodiTutti)

zticks(1:size(Ind,2))
zticklabels(2003:2021)
title('In-Degree')
set(gca,'fontsize',10,'Fontweight','bold')

subplot(1,2,2)
b= bar3h(Outd',0.5);
colormap hot           % Choose a colormap
for k = 1:length(b)
    ydata = b(k).YData;               % Use YData property to create color gradient
    b(k).CData = ydata;               % Set CData property to Ydata
    b(k).FaceColor = "interp";        % Set the FaceColor to 'interp' to enable the gradient 
end
view(-23.950731707317058,29.373883490218486)
c=colorbar('southoutside');
c.Label.String = 'Out-degree';
xticks(1:size(Ind,1))
xticklabels(NodiTutti)

zticks(1:size(Ind,2))
zticklabels(2003:2021)
title('Out-degree')
set(gca,'fontsize',10,'Fontweight','bold')

figure

subplot(1,2,1)
b= bar3h(Ins',0.5);
colormap winter           % Choose a colormap
for k = 1:length(b)
    ydata = b(k).YData;               % Use YData property to create color gradient
    b(k).CData = ydata;               % Set CData property to Ydata
    b(k).FaceColor = "interp";        % Set the FaceColor to 'interp' to enable the gradient 
end
view(-23.950731707317058,29.373883490218486)
c=colorbar('southoutside');
c.Label.String = 'In-strength';
xticks(1:size(Ind,1))
xticklabels(NodiTutti)

zticks(1:size(Ind,2))
zticklabels(2003:2021)
title('In-strength')
set(gca,'fontsize',10,'Fontweight','bold')

subplot(1,2,2)
b= bar3h(Outs',0.5);
colormap winter           % Choose a colormap
for k = 1:length(b)
    ydata = b(k).YData;               % Use YData property to create color gradient
    b(k).CData = ydata;               % Set CData property to Ydata
    b(k).FaceColor = "interp";        % Set the FaceColor to 'interp' to enable the gradient 
end
view(-23.950731707317058,29.373883490218486)
c=colorbar('southoutside');
c.Label.String = 'Out-strength';
xticks(1:size(Ind,1))
xticklabels(NodiTutti)

zticks(1:size(Ind,2))
zticklabels(2003:2021)
title('Out-strength')
set(gca,'fontsize',10,'Fontweight','bold')

figure
subplot(1,2,1)
histogram(Ind,10,'FaceColor',[0 0.4470 0.7410],'EdgeAlpha',0.5,'FaceAlpha',.7,'Normalization','probability')
hold on
histogram(Outd,10,'FaceColor',[0.8500 0.3250 0.0980],'EdgeAlpha',0.5,'FaceAlpha',.7,'Normalization','probability')
axis tight
grid on
legend('In-degree','Out-degree')
title('Degree distriution')
set(gca,'fontsize',10,'Fontweight','bold')
xlabel('Degree')
ylabel('Probability')
subplot(1,2,2)
histogram(Ins,10,'FaceColor',[0 0.4470 0.7410],'EdgeAlpha',0.5,'FaceAlpha',.7,'Normalization','probability')
hold on
histogram(Outs,10,'FaceColor',[0.8500 0.3250 0.0980],'EdgeAlpha',0.5,'FaceAlpha',.7,'Normalization','probability')
axis tight
grid on
legend('In-strength','Out-strength')
title('Strength distriution')
set(gca,'fontsize',10,'Fontweight','bold')
xlabel('Strength')
ylabel('Probability')