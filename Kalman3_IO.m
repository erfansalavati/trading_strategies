%% Importing Data
%   [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A100001:C200000');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A200001:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C100000');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C40000');
% [Data,txt,~]=xlsread('tullow-data.xlsx','A40000:C57000');
% [Data,txt,~]=xlsread('tullow-data.xlsx','A1:C40000');
% [Data,txt,~]=xlsread('tullow-data.xlsx');
% [Data,txt,~]=xlsread('BTC-ETH-17-1-1-18-2-3.csv');
[Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv');
% [Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv','A1:C200000');
% [Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv','A200000:C253327');


% Data=Data(150000:200000,:);
Y=Data(:,1);
X1=Data(:,2);
X2=Data(:,3);
Index=txt(2:end,1);
T=length(Y);

%% Kalman Filter

[beta,ZScore]=Kalman_Filter3(Y,X1,X2);

% Kalman coefficients
figure(1);
subplot(3,1,1);
plot(beta(1,:));
subplot(3,1,2);
plot(beta(2,:));
subplot(3,1,3);
plot(beta(3,:));

figure(3);
plot(ZScore);
axis([0 T -std(ZScore) std(ZScore) ]);

%% Optimize Thresholds

thPer=5000;
thRelax=1;
thWin=10000;
[thY,thX,thYcl,thXcl]=Parameter_Optimizer3(Y,X1,X2,ZScore,beta(1,:),beta(2,:),thPer,thRelax,thWin);

%% Bollinger Band Trading
K=0.002; %Transaction Cost
T=length(Y);

% Position (BTC,ETH)
% thY=thYOpt*ones(T,1); %threshold to buy BTC
% thX=thXOpt*ones(T,1); %threshold to buy ETH
% thYcl=thYclOpt*ones(T,1); %threshold to close BTC position
% thXcl=thXclOpt*ones(T,1); %%threshold to close ETH position
% thY=-3*ones(T,1); %threshold to buy BTC
% thX=3*ones(T,1); %threshold to buy ETH
% thYcl=-9*ones(T,1); %threshold to close BTC position
% thXcl=9*ones(T,1); %%threshold to close ETH position
val=100;

pos=zeros(T,3);

longY=[];
shortY=[];
close=[];
PnL=zeros(T,1);
for t=2:T
    if (ZScore(t)<thY(t))&&(ZScore(t-1)>=thY(t))&&(pos(t-1,1)<=0)
        pos(t,:)=[val/Y(t) , -val*beta(1,t)./X1(t), -val*beta(2,t)./X2(t)];
        longY=[longY , t];
    elseif (ZScore(t)>thX(t))&&(ZScore(t-1)<=thX(t))&&(pos(t-1,1)>=0) 
        pos(t,:)=[-val/Y(t) , val*beta(1,t)./X1(t) , val*beta(2,t)./X2(t)];
        shortY=[shortY , t];
    elseif (ZScore(t)<thYcl(t))&&(pos(t-1,1)>0)
        pos(t,:)=[0 , 0 , 0];
        close=[close , t];
    elseif (ZScore(t)>thXcl(t))&&(pos(t-1,1)<0)
        pos(t,:)=[0 , 0 , 0];
        close=[close , t];
    else
        pos(t,:)=pos(t-1,:);
    end
%     PnL(t)=pos(t-1,1).*(Y(t)-Y(t-1)) + pos(t-1,2).*(X(t)-X(t-1))...
%         -K/2*abs(pos(t,1)-pos(t-1,1)).*Y(t-1)-K/2*abs(pos(t,2)-pos(t-1,2)).*X(t-1);
%     curPosPnL=curPosPnL+PnL(t);
end

PnL(2:end)=pos(1:end-1,1).*(Y(2:end)-Y(1:end-1)) + pos(1:end-1,2).*(X1(2:end)-X1(1:end-1))...
+ pos(1:end-1,3).*(X2(2:end)-X2(1:end-1))-K/2*abs(pos(2:end,1)-pos(1:end-1,1)).*Y(1:end-1)...
-K/2*abs(pos(2:end,2)-pos(1:end-1,2)).*X1(1:end-1) -K/2*abs(pos(2:end,3)-pos(1:end-1,3)).*X2(1:end-1);


netVal=cumsum(PnL);
BnH=Y-Y(1);
lev=1; %Leverage
margin=[0;1/lev*(abs(pos(1:end-1,1)).*Y(2:end)+abs(pos(1:end-1,2)).*X1(2:end)...
    +abs(pos(1:end-1,3)).*X2(2:end))-min(netVal(2:end),0)];
totMargin=max(margin);
APR=netVal(end)/totMargin;
[netVal(end) totMargin APR]


sum(~((pos(2:end,1)==pos(1:end-1,1))&(pos(2:end,2)==pos(1:end-1,2)))) %Number of Transactions


%% Figures

% Kalman coefficients
figure(1);
subplot(3,1,1);
plot(beta(1,:));
subplot(3,1,2);
plot(beta(2,:));
subplot(3,1,3);
plot(beta(3,:));

% 
figure(2);
subplot(2,1,1);
plot(netVal);
subplot(2,1,2);
plot(margin);

%
figure(3);
%subplot(3,1,3);
plot(ZScore);
axis([0 T -std(ZScore) std(ZScore) ]);
hold on;
plot(thY);
hold on;
plot(thX);
hold on;
plot(thYcl);
hold on;
plot(thXcl);
hold off;
 
%
figure(4);
subplot(3,1,1);
plot(Y);
hold on;
plot(longY,Y(longY),'.','markers',12);
hold on;
plot(shortY,Y(shortY),'.','markers',12);
hold on;
plot(close,Y(close),'.','markers',12);
hold off;
subplot(3,1,2);
plot(X1);
hold on;
plot(longY,X1(longY),'.','markers',12);
hold on;
plot(shortY,X1(shortY),'.','markers',12);
hold on;
plot(close,X1(close),'.','markers',12);
hold off;
subplot(3,1,3);
plot(X2);
hold on;
plot(longY,X2(longY),'.','markers',12);
hold on;
plot(shortY,X2(shortY),'.','markers',12);
hold on;
plot(close,X2(close),'.','markers',12);
hold off;


[netVal(end) totMargin APR]


