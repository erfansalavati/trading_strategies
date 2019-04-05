%% Importing Data
%   [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A100001:C200000');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A200001:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C100000');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C40000');
% [Data,txt,~]=xlsread('tullow-data.xlsx','A40000:C57000');
% [Data,txt,~]=xlsread('tullow-data.xlsx','A1:C40000');
% [Data,txt,~]=xlsread('tullow-data.xlsx');
[Data,txt,~]=xlsread('BTC-ETH-17-1-1-18-2-3.csv');
%  [Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv');
% [Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv','A1:C200000');
% [Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv','A200000:C253327');


% Data=Data(150000:200000,:);
Y=Data(:,1);
X=Data(:,2);
Index=txt(2:end,1);


%% Kalman Filter


logY=log(Y);
logX=log(X);

T=length(Y);



% Augment x with ones to accommodate possible offset in the regression
% between y vs x.
y=logY;
x=[logX ones(T,1)];
delta=0.000001; % delta=0 allows no change (like traditional linear regression).
yhat=NaN(T,1); % measurement prediction
e=NaN(T,1); % measurement prediction error
Q=NaN(T,1); % measurement prediction error variance
%R=NaN(T,2,2); %State covariance prediction
% For clarity, we denote R(t|t) by P(t).
% initialize P and beta.
P=zeros(2);
R=zeros(2);
beta=NaN(2, T);
Vw=delta/(1-delta)*diag(ones(2, 1));
Ve=1;
% Initialize beta(:, 1) to zero
beta(:,1)=[1;0];
for t=1:T
    if (t > 1)
        beta(:, t)=beta(:, t-1); % state prediction.
        R=P+Vw; % state covariance prediction.
    end
    yhat(t)=x(t, :)*beta(:, t); % measurement prediction.
    Q=x(t, :)*R*x(t, :)'+Ve; % measurementvariance prediction.
    % Observe y(t)
    e(t)=y(t)-yhat(t); % measurement prediction error
    KG=R*x(t, :)'/Q; % Kalman gain
    beta(:, t)=beta(:, t)+KG*e(t); % State update.
    P=R-KG*x(t, :)*R; % State covariance update
end


ZScore=(logY-beta(1,:)'.*logX-beta(2,:)')./sqrt(Q);



%% Trade according to Signal from ZScore
K=0.002; %Transaction Cost

% Position (BTC,ETH)
thY=thYOpt*ones(T,1); %threshold to buy BTC
thX=thXOpt*ones(T,1); %threshold to buy ETH
thYcl=thYclOpt*ones(T,1); %threshold to close BTC position
thXcl=thXclOpt*ones(T,1); %%threshold to close ETH position
% thY=-3*ones(T,1); %threshold to buy BTC
% thX=3*ones(T,1); %threshold to buy ETH
% thYcl=-9*ones(T,1); %threshold to close BTC position
% thXcl=9*ones(T,1); %%threshold to close ETH position
val=100;

pos=zeros(T,2);

longY=[];
shortY=[];
close=[];
PnL=zeros(T,1);
for t=2:T
    if (ZScore(t)<thY(t))&&(ZScore(t-1)>=thY(t))&&(pos(t-1,1)<=0)
        pos(t,:)=[val/Y(t) , -val*beta(1,t)./X(t)];
%         pos(t,:)=[val/Y(t) , -val*abs(beta(1,t))/X(t)];
%         pos(t,:)=[val/Y(t) , -val/Y(t)*beta1(t)];
%         position(t,:)=[val/Y(t) , 0];
%         curPosPnL=0;
        longY=[longY , t];
    elseif (ZScore(t)>thX(t))&&(ZScore(t-1)<=thX(t))&&(pos(t-1,1)>=0) 
        pos(t,:)=[-val/Y(t) , val*beta(1,t)./X(t)];
%         pos(t,:)=[-val/Y(t) , val/Y(t)*beta1(t)];
%         position(t,:)=[0 , val/X(t)];
%         position(t,:)=[0 , 0];
%         position(t,:)=[0 , sign(beta(1,t))*val/ETH(t)];
%         curPosPnL=0;
        shortY=[shortY , t];
    elseif (ZScore(t)<thYcl(t))&&(pos(t-1,1)>0)
        pos(t,:)=[0 , 0];
%         curPosPnL=0;
        close=[close , t];
    elseif (ZScore(t)>thXcl(t))&&(pos(t-1,1)<0)
        pos(t,:)=[0 , 0];
%         curPosPnL=0;
        close=[close , t];
    else        
        pos(t,:)=pos(t-1,:);
    end
%     PnL(t)=pos(t-1,1).*(Y(t)-Y(t-1)) + pos(t-1,2).*(X(t)-X(t-1))...
%         -K/2*abs(pos(t,1)-pos(t-1,1)).*Y(t-1)-K/2*abs(pos(t,2)-pos(t-1,2)).*X(t-1);
%     curPosPnL=curPosPnL+PnL(t);
end

PnL(2:end)=pos(1:end-1,1).*(Y(2:end)-Y(1:end-1)) + pos(1:end-1,2).*(X(2:end)-X(1:end-1))...
-K/2*abs(pos(2:end,1)-pos(1:end-1,1)).*Y(1:end-1)-K/2*abs(pos(2:end,2)-pos(1:end-1,2)).*X(1:end-1);


netVal=cumsum(PnL);
BnH=Y-Y(1);
lev=1; %Leverage
margin=1/lev*(abs(pos(:,1)).*Y+abs(pos(:,2)).*X)-min(netVal,0);
totMargin=max(margin);
APR=netVal(end)/totMargin;
[netVal(end) totMargin APR]


sum(~((pos(2:end,1)==pos(1:end-1,1))&(pos(2:end,2)==pos(1:end-1,2)))) %Number of Transactions


%% Figures

%
figure(1);
subplot(2,1,1);
plot(beta(1,:));
subplot(2,1,2);
plot(beta(2,:));
%
figure(4);
subplot(2,1,1);
% plot((beta(1,:)-beta(2,:))/4+difflog');
plot((beta(1,:)-beta(2,:)));
subplot(2,1,2);
plot(difflog)
%
figure(2);
subplot(2,1,1);
plot(netVal);
subplot(2,1,2);
plot(ZScore);
axis([0 T -1 1 ]);
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
figure(3);
subplot(2,1,1);
plot(Y);
hold on;
plot(longY,Y(longY),'.','markers',12);
hold on;
plot(shortY,Y(shortY),'.','markers',12);
hold on;
plot(close,Y(close),'.','markers',12);
hold off;
subplot(2,1,2);
plot(X);
hold on;
plot(longY,X(longY),'.','markers',12);
hold on;
plot(shortY,X(shortY),'.','markers',12);
hold on;
plot(close,X(close),'.','markers',12);
hold off;

[netVal(end) totMargin APR]


%% Optimizimg the Parameters


%Mean_Rev_Opt(Y,X,winMA,thY,thX,thYcl,thXcl)
fitnessfcn=@(Z)(-Ret(Y,X,ZScore,Z(1),Z(2),Z(3),Z(4),beta(1,:)));

A=[-1 0 1 0; 1 -1 0 0; 0 1 0 -1];
b=[0;0;0];
LB=[-10;-10;-10;-10];
UB=[10;10;10;10];
options = optimoptions(@ga,'Display','iter');
I=ga(fitnessfcn,4,A,b,[],[],LB,UB,[],[],options);

 
thYOpt=I(1);
thXOpt=I(2);
thYclOpt=I(3);
thXclOpt=I(4);
[APR,sortinoRatio,finVal,totMargin,sharpeRatio]=Ret(Y,X,ZScore,thYOpt,thXOpt,thYclOpt,thXclOpt,beta(1,:));



%% Outputs
% delete('BTC-ETH-Output.xlsx');
outputFileName=['Output',datestr(now,'HHMMSS'),'.xlsx'];
xlswrite(outputFileName,{'Index','BTC','ETH','logBTC','logETH','difflog',...
    'Moving Average','normdiff','marketVal','positionBTC','positionETH','PnL','netVal','BnH'},1,'A1');
xlswrite(outputFileName,[Y,X,logY,logX,difflog,...
    movAverage,ZScore,marketVal,position,PnL,netVal,BnH],1,'B2');
xlswrite(outputFileName,Index,1,'A2');
