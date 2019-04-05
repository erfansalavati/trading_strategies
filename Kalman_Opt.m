%% Importing Data
%   [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A100001:C200000');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A200001:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C100000');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C40000');
% [Data,txt,~]=xlsread('tullow-data.xlsx','A40000:C57000');
% [Data,txt,~]=xlsread('tullow-data.xlsx','A1:C40000');
 [Data,txt,~]=xlsread('tullow-data.xlsx');
% [Data,txt,~]=xlsread('BTC-ETH-17-1-1-18-2-3.csv');
%[Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv');
% [Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv','A1:C200000');
% [Data,txt,~]=xlsread('BTC-ETH-BCH-17-8-1-18-1-24.csv','A200000:C253327');


% Data=Data(150000:200000,:);
Y=Data(:,1);
X=Data(:,2);
Index=txt(2:end,1);

%%
K=0.002; %Transaction Cost
% beta1=1;
logY=log(Y);
logX=log(X);
% ratio=BTC./ETH;
%window=2040;
% winReg=100;
T=length(Y);

% movAverage=movmean(difflog,[winMA 0]);
% movStD=movstd(difflog,[winMA 0]);
% ZScore=(difflog-movAverage)./movStD;
% ZScore=[zeros(winMA,1);ZScore(winMA+1:end)]; 

% marketVal=abs(normdiff);

% Position (BTC,ETH)


thY=-1*ones(T,1); %threshold to buy BTC
thX=1*ones(T,1); %threshold to buy ETH
thYcl=-5*ones(T,1); %threshold to close BTC position
thXcl=5*ones(T,1); %%threshold to close ETH position
val=100;   %market value of each position


pos=zeros(T,2);

longY=[];
shortY=[];
close=[];
PnL=zeros(T,1);

thPer=5000;
thRelax=1;
thWin=10000;

%initialize Optimization
A=[-1 0 1 0; 1 -1 0 0; 0 1 0 -1];
b=[0;0;0];
LB=[-10;-10;-10;-10];
UB=[10;10;10;10];
options = optimoptions(@ga,'Display','iter');

% %Initialize Kalman Filter
% % Augment x with ones to accommodate possible offset in the regression
% % between y vs x.
% y=logY;
% x=[logX ones(T,1)];
% yhat=NaN(T,1); % measurement prediction
% e=NaN(T,1); % measurement prediction error
% Q=NaN(T,1); % measurement prediction error variance
% %R=NaN(T,2,2); %State covariance prediction
% % For clarity, we denote R(t|t) by P(t).
% % initialize P and beta.
% P=zeros(2);
% R=zeros(2);
% beta=NaN(2, T);
% delta=0.000001; % delta=0 allows no change (like traditional linear regression).
% Ve=1;
% Vw=delta/(1-delta)*diag(ones(2, 1));
% % Initialize beta(:, 1) to zero
% beta(:,1)=[1;0];
%     yhat(1)=x(1, :)*beta(:, 1); % measurement prediction.
%     Q=x(1, :)*R*x(1, :)'+Ve; % measurementvariance prediction.
%     % Observe y(t)
%     e(1)=y(1)-yhat(1); % measurement prediction error
%     KG=R*x(1, :)'/Q; % Kalman gain
%     beta(:, 1)=beta(:, 1)+KG*e(1); % State update.
%     P=R-KG*x(1, :)*R; % State covariance update

for t=2:T
%       Optimize thresholds
        if (mod(t,thPer)==0)
            t
            fitnessfcn=@(Z)(-Ret(Y(max(1,t-thWin+1):t),X(max(1,t-thWin+1):t),ZScore(max(1,t-thWin+1):t),Z(1),Z(2),Z(3),Z(4),beta(1,max(1,t-thWin+1):t)));
            I=ga(fitnessfcn,4,A,b,[],[],LB,UB,[],[])%,options)
            thY(t:min(t+thPer-1,T))=I(1);
            thX(t:min(t+thPer-1,T))=I(2);
            thYcl(t:min(t+thPer-1,T))=I(3);
            thXcl(t:min(t+thPer-1,T))=I(4);
            LB=[I(1)-abs(I(1))*thRelax;I(2)-abs(I(2))*thRelax;I(3)-abs(I(3))*thRelax;I(4)-abs(I(4))*thRelax];
            UB=[I(1)+abs(I(1))*thRelax;I(2)+abs(I(2))*thRelax;I(3)+abs(I(3))*thRelax;I(4)+abs(I(4))*thRelax];
            LB=max(LB,[-10;-10;-10;-10]);
            UB=min(UB,[10;10;10;10]);
%             thY(t:min(t+thPer-1,T))=thRelax*I(2)+(1-thRelax)*thY(t-1);
%             thX(t:min(t+thPer-1,T))=thRelax*I(3)+(1-thRelax)*thX(t-1);
%             thYcl(t:min(t+thPer-1,T))=thRelax*I(4)+(1-thRelax)*thYcl(t-1);
%             thXcl(t:min(t+thPer-1,T))=thRelax*I(5)+(1-thRelax)*thXcl(t-1);
        end;
        
        
% %        Kalman Filter
%         beta(:, t)=beta(:, t-1); % state prediction.
%         R=P+Vw; % state covariance prediction.
%         yhat(t)=x(t, :)*beta(:, t); % measurement prediction.
%         Q=x(t, :)*R*x(t, :)'+Ve; % measurementvariance prediction.
%         % Observe y(t)
%         e(t)=y(t)-yhat(t); % measurement prediction error
%         KG=R*x(t, :)'/Q; % Kalman gain
%         beta(:, t)=beta(:, t)+KG*e(t); % State update.
%         P=R-KG*x(t, :)*R; % State covariance update


%Trade
        
    if (ZScore(t)<thY(t))&&(ZScore(t-1)>=thY(t))&&(pos(t-1,1)<=0)
        pos(t,:)=[val/Y(t) , -val/X(t)];
%         pos(t,:)=[val/Y(t) , -val*abs(beta(1,t))/X(t)];
%         pos(t,:)=[val/Y(t) , -val/Y(t)*beta1(t)];
%         position(t,:)=[val/Y(t) , 0];
%         curPosPnL=0;
        longY=[longY , t];
    elseif (ZScore(t)>thX(t))&&(ZScore(t-1)<=thX(t))&&(pos(t-1,1)>=0) 
        pos(t,:)=[-val/Y(t) , val/X(t)];
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

noTrad=2*thPer;
pos(1:noTrad,:)=0;
portf(1:noTrad,:)=0;

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
figure(1);
subplot(2,1,1);
plot(netVal);
subplot(2,1,2);
plot(ZScore);
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
figure(2);
subplot(2,1,1);
plot(netVal);
subplot(2,1,2);
plot(difflog);
hold on;
plot(movAverage+thY.*movStD);
hold on;
plot(movAverage+thX.*movStD);
hold on;
plot(movAverage+thYcl.*movStD);
hold on;
plot(movAverage+thXcl.*movStD);
% axis([1 T 1 2]);
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



%% Outputs
% delete('BTC-ETH-Output.xlsx');
outputFileName=['Output',datestr(now,'HHMMSS'),'.xlsx'];
xlswrite(outputFileName,{'Index','BTC','ETH','logBTC','logETH','difflog',...
    'Moving Average','normdiff','marketVal','positionBTC','positionETH','PnL','netVal','BnH'},1,'A1');
xlswrite(outputFileName,[Y,X,logY,logX,difflog,...
    movAverage,ZScore,marketVal,pos,PnL,netVal,BnH],1,'B2');
xlswrite(outputFileName,Index,1,'A2');
