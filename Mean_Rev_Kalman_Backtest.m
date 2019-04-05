%% Importing Data
%  [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C240861');
[Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A50000:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C40000');

% Data=Data(150000:200000,:);
BTC=Data(:,1);
ETH=Data(:,2);
Index=txt(2:end,1);

%%

cost=0.005; %Transaction Cost

logBTC=log(BTC);
logETH=log(ETH);

T=length(BTC);



% Augment x with ones to accommodate possible offset in the regression
% between y vs x.
y=logBTC;
x=[logETH ones(T,1)];
delta=0.0002; % delta=0 allows no change (like traditional linear regression).
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
Ve=0.001;
% Initialize beta(:, 1) to zero
beta(:,1)=0;
ZScore=NaN(T,1);
for t=1:T
    if (t > 1)
        beta(:, t)=beta(:, t-1); % state prediction.
        R=P+Vw; % state covariance prediction.
    end
    yhat(t)=x(t, :)*beta(:, t); % measurement prediction.
    Q=x(t, :)*R*x(t, :)'+Ve; % measurement variance prediction.
    % Observe y(t)
    e(t)=y(t)-yhat(t); % measurement prediction error
    K=R*x(t, :)'/Q; % Kalman gain
    beta(:, t)=beta(:, t)+K*e(t); % State update.
    P=R-K*x(t, :)*R; % State covariance update
    ZScore(t)=(logBTC(t)-beta(1,t)'.*logETH(t)-beta(2,t)')./sqrt(Q);
end




% marketVal=abs(normdiff);

% Position (BTC,ETH)

thBTC=-0.03; %threshold to buy BTC
thETH=0.03; %threshold to buy ETH
val=100;   %market value of each position

position=zeros(T,2);
for t=2:T
    if (ZScore(t)<thBTC)&&(position(t-1,1)<=0)
%         position(t,:)=[val/BTC(t) , -val/ETH(t)];
%         position(t,:)=[val/BTC(t) , -val*abs(beta(1,t))/ETH(t)];
        position(t,:)=[val/BTC(t) , 0];
    elseif (ZScore(t)>thETH)&&(position(t-1,1)>=0) 
%         position(t,:)=[-val/BTC(t) , val/ETH(t)];
%         position(t,:)=[-val/BTC(t) , val*abs(beta(1,t))/ETH(t)];
        position(t,:)=[0 , val/ETH(t)];
%         position(t,:)=[0 , 0];
%         position(t,:)=[0 , sign(beta(1,t))*val/ETH(t)];
    else
        position(t,:)=position(t-1,:);
    end
end

% position=[-normdiff./BTC , normdiff./ETH];

marketVal=repmat(val,T,1);

PnL=position(1:end-1,1).*(BTC(2:end)-BTC(1:end-1)) + position(1:end-1,2).*(ETH(2:end)-ETH(1:end-1))...
    -cost/2*abs(position(2:end,1)-position(1:end-1,1)).*BTC(2:end)-cost/2*abs(position(2:end,2)-position(1:end-1,2)).*ETH(2:end);
PnL=[0;PnL];




% netVal=cumsum(PnL)+position(:,1).*BTC+position(:,2).*ETH;
netVal=cumsum(PnL);
BnH=BTC-BTC(1);

netVal(end)
sum(~((position(2:end,1)==position(1:end-1,1))&(position(2:end,2)==position(1:end-1,2))))

subplot(5,1,1);
plot(netVal);
subplot(5,1,2);
plot(BTC);
subplot(5,1,3);
plot(ETH);
subplot(5,1,4);
plot(beta(1,:));
axis([0 T -.2 .5]);
subplot(5,1,5);
axis([0 T -.1 .1]);
 
% figure;
% plot(movstd(ZScore(1000:end),[2000 0]));

%% Outputs
% delete('BTC-ETH-Output.xlsx');
outputFileName=['Output',datestr(now,'HHMMSS'),'.xlsx'];
xlswrite(outputFileName,{'Index','BTC','ETH','logBTC','logETH','difflog',...
    'Moving Average','normdiff','marketVal','positionBTC','positionETH','PnL','netVal','BnH'},1,'A1');
xlswrite(outputFileName,[BTC,ETH,logBTC,logETH,difflog,...
    movAverage,ZScore,marketVal,position,PnL,netVal,BnH],1,'B2');
xlswrite(outputFileName,Index,1,'A2');

%% Optimizimg the Mean Reverting Strategy

% Y=fmincon(@(Y) (net_Value(BTC,ETH,1000,Y(1),Y(2))),[-0.5,0.5]);

thBTC=-2:.2:2;
thETH=-2:.2:2;
[thBTC,thETH]=meshgrid(thBTC,thETH);

V=arrayfun(@(x,y) (net_Value(BTC,ETH,1000,x,y)),thBTC,thETH);

figure;
surf(thBTC,thETH,V);

[~,I]=max(V(:));
[thBTC(I),thETH(I)]

% -0.6 , 0

