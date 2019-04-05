%% simulated prices

sT=1000;

delta=0.0001;
% Vw=delta/(1-delta)*diag(ones(2, 1));
Vw1=0.0001;
Vw2=0.000001;
Ve=0.001;

sbeta=NaN(2, sT);

sbeta(1,:)=cumsum(sqrt(Vw1)*randn(1,sT));
sbeta(2,:)=cumsum(sqrt(Vw2)*randn(1,sT));

sdiffX=randn(sT,1);

sdiffY= sbeta(1,:)'.*sdiffX + sbeta(2,:)'+sqrt(Ve)*randn(sT,1);

sX=cumsum(sdiffX);
sY=cumsum(sdiffY);

figure(3);

subplot(4,1,1);
plot(sbeta(1,:));
subplot(4,1,2);
plot(sbeta(2,:));
subplot(4,1,3);
plot(sX);
subplot(4,1,4);
plot(sY);


%% Kalman Filtering on Simulated data



% Augment x with ones to accommodate possible offset in the regression
% between y vs x.
y=sdiffY;
x=[sdiffX ones(sT,1)];
yhat=NaN(sT,1); % measurement prediction
e=NaN(sT,1); % measurement prediction error
Q=NaN(sT,1); % measurement prediction error variance
%R=NaN(T,2,2); %State covariance prediction
% For clarity, we denote R(t|t) by P(t).
% initialize P and beta.
P=zeros(2);
R=zeros(2);
beta=NaN(2, sT);
Vw=diag([Vw1 Vw2]);
% Ve=0.001;
% Initialize beta(:, 1) to zero
beta(:,1)=0;
ZScore=NaN(sT,1);
for t=1:sT
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
%     noise(t)=(Y(t)-beta(1,t)*X(t)-beta(2,t))/sqrt(Q);
    ZScore(t)=(diffY(t)-beta(1,t)*diffX(t)-beta(2,t))/sqrt(Q);
end

figure(4);
subplot(2,1,1);
plot(sbeta(1,:));
hold on;
plot(beta(1,:));
hold off;
subplot(2,1,2);
plot(sbeta(2,:));
hold on;
plot(beta(2,:));
hold off;


%% Importing Data
%  [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A50000:C240861');
% [Data,txt,~]=xlsread('BTC-ETH-Data.xlsx','A1:C40000');
[Data,txt,~]=xlsread('tullow-data.xlsx');
% [Data,txt,~]=xlsread('tullow-data.xlsx','A1:C1000');

% Data=Data(150000:200000,:);
Y=Data(:,2);
X=Data(:,1);
Index=txt(2:end,1);

T=length(Y);

%%

cost=0.005; %Transaction Cost

% logY=log(Y);
% logX=log(X);

T=length(Y);


diffY=[0;Y(2:end)-Y(1:end-1)];

diffX=[0;X(2:end)-X(1:end-1)];



% Augment x with ones to accommodate possible offset in the regression
% between y vs x.
y=diffY;
x=[diffX ones(T,1)];
delta=0.5;    %0.0001 % delta=0 allows no change (like traditional linear regression).
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
%     noise(t)=(Y(t)-beta(1,t)*X(t)-beta(2,t))/sqrt(Q);
    ZScore(t)=(diffY(t)-beta(1,t)*diffX(t)-beta(2,t))/sqrt(Q);
end

% winMA=1000;
% ZScore=movmean(noise,[winMA 0]);



% marketVal=abs(normdiff);

% Position (BTC,ETH)

thY=-0.003*ones(T,1); %threshold to buy BTC
thX=0.003*ones(T,1); %threshold to buy ETH
val=1;   %market value of each position

pos=zeros(T,2);
thPer=4000;
thWin=4000;
thRelax=0.5;
longY=[];
shortY=[];
close=[];
notrade=2;
for t=notrade:T
    %Optimize Thresholds
%         if (mod(t,thPer)==0)
%             t
%             thYarray=-0.05:0.005:0.05;
%             thXarray=-0.05:0.005:0.05;
%             [thYarray,thXarray]=meshgrid(thYarray,thXarray);
%             V=arrayfun(@(x,y) (net_Value_Kal(ZScore(t-thWin+1:t),Y(t-thWin+1:t),X(t-thWin+1:t)...
%                 ,x,y)),thYarray,thXarray);
%             [~,I]=max(V(:));
%             thY(t:min(t+thPer-1,T))=thRelax*thYarray(I)+(1-thRelax)*thY(t-1);
%             thX(t:min(t+thPer-1,T))=thRelax*thXarray(I)+(1-thRelax)*thX(t-1);
%         end;

    if (ZScore(t)<thY(t))&&(pos(t-1,1)<=0)
%         pos(t,:)=[val/Y(t) , -val/X(t)];
%         pos(t,:)=[val/Y(t) , -val*abs(beta(1,t))/X(t)];
        pos(t,:)=[val/Y(t) , -val/Y(t)*beta(1,t)];
%         position(t,:)=[val/Y(t) , 0];
        longY=[longY , t];
    elseif (ZScore(t)>thX(t))&&(pos(t-1,1)>=0) 
%         pos(t,:)=[-val/Y(t) , val/X(t)];
        pos(t,:)=[-val/Y(t) , val/Y(t)*beta(1,t)];
%         position(t,:)=[0 , val/X(t)];
%         position(t,:)=[0 , 0];
%         position(t,:)=[0 , sign(beta(1,t))*val/ETH(t)];
        shortY=[shortY , t];
    else
        pos(t,:)=pos(t-1,:);
    end
end

% position=[-normdiff./BTC , normdiff./ETH];

marketVal=repmat(val,T,1);

PnL=pos(1:end-1,1).*(Y(2:end)-Y(1:end-1)) + pos(1:end-1,2).*(X(2:end)-X(1:end-1))...
    -cost/2*abs(pos(2:end,1)-pos(1:end-1,1)).*Y(2:end)-cost/2*abs(pos(2:end,2)-pos(1:end-1,2)).*X(2:end);
PnL=[0;PnL];




% netVal=cumsum(PnL)+pos(:,1).*Y+pos(:,2).*X;
netVal=cumsum(PnL);
BnH=Y-Y(1);
lev=10; %Leverage
margin=1/lev*(max(pos(:,1),0).*Y+max(pos(:,2),0).*X)-min(netVal,0);

[netVal(end) , max(margin)]
sum(~((pos(2:end,1)==pos(1:end-1,1))&(pos(2:end,2)==pos(1:end-1,2))))

%  close all;
figure(1);
subplot(3,1,1);
plot(netVal);
subplot(3,1,2);
plot(Y);
hold on;
plot(longY,Y(longY),'.');
hold off;
subplot(3,1,3);
plot(X);
hold on;
plot(shortY,X(shortY),'.');
hold off;
% subplot(5,1,4);
% plot(beta(1,:));
% axis([0 T -.2 .5]);
% subplot(4,1,4);
figure(2);
subplot(3,1,1);
plot(ZScore);
% axis([0 T -.001 .001]);
hold on;
plot(thY);
hold on;
plot(thX);
hold off;
subplot(3,1,2);
plot(beta(1,:));
subplot(3,1,3);
plot(beta(2,:));
% figure;
% plot(movstd(ZScore(1000:end),[2000 0]));
% 
% figure;
% subplot(2,1,1);
% plot(mm);
% axis([0 T -0.02 0.02])
% subplot(2,1,2);
% plot(log(Y./X))


%% Outputs
% delete('BTC-ETH-Output.xlsx');
outputFileName=['Output',datestr(now,'HHMMSS'),'.xlsx'];
xlswrite(outputFileName,{'Index','BTC','ETH','logBTC','logETH','difflog',...
    'Moving Average','normdiff','marketVal','positionBTC','positionETH','PnL','netVal','BnH'},1,'A1');
xlswrite(outputFileName,[Y,X,logY,logX,difflog,...
    movAverage,ZScore,marketVal,pos,PnL,netVal,BnH],1,'B2');
xlswrite(outputFileName,Index,1,'A2');

%% Optimizimg the Mean Reverting Strategy

% Y=fmincon(@(Y) (net_Value(BTC,ETH,1000,Y(1),Y(2))),[-0.5,0.5]);

thY=-2:.2:2;
thX=-2:.2:2;
[thY,thX]=meshgrid(thY,thX);

V=arrayfun(@(x,y) (net_Value(Y,X,1000,x,y)),thY,thX);

figure;
surf(thY,thX,V);

[~,I]=max(V(:));
[thY(I),thX(I)]

% -0.6 , 0

