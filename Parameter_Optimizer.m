%Threshold_Optimizer(Y,X,ZScore,beta)
function [thY,thX,thYcl,thXcl]=Parameter_Optimizer(Y,X,ZScore,beta1,thPer,thRelax,thWin)
K=0.002; %Transaction Cost
T=length(Y);

% switch nargin
%     case 3
%         beta1=ones(T,1);
%     case 7
%         
% end

logY=log(Y);
logX=log(X);



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

% thPer=5000;
% thRelax=1;
% thWin=10000;

%initialize Optimization
A=[-1 0 1 0; 1 -1 0 0; 0 1 0 -1];
b=[0;0;0];
LB=[-10;-10;-10;-10];
UB=[10;10;10;10];
options = optimoptions(@ga,'Display','iter');


for t=2:T
    %Optimize thresholds
    if (mod(t,thPer)==0)
        t
        fitnessfcn=@(Z)(-Bollinger(Y(max(1,t-thWin+1):t),X(max(1,t-thWin+1):t),ZScore(max(1,t-thWin+1):t),Z(1),Z(2),Z(3),Z(4),beta1(max(1,t-thWin+1):t)));
        [I,fval]=ga(fitnessfcn,4,A,b,[],[],LB,UB,[],[])%,options)
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
    
%     %Trade
%     if (ZScore(t)<thY(t))&&(ZScore(t-1)>=thY(t))&&(pos(t-1,1)<=0)
%         pos(t,:)=[val/Y(t) , -val/X(t)*beta1(t)];
%         %         pos(t,:)=[val/Y(t) , -val*abs(beta(1,t))/X(t)];
%         %         pos(t,:)=[val/Y(t) , -val/Y(t)*beta1(t)];
%         %         position(t,:)=[val/Y(t) , 0];
%         %         curPosPnL=0;
%         longY=[longY , t];
%     elseif (ZScore(t)>thX(t))&&(ZScore(t-1)<=thX(t))&&(pos(t-1,1)>=0)
%         pos(t,:)=[-val/Y(t) , val/X(t)*beta1(t)];
%         %         pos(t,:)=[-val/Y(t) , val/Y(t)*beta1(t)];
%         %         position(t,:)=[0 , val/X(t)];
%         %         position(t,:)=[0 , 0];
%         %         position(t,:)=[0 , sign(beta(1,t))*val/ETH(t)];
%         %         curPosPnL=0;
%         shortY=[shortY , t];
%     elseif (ZScore(t)<thYcl(t))&&(pos(t-1,1)>0)
%         pos(t,:)=[0 , 0];
%         %         curPosPnL=0;
%         close=[close , t];
%     elseif (ZScore(t)>thXcl(t))&&(pos(t-1,1)<0)
%         pos(t,:)=[0 , 0];
%         %         curPosPnL=0;
%         close=[close , t];
%     else
%         pos(t,:)=pos(t-1,:);
%     end
%     %     PnL(t)=pos(t-1,1).*(Y(t)-Y(t-1)) + pos(t-1,2).*(X(t)-X(t-1))...
%     %         -K/2*abs(pos(t,1)-pos(t-1,1)).*Y(t-1)-K/2*abs(pos(t,2)-pos(t-1,2)).*X(t-1);
%     %     curPosPnL=curPosPnL+PnL(t);
end

% noTrad=2*thPer;
% pos(1:noTrad,:)=0;
% portf(1:noTrad,:)=0;
% 
% PnL(2:end)=pos(1:end-1,1).*(Y(2:end)-Y(1:end-1)) + pos(1:end-1,2).*(X(2:end)-X(1:end-1))...
% -K/2*abs(pos(2:end,1)-pos(1:end-1,1)).*Y(1:end-1)-K/2*abs(pos(2:end,2)-pos(1:end-1,2)).*X(1:end-1);
% 
% 
% netVal=cumsum(PnL);
% BnH=Y-Y(1);
% lev=1; %Leverage
% margin=1/lev*(abs(pos(:,1)).*Y+abs(pos(:,2)).*X)-min(netVal,0);
% totMargin=max(margin);
% APR=netVal(end)/totMargin;
% [netVal(end) totMargin APR]


% sum(~((pos(2:end,1)==pos(1:end-1,1))&(pos(2:end,2)==pos(1:end-1,2)))) %Number of Transactions

