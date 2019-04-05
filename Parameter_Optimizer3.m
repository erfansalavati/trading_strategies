%Threshold_Optimizer(Y,X,ZScore,beta)
function [thY,thX,thYcl,thXcl]=Parameter_Optimizer3(Y,X1,X2,ZScore,beta1,beta2,thPer,thRelax,thWin)
K=0.002; %Transaction Cost
T=length(Y);
 
% switch nargin
%     case 3
%         beta1=ones(T,1);
%     case 7
%         
% end

logY=log(Y);
logX=log(X1);



thY=-1*ones(T,1); %threshold to long Y
thX=1*ones(T,1); %threshold to short Y
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
        fitnessfcn=@(Z)(-Bollinger3(Y(max(1,t-thWin+1):t),X1(max(1,t-thWin+1):t)...
            ,X2(max(1,t-thWin+1):t),ZScore(max(1,t-thWin+1):t),Z(1),Z(2),Z(3),Z(4)...
            ,beta1(max(1,t-thWin+1):t),beta2(max(1,t-thWin+1):t)));
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
end