function V=net_Value_Kal(ZScore,Y,X,thY,thX)
cost=0.005; %Transaction Cost


T=length(Y);



val=100;   %market value of each position

pos=zeros(T,2);

for t=2:T
    if (ZScore(t)<thY)&&(pos(t-1,1)<=0)
%         position(t,:)=[val/Y(t) , -val/X(t)];
%         position(t,:)=[val/Y(t) , -val*abs(beta(1,t))/X(t)];
        pos(t,:)=[val/Y(t) , -val/Y(t)*beta(1,t)];
%         position(t,:)=[val/Y(t) , 0];
    elseif (ZScore(t)>thX)&&(pos(t-1,1)>=0) 
%         position(t,:)=[-val/Y(t) , val/X(t)];
%         position(t,:)=[-val/Y(t) , val*abs(beta(1,t))/X(t)];
        pos(t,:)=[-val/Y(t) , val/Y(t)*beta(1,t)];
%         position(t,:)=[0 , val/X(t)];
    else
        pos(t,:)=pos(t-1,:);
    end
end


PnL=pos(1:end-1,1).*(Y(2:end)-Y(1:end-1)) + pos(1:end-1,2).*(X(2:end)-X(1:end-1))...
    -cost/2*abs(pos(2:end,1)-pos(1:end-1,1)).*Y(2:end)-cost/2*abs(pos(2:end,2)-pos(1:end-1,2)).*X(2:end);
PnL=[0;PnL];


V=sum(PnL);

end