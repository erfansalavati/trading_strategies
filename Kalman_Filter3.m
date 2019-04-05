function [beta,ZScore]=Kalman_Filter3(Y,X1,X2)

logY=log(Y);
logX1=log(X1);
logX2=log(X2);

T=length(Y);




% Augment x with ones to accommodate possible offset in the regression
% between y vs x.
y=logY;
x=[logX1 logX2 ones(T,1)];
delta=0.000001; % delta=0 allows no change (like traditional linear regression).
yhat=NaN(T,1); % measurement prediction
e=NaN(T,1); % measurement prediction error
Q=NaN(T,1); % measurement prediction error variance
%R=NaN(T,2,2); %State covariance prediction
% For clarity, we denote R(t|t) by P(t).
% initialize P and beta.
P=zeros(3);
R=zeros(3);
beta=NaN(3, T);
Vw=delta/(1-delta)*diag(ones(3, 1));
Ve=100;
% Initialize beta(:, 1) to zero
beta(:,1)=[1;0;0];
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


ZScore=(logY-beta(1,:)'.*logX1-beta(2,:)'.*logX2-beta(3,:)')./sqrt(Q);

