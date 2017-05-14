clc; clear; close
set(groot,'DefaultLineLineWidth',2,'DefaultLineMarkerSize',10)
rng(2017)
K = 2;
mu1 = [0,1]; mu2 = [2,3];
SIGMA1 = eye(2); SIGMA2 = .5*eye(2);
N1 = 10; N2 = 8;
X1 = mvnrnd(mu1,SIGMA1,N1);
X2 = mvnrnd(mu2,SIGMA2,N2);
M = 10;
t = [randi(6)-3,randi(6)-3]; % translation
Y = zeros(M,2); 
L = zeros(M,1); % true label
for j = 1:10
    L(j) = randi(2);
    if L(j) == 1
        Y(j,:) = t+mvnrnd(mu1,SIGMA1)+ .1*[randn,randn];
    else
        Y(j,:) = t+mvnrnd(mu2,SIGMA2)+ .1*[randn,randn];
    end
end
plot(X1(:,1),X1(:,2),'r.')
hold on
plot(X2(:,1),X2(:,2),'b.')
plot(Y(:,1),Y(:,2),'ko')

% Initial
P0 = .5*ones(M,1); 
w1 = P0./sum(P0);
w2 = (1-P0)./sum(1-P0);
% solve for W1
C1 = 
% min N1*W1^2 + N2*W2^2