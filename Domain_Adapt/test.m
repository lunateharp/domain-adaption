clc; clear; close
rng(2017)
K = 3;
mu = cell(K,1); SIGMA = cell(K,1);
N = zeros(K,1);
X = cell(K,1);
for k = 1:K
    mu{k} = 5*[rand;rand];
    N(k) = randi([10,20]);
    SIGMA{k} = randn(2);
    SIGMA{k} = SIGMA{k}*SIGMA{k}';
    X{k} = mvnrnd(mu{k},SIGMA{k},N(k));
    plot(X{k}(:,1),X{k}(:,2),'.')
    hold on
end

M = 7;
t = [randi(6)-3,randi(6)-3]; % translation
Y = zeros(M,2); 
L = zeros(M,1); % true label
for j = 1:M
    L(j) = randi(K);
    Y(j,:) = t+mvnrnd(mu{L(j)},SIGMA{L(j)})+ .1*[randn,randn];
end

plot(Y(:,1),Y(:,2),'ko')
P = Knowtran_D(X,Y,0);

