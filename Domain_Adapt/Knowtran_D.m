function P = Knowtran_D(X,Y,alpha,tol,P0)
% Domain adaptation via optimal transport
% by solving min_P sum_k alpha_k*W_k^2
%               s.t P is probability,
%                   W_k the Wasserstein distance

% Input:
%   X - cell arrays of data sets with label 1:K
%   Y - target data sets to be labeled
%   alpha - parameters controlling the weights of the kth Wasserstein
%           distance with respect to the kth sample size
%   P0 - initial distribution if good guess is available

if nargin < 2 || nargin > 5
    error('Knowtran(X,Z,alpha,tol,P0)')
end

K = length(X);
M = size(Y,1); % the number of points in Y
N = zeros(K,1);
for k = 1:K  % the number of points in X_k
    N(k) = size(X{k},1);
end
if nargin == 2
    alpha = 1;
    tol = 1e-10;
    P0 = repmat(ones(1,K)/K,M,1);
elseif nargin == 3
    tol = 1e-10;
    P0 = repmat(ones(1,K)/K,M,1);
elseif nargin == 4
    P0 = repmat(ones(1,K)/K,M,1);
end

%% compute initial point
W0 = zeros(K,1); % the k Wasserstein distances
w0 = P0./repmat(sum(P0),M,1); % weights in each group
A = cell(K,1); % constraint matrix Ax=b
b = cell(K,1);
C = cell(K,1); % distance matrix C_ij
y0 = cell(K,1); % full dual variable
v0 = zeros(M,K); % dual variable corresponding to the jth constraint
x0 = cell(K,1); % primal variable
for k = 1:K % compute the initial OT
    C{k} = zeros(N(k),M);
    for i = 1:N(k)
        for j = 1:M
            C{k}(i,j) = (X{k}(i,:)-Y(j,:))*(X{k}(i,:)-Y(j,:))';
        end
    end
    A{k} = [kron(ones(1,M),eye(N(k)));
        kron(eye(M),ones(1,N(k)))];
    b{k} = [ones(N(k),1)/N(k);w0(:,k)];
    clear model
    model.A = sparse(A{k});
    model.obj = C{k}(:);
    model.sense = '=';
    model.rhs = b{k};
    result = gurobi(model);
    x0{k} = result.x;
    W0(k) = result.objval;
    y0{k} = result.pi;
    v0(:,k) = y0{k}(N(k)+1:end);
end
F0 = .5*(N.^alpha)'*(W0.^2);
iter = 0; maxiter = 100;
%% compute gradient of P
gradF = zeros(M,K);
    for k = 1:K
        gradF(:,k) = -(N(k)^alpha)*W0(k)*...
            (repmat(v0(:,k)'*P0(:,k),1,M)/sum(P0(:,k))^2-v0(:,k)'/sum(P0(:,k)))';
    end

%% itertion
while iter < maxiter 
    iter = iter + 1;
    % search for the next point
    % [P0,W0,v0] = main_loop(P0,DP,);
    x1 = cell(K,1); y1 = cell(K,1);
    dt = .01;
    P1 = P0 - dt*gradF;
    for j = 1:M
        P1(j,:) = projsplx(P1(j,:));
    end
    w1 = P1./repmat(sum(P1),M,1);
    W1 = zeros(K,1);
    for k = 1:K  % solve W1 for each k
        b{k}(N(k)+1:end) = w1(:,k);
        [W1(k),x1{k},y1{k}] = PLP_EE(C{k}(:),A{k},b{k},[x0{k};y0{k}]);
    end
    F1 = .5*(N.^alpha)'*(W1.^2);
    while F1 > F0
        dt = dt/2; 
        if dt < tol
            P = P1;
            return
        end
        P1 = P0 - dt*gradF;
        for j = 1:M
            P1(j,:) = projsplx(P1(j,:));
        end
        w1 = P1./repmat(sum(P1),M,1);
        W1 = zeros(K,1);
        for k = 1:K  % solve W1 for each k
            b{k}(N(k)+1:end) = w1(:,k);
            [W1(k),x1{k},y1{k}] = PLP_EE(C{k}(:),A{k},b{k},[x0{k};y0{k}]);
        end
        F1 = .5*(N.^alpha)'*(W1.^2);
    end
    P0 = P1; 
    W0 = W1; F0 = F1;
    x0 = x1; 
    y0 = y1;
    for k = 1:K
        v0(:,k) = y0{k}(N(k)+1:end);
        gradF(:,k) = (N(k)^alpha)*W0(k)*...
            (repmat(v0(:,k)'*P0(:,k),1,M)-v0(:,k)'/sum(P0(:,k))');
    end
end
P = P0; 
end