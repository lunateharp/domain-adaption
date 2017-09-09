k = randi(K);
i = randi(N(k)); j = randi(M);
e_jk = zeros(M,K);
e_jk(j,k) = 1;
h = 1e-6;
P1 = P0 + h*e_jk;
w1 = P1./repmat(sum(P1),M,1);
W1 = zeros(K,1);
for k = 1:K  % solve W1 for each k
    b{k}(N(k)+1:end) = w1(:,k);
    W1(k) = PLP_EE(C{k}(:),A{k},b{k},[x0{k};y0{k}]);
end
F1 = .5*(N.^alpha)'*(W1.^2);
(F1-F0)/h
GP(j,k)