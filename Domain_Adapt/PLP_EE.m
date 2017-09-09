function [f_opt,x_opt,y_opt] = PLP_EE(c,A,b,z0,tol,maxiter)

% Using event-driven algorithm to solve the perturbed linear 
% programming in standard form exactly
% z0 = [x0;y0] is the initial guess


if nargin < 3
    error('PLP_EE(c,A,b,z0)')
end
if nargin > 6
    error('PLP_EE(c,A,b,z0,tol,maxiter)')
end
[m,n] = size(A);
switch nargin
    case 3
        z0 = ones(m+n,1);
        tol = 1e-10;
        maxiter = 500;
    case 4
        tol = 1e-10;
        maxiter = 500;
    case 5
        maxiter = 500;
end

% initial setup
f_opt = NaN;
x_opt = NaN; y_opt = NaN;
z0 = z0(:);
x0 = z0(1:n); y0 = z0(n+1:end);
I0 = x0 > 0; x0_a = x0(I0);
A_a = A(:,I0);
xi0 = A'*y0 - c; xi0_a = xi0(I0);
eta0 = A*x0 - b; 
Lya = @(x) norm(x);
lya0 = Lya([xi0_a;eta0]);
dt0 = 1e3;
dt = dt0;
thresh = 1.1; gamma = 2;
k = 0;
while k<maxiter
    k = k+1;
    while true
        x1 = zeros(n,1); 
        I1 = I0; 
        [~,n_a] = size(A_a);
        z_hat = [eye(n_a),-dt*A_a';...
            dt*A_a,eye(m)]\...
            [x0_a-c(I0)*dt;y0+b*dt];
        x1(I1) = z_hat(1:n_a);
        y1 = z_hat(n_a+1:end);
        event = 0;
        for i = 1:n
            if x1(i) < 0
               event = event+1;
               I1(i) = false; x1(i) = 0;
            end
        end
        xi1 = A'*y1-c;
        eta1 = A*x1-b;
        for i = 1:n
            if I1(i) == false && xi1(i) >= 0
                event = event + 1;
                I1(i) = true;
            end
        end
        xi1_a = xi1(I1); 
        lya1 = Lya([xi1_a;eta1]);
        if event == 0
            if lya1 < tol
                x_opt = x1; y_opt = y1;
                f_opt = c'*x_opt;
                return                
            elseif any([x1;y1]>1e4)
                f_opt = inf;
                disp('infeas or ubd')
                return;
            else
                dt = dt*gamma;
                x0 = x1; y0 = y1;
                I0 = I1; 
                x0_a = x0(I0); 
                A_a = A(:,I0); lya0 = lya1;
            end
            break
        else
            if lya1 < lya0 * thresh % accept
                x0 = x1; y0 = y1;
                I0 = I1; 
                x0_a = x0(I0); 
                A_a = A(:,I0); lya0 = lya1;
                if dt == dt0
                    dt = 1;
                end
                break
            else % reject
                if dt == 1e6
                    dt = 1;
                else
                  dt = dt/gamma;
                end
            end
        end
    end
end
end