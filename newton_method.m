function [iter, x_opt, convergence] = newton_method(f, g, H, x, h)

% Initialization
epsilon= 1e-5;
beta= 0.5;  
sigma= 0.25;
a0= 1; 
iter= 0;

f_val= f(x);
g_val= g(x);
H_val= H(x);

convergence= struct;
convergence.x= zeros(size(x));convergence.f= 0;convergence.g_norm= 0;convergence.max_h= 0;

while norm(g_val) > epsilon
    iter = iter + 1;
    
    % find d by using cholseky factorization
    d= solve_for_d(g_val,H_val); 
    d= d/norm(d);

    % armijo rule
    alpha= a0;
    while (f(x + alpha*d) > f_val + alpha*g_val'*d*sigma)
        alpha = alpha*beta;
    end
    
    % update 
    x = x + alpha*d;
    
    f_val= f(x);
    g_val= g(x);
    H_val= H(x);
    
    m = length(convergence.f);
    if m < iter
        convergence.f(end+1 : end+2*m)= 0;
        convergence.g_norm(end+1 : end+2*m)= 0;
        convergence.max_h(end+1 : end+2*m)= 0;
        convergence.x(:, end+1 : end+2*m)= 0;       
    end
    
    convergence.x(:, iter)= x;
    convergence.f(iter)= f_val;
    convergence.g_norm(iter)= norm(g_val);
    convergence.max_h(iter)= max(h(x));
end

x_opt = x;
convergence.f= convergence.f(1:iter);
convergence.g_norm= convergence.g_norm(1:iter);
convergence.max_h= convergence.max_h(1:iter);
convergence.x= convergence.x(:,1:iter);

end


function [ d ] = solve_for_d( g,H )
    n=length(g);
    % cholesky factorization
    [L,diagonal]=mcholmz(H); % given to us
    D=diag(diagonal);
    % step 1: foward substitution
    y=zeros(n,1);
    y(1,1)=-g(1)/L(1,1);
    for i=2:n
        y(i)=((-g(i)-L(i,1:i-1)*y(1:i-1))/L(i,i));
    end
    % step 2
    m=D\y;
    % step 3: back substitution
    L_t=L';
    d=zeros(n,1);
    d(n,1)=m(n)/L_t(n,n);
    for i=1:n-1 
        d(n-i)=(m(n-i)-L_t(n-i,n-i+1:n)*d(n-i+1:n))/L_t(n-i,n-i);
    end
end
