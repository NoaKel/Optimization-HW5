function [Total_NM_iter,convergence] = augmented_lagrangian(f, g, H, h, dh, ddh, x)
% f- objective function value
% g - objective function gradient
% H - objective function hessian
% h - constrains functions
% dh - constrians function graident
% ddh - constrains function hessian
% x - initial x


% Initializing
p=1;
alpha=3;
p_max=10;
lambda=ones(length(h(x)),1);
epsilon=10^-5;
diff=inf;
iter=0;
max_iter=10^5;

convergence=struct;convergence.x= [];convergence.f= [];convergence.g_norm= [];convergence.max_h= [];convergence.lambda= zeros(length(lambda), max_iter);
Total_NM_iter=zeros(1,max_iter);

while (diff > epsilon || max(h(x)) > epsilon)
    iter = iter + 1;
    
    % Augmented Lagrangian Function 
    ALf = @(x) f(x) + lambda'*phi(p, h(x));
    ALg = @(x) g(x) + dh(x)'*(lambda.*dphi(p ,h(x)));
    ALH = @(x) H(x) + dh(x)'*(diag(lambda).*diag(ddphi(p, h(x))))*dh(x) + sum( bsxfun(@times, reshape(dphi(p, h(x)), 1, 1, length(dphi(p, h(x)))), ddh(x)), 3);
          
    [NM_iter, x, NM_convergence] = newton_method(ALf, ALg, ALH, x, h);
    
    % update lambda with safe guard    
    lambda_new=lambda.*dphi(p, h(x));
    lambda=3*lambda.*(lambda_new > 3*lambda)+1/3*lambda.*(lambda_new < 1/3*lambda)+lambda_new.*(lambda_new>=1/3*lambda & lambda_new<=3*lambda);
    
    % update p
    p = min(p_max,alpha*p);
    
    m = length(NM_convergence.f);
    convergence.x(:, end+1 : end+m)= NM_convergence.x;
    convergence.f(end+1 : end+m)= NM_convergence.f;
    convergence.g_norm(end+1 : end+m)= NM_convergence.g_norm;
    convergence.max_h(end+1 : end+m)= NM_convergence.max_h;
    convergence.lambda(:, iter)= lambda;
    
    if iter>1
        Total_NM_iter(:,iter)=Total_NM_iter(:,iter-1)+NM_iter;
    else
        Total_NM_iter(:,iter)=NM_iter;
    end
    
    f_prev = NM_convergence.f(end);
    diff = abs(f(x) - f_prev);
end

    convergence.lambda= convergence.lambda(:,1:iter);
    Total_NM_iter= Total_NM_iter(:,1:iter);
    
end

%% quadratic logaritmic 
 
function phi = phi(p, x)
phi=zeros(length(x),1);

    for i=1:length(x)
        t=p*x(i);
        if t >= -0.5
            phi(i,1) = (1/p)*(0.5*t^2+t);
        else
            phi(i,1) = (1/p)*(-0.25*log(-2*t)-3/8);
        end
    end

end

function dphi = dphi(p, x)
dphi=zeros(length(x),1);

    for i=1:length(x)
        t=p*x(i);
        if t >= -0.5
           dphi(i,1) = t+1;
        else
           dphi(i,1) = -1/(t*4);
        end
    end
 
end

function ddphi = ddphi(p, x)

ddphi=zeros(length(x),1);

    for i=1:length(x)
        t=p*x(i);
        if t >= -0.5
             ddphi(i,1) = p;
        else
             ddphi(i,1) = p/(4*(t)^2);
        end
    end
 

end
