%% Augmented Lagrangian HW5
close all; clear all;

% plot function feasible area (with its constrains)
plot_func;

% Initialization
x_start = [-1;2];

% True values calculated analyticaly
x_true=[2/3;2/3];
f_true=37+(2/3);
lambda_true=[12;11+(1/3);0];

% Function
Q=[4,0;0,2];
d=[-20;-2];
e=51;

f= @(x) 0.5*x'*Q*x + d'*x + e;
g= @(x) Q*x + d;
H= @(x) Q;

% Constrains
A=[0.5,1;1,-1;-1,-1];
b=[-1;0;0];
 
h= @(x) A*x + b;
dh= @(x) A;
ddh= @(x) zeros(2, 2, 3);

%Augmented Lagrangian Method 
[Total_NM_iter,convergence] = augmented_lagrangian(f, g, H, h, dh, ddh, x_start);

%% Plotting:
n= length(convergence.f);
m= length(convergence.lambda);
f_diff= zeros(1,n);
x_diff= zeros(length(x_start),n);
lambda_diff= zeros(1,m);

for i=1:n
    f_diff(:,i)=abs(convergence.f(i)-f_true);
    x_diff(:,i)=norm(convergence.x(:,i)-x_true);
end

for i=1:m
    lambda_diff(:,i)=norm(convergence.lambda(:,i)-lambda_true);
end

% plot
figure
subplot(2,2,1);
semilogy(1:1:length(convergence.g_norm),convergence.g_norm);
title('Gradient of the Augmented Lagrangian aggregate')
xlabel('NM iterations');
ylabel('||g||');
axis tight;

subplot(2,2,2);
semilogy(1:1:length(convergence.max_h),convergence.max_h);
title('Maximal constraint violation ')
xlabel('NM iterations');
ylabel('maximal violation');
axis tight;

subplot(2,2,3)
semilogy(1:1:length(f_diff),f_diff);
title('Residual in the objective function |f(x) -f(x*)| ')
xlabel('NM iterations');
ylabel('|f(x) -f(x*)|');
axis tight;

subplot(2,2,4)
semilogy(1:1:length(x_diff),x_diff, Total_NM_iter,lambda_diff);
title('Distance to the optimal point ||x -x*||  and to the optimal multipliers ||lambda – lambda*|| ');
% legend('||x -x*||','||lambda – lambda*||');
xlabel('NM iterations');
ylabel('||x -x*||,||lambda – lambda*||');
axis tight;

