function[]= plot_func()
%Plot function with constrains

dx1= 0.01;
dx2= 0.01;

[x1m,x2m]= meshgrid(-2:dx1:2, -2:dx2:2);

f_val=NaN(size(x1m));

validpoints= (x2m <= (1-0.5*x1m)) &(x2m >= x1m) & (x2m >= -x1m);
f_val(validpoints)= 2*(x1m(validpoints)-5).^2 + (x2m(validpoints)-1).^2;

f= 2*(x1m-5).^2 + (x2m-1).^2;

line1= 1-0.5*x1m;
line2= x1m;
line3= -x1m;

% Surf plot the function - NaNs aren't displayed

% figure; surf(x1m, x2m, f_val,  'EdgeColor', 'none');
% title('Feasible Area on 3D map');

figure; hold on;
contour(x1m, x2m, f);
contour(x1m, x2m, f_val,'LineColor','k');
plot(x1m,line1,'k',x1m,line2,'k',x1m, line3,'k');
xlabel('x1');
ylabel('x2');
legend('Contour','Feasible Area','x2<=1-0.5x1,x2>=x1,x2>=x1');
title('Feasible Area on Contour Map');

end