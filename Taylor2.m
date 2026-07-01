clc;
clear;
fprintf("\n Euler method to solve y'=te^(3t)-2y, y(t0)=y0. :\n");

% Initialization
t0=input("\n Enter the initial time of the interval: ");
T=input("\n Enter the final time of the interval: ");
y0=input("\n Enter the intial value of the solution: ");
% h=input("\n Enter the step length of each subinterval: ");
N=input("\n Enter the number of partitions of the interval: ");
h=(T-t0)/N;
% Nref=input("\n Enter the number of partitions of the interval for reference solution: ");
% href=(T-t0)/Nref;
href=(T-t0)/(2^9);
% N=(T-t_0)/h;

% RHS side of the euqation
f=@(t,y) t.*exp(3*t)-2*y;

% Partial derivatives of f
f_t=@(t,y) (1+3.*t).*exp(3*t);
f_y=@(t,y) -2;

% Actual solution
g=@(t, t0, y0) exp(-2*(t-t0))*y0+((5*t-1).*exp(3*t))./25 -((5*t0-1)*exp(5*t0-2.*t))./25;

% Taylor method of order two
t=t0:h:T;
w=zeros(1,N+1);
w(1)=y0;
for i=1:N
    F=f(t(i), w(i))+(h/2)*(f_t(t(i), w(i))+f_y(t(i), w(i)).*f(t(i), w(i)));
    w(i+1)=w(i)+h*F;
end

% Error computation
 y=g(t,t0,y0);
 er=abs(y-w);
 fprintf("\n Numerical values for exact and approximated soluton and corresponding erros:\n");
 D=['         t     ' '     Exact(t)  ' '    Taylor soln.    ' '    Error.   ' ];
 disp(D);
 disp(cell2mat(compose('%14.4f', [t; y; w; er]')));
 
 
 % Plotting actual solution
 tref=t0:href:T;
 z=g(tref,t0,y0);
 plot (tref,z, '-k');
 hold on

 % Plotting Euler solution
 plot(t, y, 'ko');
 plot(t, w, '--m*');

 title('Taylor method of order two', 11)
 xlabel('t-->')
 ylabel('y(t)-->')
 % xlim([-0.2 1.2])
 % ylim([-.3 4])
legend('Exact Soln.','Data points on exact soln.','Taylor(order two) Approx. soln.', 'Location', 'northwest')
grid on