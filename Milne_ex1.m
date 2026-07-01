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

% Actual solution
  g=@(t, t0, y0) exp(-2*(t-t0))*y0+((5*t-1).*exp(3*t))./25 -((5*t0-1)*exp(5*t0-2.*t))./25;

% RK algorithm of order 4 used for intial steps
t=t0:h:T;
w=zeros(1,N+1);
w(1)=y0;
for i=1:N
    k1=f(t(i),w(i));
    k2=f(t(i)+(h/2),w(i)+(h/2)*k1);
    k3=f(t(i)+(h/2),w(i)+(h/2)*k2);
    k4=f(t(i)+h,w(i)+h*k3);
    w(i+1)=w(i)+(h/6)*(k1 + 2*k2 + 2*k3 + k4);
end

% Intialization using exact values of solution
% for i=1:4   
%    w(i)=g(t(i),t0,y0);
% end

% Milnes method 4-step
w(1)=y0;
for i=4:N
    w(i+1)=w(i-3)+(4*h/3)*(2*f(t(i),w(i))-f(t(i-1),w(i-1))+2*f(t(i-2),w(i-2)));
end

% Error computation
 y=g(t,t0,y0);
 er=abs(y-w);
 fprintf("\n Numerical values for exact and approximated soluton and corresponding erros:\n");
 D=['         t     ' '     Exact(t)  ' '    Milne soln.    ' '    Error.   ' ];
 disp(D);
 disp(cell2mat(compose('%14.7f', [t; y; w; er]')));
 
 
 % Plotting actual solution
 tref=t0:href:T;
 z=g(tref,t0,y0);
 plot (tref,z, '-k');
 hold on

 % Plotting RK solution
 plot(t, y, 'ko');
 hold on
 plot(t, w, '--r*');

 title('Milnes method', 11)
 xlabel('t-->')
 ylabel('y(t)-->')
 % xlim([-0.2 1.2])
 % ylim([-.3 4])
legend('Exact Soln.','Data points on exact soln.','Mile Approx. soln.', 'Location', 'northwest')
grid on