clc;
clear;
fprintf("\n Finite difference method to solve y''=p(t)y'+q(t)y+r(t), y(a)=ya, y(b)=yb. :\n");

% Initialization
a=input("\n Enter the left end point of the interval: ");
b=input("\n Enter the right end point of the interval: ");
ya=input("\n Enter the value of the solution at left end point of the interval: ");
yb=input("\n Enter the value of the solution at right end point of the interval: ");
% h=input("\n Enter the step length of each subinterval: ");
N=input("\n Enter the number of partitions of the interval: ");
h=(b-a)/N;
% Nref=input("\n Enter the number of partitions of the interval for reference solution: ");
% href=(T-t0)/Nref;
href=(b-a)/(2^9);
% N=(T-t_0)/h;

% Coefficients of the euqation
 p=@(t) 0;
 q=@(t) 4;
 r=@(t)  -4*t;

% Actual solution
 g=@(t) (exp(2)/(exp(4)-1))*(exp(2*t)-exp(-2*t))+t;

% To solve two linear system
w=zeros(N+1,1);
% Intialization of two systems
w(1)=ya;
w(N+1)=yb;

% RK algorithm of order 4 to solve two linear systems
t=a:h:b;
d=diag((2+h^2.*q(t(2:N))));
du=diag(1-(h/2).*p(t(2:N-1)),1);
dl=diag(1+(h/2).*p(t(3:N)),-1);
A=dl-d+du;
%  disp(A);

B=zeros(N-1,1);
B(1)=h^2*r(t(2))-(1+(h/2)*p(t(2)))*ya;
B(2:N-2)=h^2.*r(t(3:N-1));
B(N-1)=h^2*r(t(N))-(1-(h/2)*p(t(N)))*yb;
% disp(B);

% Solving Tridiagonal system: One can use Crouts method here to solve
 u=A\B;
 w(2:N)=u;

% Error computation
 phi=g(t);
 er=abs(phi-w');
 fprintf("\n Numerical values for exact and approximated soluton and corresponding erros:\n");
 D=['         t    ' '    FD sol(1st)  '  ' Exact Soln'  '   Error.(2nd)' ];
 disp(D); 
 disp(cell2mat(compose('%14.8f', [t; w'; phi; er]')));
 
 % Plotting actual solution
 tref=a:href:b;
 z=g(tref);
 plot (tref,z, '-k');
 hold on

 % Plotting RK solution
 plot(t, phi, 'ko');
 plot(t, w, '--r*');

 title('Runge-Kutta method of order four', 11)
 xlabel('t-->')
 ylabel('y1(t)-->')
 % xlim([-0.2 1.2])
 % ylim([-.3 4])
legend('Exact Soln.','Data points on exact soln.','FD Approx. (1st) soln.', 'Location', 'northwest')
grid on