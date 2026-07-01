clc;
clear;
fprintf("\n Shooting method to solve y''=f(t,y,y'), y(a)=ya, y(b)=yb. :\n");

% Initialization
a=input("\n Enter the left end point of the interval: ");
b=input("\n Enter the right end point of the interval: ");
ya=input("\n Enter the value of the solution at left end point of the interval: ");
yb=input("\n Enter the value of the solution at right end point of the interval: ");
M=input("\n Enter the maixmum number of iteration for associated Newton: ");
Tol=input("\n Enter the tolerance of approx solution: ");
% h=input("\n Enter the step length of each subinterval: ");
N=input("\n Enter the number of partitions of the interval: ");

% Nref=input("\n Enter the number of partitions of the interval for reference solution: ");
% href=(T-t0)/Nref;
href=(b-a)/(2^9);
% N=(T-t_0)/h;

% Coefficients of the euqation
 f=@(t, y, dy) -dy.^2-y+log(t);

 fy=@(t, y, dy) -1;
 
 fdy=@(t, y, dy) -2*dy;


% Actual solution
 g=@(t) log(t);


% To solve two linear system
w1=zeros(1,N+1);
w2=zeros(1,N+1);

h=(b-a)/N;
TK=(yb-ya)/(b-a);
for k=1:M
% Intialization of two systems
w1(1)=ya;
w2(1)=TK;
u1=0;
u2=1;

% RK algorithm of order 4 to solve two linear systems
t=a:h:b;
for i=1:N
    % solution of 1st system
    k11=w2(i);
    k12=f(t(i), w1(i), w2(i));
    k21=w2(i)+ (h/2)*k12;
    k22=f(t(i)+(h/2), w1(i)+(h/2)*k11,w2(i)+(h/2)*k12);
    k31=w2(i)+ (h/2)*k22;
    k32=f(t(i)+(h/2), w1(i)+(h/2)*k21,w2(i)+(h/2)*k22);
    k41=w2(i)+ h*k32;
    k42=f(t(i)+h, w1(i)+h*k31,w2(i)+h*k32);
    w1(i+1)=w1(i)+(h/6)*(k11 + 2*k21 + 2*k31 + k41);
    w2(i+1)=w2(i)+(h/6)*(k12 + 2*k22 + 2*k32 + k42);
    % Solution of 2nd syhstem
    J11=u2;
    J12=fy(t(i), w1(i),w2(i))*u1+ fdy(t(i), w1(i),w2(i))*u2;
    J21=u2+ (h/2)*J12;
    J22=fy(t(i)+(h/2), w1(i),w2(i))*(u1+(h/2)*J11)+ fdy(t(i)+(h/2), w1(i),w2(i))*(u2+(h/2)*J12);
    J31=u2+ (h/2)*J22;
    J32=fy(t(i)+(h/2), w1(i),w2(i))*(u1+(h/2)*J21)+ fdy(t(i)+(h/2), w1(i),w2(i))*(u2+(h/2)*J22);
    J41=u2+ h*J32;
    J42=fy(t(i)+h, w1(i),w2(i))*(u1+h*J31)+ fdy(t(i)+h, w1(i),w2(i))*(u2+h*J32);
    u1=u1+(h/6)*(J11 + 2*J21 + 2*J31 + J41);
    u2=u2+(h/6)*(J12 + 2*J22 + 2*J32 + J42);
end
if abs(w1(N+1)-yb)< Tol
    k=M+1;
end

% Newton iteration
TK=TK-((w1(N+1)-yb)/u1);
end
% Error computation
 phi=g(t);
 er=abs(phi-w1);
 fprintf("\n Numerical values for exact and approximated soluton and corresponding erros:\n");
 D=['         t    ' '  Shooting sol' '   Exact Soln '   '   Error' ];
 disp(D); 
 disp(cell2mat(compose('%14.6f', [t; w1; phi; er]')));
 
 % Plotting actual solution
 tref=a:href:b;
 z=g(tref);
 plot (tref,z, '-k');
 hold on

 % Plotting RK solution
 plot(t, phi, 'ko');
 plot(t, w1, '--b*');

 title('Runge-Kutta method of order four', 11)
 xlabel('t-->')
 ylabel('y(t)-->')
 % xlim([-0.2 1.2])
 % ylim([-.3 4])
legend('Exact Soln.','Data points on exact soln.','RK Approx. (1st) soln.', 'Location', 'northwest')
grid on