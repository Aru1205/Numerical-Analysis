clc;
clear;
fprintf("\n Shooting method to solve y''=p(t)y'+q(t)y+r(t), y(a)=ya, y(b)=yb. :\n");

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
u1=zeros(1,N+1);
u2=zeros(1,N+1);
v1=zeros(1,N+1);
v2=zeros(1,N+1);

% Intialization of two systems
u1(1)=ya;
u2(1)=0;
v1(1)=0;
v2(1)=1;

% RK algorithm of order 4 to solve two linear systems
t=a:h:b;
for i=1:N
    % solution of 1st system
    k11=u2(i);
    k12=p(t(i))*u2(i)+q(t(i))*u1(i)+r(t(i));
    k21=u2(i)+ (h/2)*k12;
    k22=p(t(i)+(h/2))*(u2(i)+(h/2)*k12)+q(t(i)+(h/2))*(u1(i)+(h/2)*k11)+r(t(i)+(h/2));
    k31=u2(i)+ (h/2)*k22;
    k32=p(t(i)+(h/2))*(u2(i)+(h/2)*k22)+q(t(i)+(h/2))*(u1(i)+(h/2)*k21)+r(t(i)+(h/2));
    k41=u2(i)+ h*k32;
    k42=p(t(i)+h)*(u2(i)+h*k32)+q(t(i)+h)*(u1(i)+h*k31)+r(t(i)+h);
    u1(i+1)=u1(i)+(h/6)*(k11 + 2*k21 + 2*k31 + k41);
    u2(i+1)=u2(i)+(h/6)*(k12 + 2*k22 + 2*k32 + k42);
    % Solution of 2nd syhstem
    J11=v2(i);
    J12=p(t(i))*v2(i)+q(t(i))*v1(i);
    J21=v2(i)+ (h/2)*J12;
    J22=p(t(i)+(h/2))*(v2(i)+(h/2)*J12)+q(t(i)+(h/2))*(v1(i)+(h/2)*J11);
    J31=v2(i)+ (h/2)*J22;
    J32=p(t(i)+(h/2))*(v2(i)+(h/2)*J22)+q(t(i)+(h/2))*(v1(i)+(h/2)*J21);
    J41=v2(i)+ h*J32;
    J42=p(t(i)+h)*(v2(i)+h*J32)+q(t(i)+h)*(v1(i)+h*J31);
    v1(i+1)=v1(i)+(h/6)*(J11 + 2*J21 + 2*J31 + J41);
    v2(i+1)=v2(i)+(h/6)*(J12 + 2*J22 + 2*J32 + J42);
end
% Computing the solution of the BVP
w=zeros(1,N+1);
w(1)=ya;
for i=2:N+1
    w(i)=u1(i)+((yb-u1(N+1))/v1(N+1))*v1(i);
end

% Error computation
 phi=g(t);
 er=abs(phi-w);
 fprintf("\n Numerical values for exact and approximated soluton and corresponding erros:\n");
 D=['         t    ' '    RK sol(1st) ' '  RK sol(2nd) ' ' Shooting sol' '   Exact Soln'  '   Error.(2nd)' ];
 disp(D); 
 disp(cell2mat(compose('%14.8f', [t; u1; v1; w; phi; er]')));
 
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
 ylabel('y(t)-->')
 % xlim([-0.2 1.2])
 % ylim([-.3 4])
legend('Exact Soln.','Data points on exact soln.','RK Approx. (1st) soln.', 'Location', 'northwest')
grid on