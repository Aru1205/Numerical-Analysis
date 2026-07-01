clc;
clear;
fprintf("\n Runge-Kutta method of order 4 to solve \n");
fprintf("y_1'=-4y_1+3y_22+6,\n");
fprintf("y_2'=-2.4y_1+1.6y_2+3.6,\n");
fprintf("y_(t0)=y_10, y_2(t0)=y_20 :\n");

% Initialization
t0=input("\n Enter the initial time of the interval: ");
T=input("\n Enter the final time of the interval: ");
y10=input("\n Enter the intial value of the 1st conponent of the solution: ");
y20=input("\n Enter the intial value of the 2nd conponent of the solution: ");
% h=input("\n Enter the step length of each subinterval: ");
N=input("\n Enter the number of partitions of the interval: ");
h=(T-t0)/N;
% Nref=input("\n Enter the number of partitions of the interval for reference solution: ");
% href=(T-t0)/Nref;
href=(T-t0)/(2^9);
% N=(T-t_0)/h;

% RHS side of the euqation
f1=@(t,y1,y2) -4*y1-2*y2+cos(t)+4*sin(t);
f2=@(t,y1,y2) 3*y1+y2-3*sin(t);

% Actual solution
g1=@(t) 2*exp(-t)-2*exp(-2*t)+sin(t);
g2=@(t)  -3*exp(-t)+2*exp(-2*t);
% To solve two linear system
w1=zeros(1,N+1);
w2=zeros(1,N+1);

% Intialization of two systems
w1(1)=y10;
w2(1)=y20;

% RK algorithm of order 4 to solve two linear systems
t=t0:h:T;
for i=1:N
    k11=f1(t(i),w1(i),w2(i));
    k12=f2(t(i),w1(i),w2(i));
    k21=f1(t(i)+(h/2),w1(i)+(h/2)*k11,w2(i)+(h/2)*k12);
    k22=f2(t(i)+(h/2),w1(i)+(h/2)*k11,w2(i)+(h/2)*k12);
    k31=f1(t(i)+(h/2),w1(i)+(h/2)*k21,w2(i)+(h/2)*k22);
    k32=f2(t(i)+(h/2),w1(i)+(h/2)*k21,w2(i)+(h/2)*k22);
    k41=f1(t(i)+h,w1(i)+h*k31,w2(i)+h*k32);
    k42=f2(t(i)+h,w1(i)+h*k31,w2(i)+h*k32);
    w1(i+1)=w1(i)+(h/6)*(k11 + 2*k21 + 2*k31 + k41);
    w2(i+1)=w2(i)+(h/6)*(k12 + 2*k22 + 2*k32 + k42);
end

% Error computation
 u1=g1(t);
 u2=g2(t);
 er1=abs(u1-w1);
 er2=abs(u2-w2);
 fprintf("\n Numerical values for exact and approximated soluton and corresponding erros:\n");
 D=['         t     ' 'RK soln.(1st) ' '  RK soln.(2nd)  ' '  Error.(1st)  ' '  Error.(2nd)' ];
 disp(D); 
 disp(cell2mat(compose('%14.9f', [t; w1; w2; er1; er2]')));
 
 
 % Plotting actual solution
 tref=t0:href:T;
 z1=g1(tref);
 plot (tref,z1, '-k');
 hold on

 % Plotting RK solution
 plot(t, u1, 'ko');
 plot(t, w1, '--b*');

 title('Runge-Kutta method of order four', 11)
 xlabel('t-->')
 ylabel('y1(t)-->')
 % xlim([-0.2 1.2])
 % ylim([-.3 4])
legend('Exact Soln.','Data points on exact soln.','RK Approx. (1st) soln.', 'Location', 'northwest')
grid on