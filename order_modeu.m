clc;
clear;
fprintf("\n Euler method to solve y'=te^(3t)-2y, y(t0)=y0. :\n");

% Initialization
t0=input("\n Enter the initial time of the interval: ");
T=input("\n Enter the final time of the interval: ");
y0=input("\n Enter the intial value of the solution: ");
% h=input("\n Enter the step length of each subinterval: ");
N=input("\n Enter the numbers of partitions of the interval as a vector: ");

% Different step sizes
h=(T-t0)./N;

% Nref=input("\n Enter the number of partitions of the interval for reference solution: ");
% href=(T-t0)/Nref;
href=(T-t0)/(2^9);
% N=(T-t_0)/h;

% RHS side of the euqation
f=@(t,y) t.*exp(3*t)-2*y;

% Actual solution
g=@(t, t0, y0) exp(-2*(t-t0))*y0+((5*t-1).*exp(3*t))./25 -((5*t0-1)*exp(5*t0-2.*t))./25;

% Perform iteration for different step sizes
L=length(h);
ERR=zeros(1,L);
for n=1:L
    t=t0:h(n):T;
    w=zeros(1,N(n)+1);
    w(1)=y0;

    % Modified Euler iteration
    for i=1:N(n)
        k1=f(t(i),w(i));
        k2=f(t(i)+h(n),w(i)+h(n)*k1);
        w(i+1)=w(i)+(h(n)/2)*(k1 + k2);
    end
    
    % Error computation
    y=g(t,t0,y0);
    er=abs(y-w);
    ERR(n)=max(er);
end

%Plotting order of convergence
grid on
loglog(h,ERR, '--ko')
hold on
loglog(h, h.^2, '-b')
grid on

title('Order of convergence for Modified Euler method', 11)
 xlabel('h-->')
 ylabel('Global Error(h)-->')
 % xlim([-0.2 1.2])
 % ylim([-.3 4])
legend('Order of convergence line','Reference line of slope 2', 'Location', 'northwest')