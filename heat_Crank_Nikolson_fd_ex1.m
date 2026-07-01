clc;
clear;
fprintf("\n Finite-Difference method to solve u_t-u_xx=f(t,x), u(a)=Ua, u(b)=Ub, u(0,x)=g(x). :\n");

% Initialization
a=input("\n Enter the left end point of the spatial interval: ");
b=input("\n Enter the right end point of the sparial interval: ");
ya=input("\n Enter the value of the solution at left end point of the interval: ");
yb=input("\n Enter the value of the solution at right end point of the interval: ");
T=input("\n Enter the final time: ");
c=input("\n Enter the diffusion coefficient: ");
% M=input("\n Enter the maixmum number of iteration for associated Newton: ");
% Tol=input("\n Enter the tolerance of approx solution: ");
% h=input("\n Enter the step length of each subinterval: ");
 n=input("\n Enter the number of partitions of the space interval: ");
 m=input("\n Enter the number of partitions of the time interval: ");

% Nref=input("\n Enter the number of partitions of the interval for reference solution: ");
% href=(T-t0)/Nref;
href=(b-a)/(2^9);
% N=(T-t_0)/h;

% RHS of the heat euqation
 f=@(t,x) 0;
 
% Intial conditon
g=@(x) sin(pi.*x/2);

% Actual solution
 u=@(t,x) exp(-((pi^2/4)*t)).*sin (pi*x/2);


% Backward Finite Difefrence
 h=(b-a)/n;
 x=a:h:b;
 k=T/m;
 t=0:k:T;
 lam=c*k/h^2;
 A=diag((1+lam)*ones(n-1,1))-diag((lam/2)*ones(n-2,1),1)-diag((lam/2)*ones(n-2,1),-1);
 B=diag((1-lam)*ones(n-1,1))+diag((lam/2)*ones(n-2,1),1)+diag((lam/2)*ones(n-2,1),-1);
 disp(A);
 w=zeros(m+1,n+1);
 w(:,1)=0; % Dirichlet boundary at left end point
 w(:,n+1)=0; % Dirichlet boundary at right end point
 w(1,:)=g(x); % Initial Condition
 z=g(x(2:n))';
 for i=1:m
     z= A\B*z; % One can use Crouts method here
     w(i+1,2:n)=z;
 end

% Error computation
 phi=u(t(m+1),x);
 er=abs(phi-w(m+1,:));
 fprintf("\n Numerical values for exact and approximated soluton and corresponding erros:\n");
 D=['         x    ' '    FD sol(1st)  '  ' Exact Soln'  '   Error.(2nd)' ];
 disp(D); 
 disp(cell2mat(compose('%14.8f', [x; w(m+1,:); phi; er]')));

 surf(x, t', w);






 