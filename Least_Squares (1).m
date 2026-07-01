clc;
clear;
% Taking the data inputs
t=input("Enter the x-coordinates of the data points as row vector: ");
y=input("Enter the y-coordinates of the data points as row vector: ");
N1=length(t);
N2=length(y);
p=0;
while (p==0)
    if N1~=N2
        fprintf("Please check the number of data points and re-enter.\n");
        t=input("Enter the x-coordinates of the data points as row vector: ");
        y=input("Enter the y-coordinates of the data points as row vector: ");
        N1=length(t);
        N2=length(y);
    else 
        m=N1;
        p=1;
        fprintf("\n The data is given in a table as: \n\n ");
        disp(cell2mat(compose('%11.5f', [t ; y])));
    end
end
d=input("Enter the degree of the polynomial to fit the above data: ");
n=d+1;
A=zeros(m,n); % Finding the coefficient matrix
for i=1:n
    A(:,i)=t'.^(i-1);
end
fprintf("Matrix A is given as: ");
disp(A);
b=A'*y'; % RHS of the normal system
fprintf("The augmented Normal system of equations associated to the given data is: \n");
disp([A'*A b]);
% Finding Cholesky decomposition of A^T*A
fprintf("Lower triangular matrix of the Cholesky decomposition of the normal system is given as: \n ");
a=A'*A;
L=zeros(n);
L(1,1)=sqrt(a(1,1));
for j=2:n
    L(j,1)=a(j,1)/L(1,1);
end
for i=2:n-1
    s=0;
    for k=1:i-1
        s=s+L(i,k)^2;
    end
    L(i,i)=sqrt((a(i,i)-s));
    for j=i+1:n
        s=0;
        for k=1:i-1
            s=s+L(j,k)*L(i,k);
        end
        L(j,i)=(a(j,i)-s)/L(i,i);
    end
end
s=0;
for k=1:n-1
    s=s+L(n,k)^2;
end
L(n,n)=sqrt((a(n,n)-s));
% L=chol(A'*A,'lower');
disp(L);
% Solving the system L*L^Tx=A^T*y
z=zeros(1,n);
x=zeros(1,n);
z(1)=b(1)/L(1,1);
for i=2:n
    s=0;
    for j=1:i-1
        s=s+L(i,j)*z(j);
    end
    z(i)=(b(i)-s)/L(i,i);
end
x(n)=z(n)/L(n,n);
for i=n-1:-1:1
    s=0;
    for j=i+1:n
        s=s+L(j,i)*x(j);
    end
    x(i)=(z(i)-s)/L(i,i);
end
fprintf("Solution of the Normal equation is given as: ");
% disp(x);
disp(cell2mat(compose('%14.7f', x)));
% disp(y'-A*x');
e=norm(y'-A*x',2);
E=e^2;
fprintf("Minimized Error for this approximation is: %5.10f \n\n", E);
% Plotting the best fit polynomial with given data set
grid on
plot(t, y, "k*");
hold on;
h=(min(t)-1):.01: (max(t)+1);
w=0;
for i=1:n
    w=w+x(i)*h.^(i-1);  
end
plot(h, w, "-");
title('Best fit plynomial of degree d for the given data set', 'Fontsize', 11)
xlabel('t --->')
ylabel('y --->')
xlim([(min(t)-1) (max(t)+1)])
ylim([(min(y)-1) (max(y)+1)])
grid on