clc;
clear;
% Taking the data inputs
A=input("Enter the mtrix : ");
fprintf("The given matrix A is:\n ");
disp(A);
[m, n]=size(A);
Q=eye(m);
for k=1: min(n,m-1)
    a=A(k:m,k);
    e=eye(m+1-k,1);
    if norm(a)~=0
        if a(1)>=0
            alpha=1;
        else
            alpha=-1;
        end
        v=a-alpha*norm(a)*e;
        Htilde=eye(m+1-k)-(2/norm(v)^2)*(v*v');
        H=blkdiag(eye(k-1), Htilde);
        fprintf('The Houser matrix H(%d) is:\n ', k);
        disp(H);
        A=H*A;
        Q=H*Q;
    end
end
fprintf("The matrix Q in QR factorization is given as:\n ");
disp(Q');
fprintf("The matrix R in QR factorization is given as:\n ");
disp(A);