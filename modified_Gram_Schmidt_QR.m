clc;
clear;
% Taking the data inputs
A=input("Enter the mtrix of full column rank : ");
fprintf("The given matrix A is:\n ");
disp(A);
[m, n]=size(A);
Q=zeros(m,n);
R=zeros(n,n);
for k=1:n
    R(k,k)=norm(A(:,k));
    if R(k,k)==0
        k=n+1;
        fprintf("A is not of full column rank.\n");
    else
        Q(:,k)=(1/R(k,k))*A(:,k);
        for j=k+1:n
            R(k,j)=Q(:,k)'*A(:,j);
            A(:,j)=A(:,j)-R(k,j)*Q(:,k);
        end
    end  
end
fprintf("The matrix Q_1 in QR factorization is given as:\n ");
disp(Q);
fprintf("The matrix R in QR factorization is given as:\n ");
disp(R);