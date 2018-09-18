clear all; clc ;close all;

n     =  30;
m     =  8;
start = -2;
stop  =  2;
x     =  linspace(-2, 2, 30);

eps = 1;
rng(1);
r = rand(1,n) * eps;
y_1 = x.*(cos(r+0.5*x.^3)+sin(0.5*x.^3));                     % equation (1)
y_2 = 4*x.^5 - 5*x.^4 - 20*x.^3 + 10*x.^2 + 40*x + 10 + r;    % equation (2)

%    making a column vector from y_2
%    AX = f
y1 = y_1 .';
y2 = y_2 .';

A = ones(n , m);
      
% making matrix A which is Vandermonde matrix

for j = 2 : m + 1
    for i = 1 : n
        
        A(i , j) = x (1, i)^(j-1);
        
    end
end

   

[Q , R] = qr(A, 0);
b2      = (Q') * y2;        %  for equation (2)
b1      = (Q') * y1;        %  for equation (1)

% ************* Separation of zero part of R. Actually what we ge from QR
% factorization is R = [R1 , 0] and for polynomial approximation, R1 is
% required
% This command separates zero part of R i.e.  "R1 = R(any(R,2),:)"

% ************ Using Back Substitution *********

   c2 = backsubst(R , b2);
   p2 = A * c2;              % approximation of equation (2)
   plot (x, p2)
   hold on
   plot (x, y2,'o');
   
   
   c1 = backsubst(R , b1);
   p1 = A * c1;              % approximation of equation (1)
   figure
   plot (x, p1)
   hold on
   plot (x, y1,'o');

% ************ Using Cholesky *********






