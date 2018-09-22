function [] = assignment_1(m)
clc ;
close all;
n     =  30;
start = -2;
stop  =  2;
x     =  linspace(-2, 2, 30);

eps = 1;
rng(1);
r = rand(1,n) * eps;
y_1 = x.*(cos(r+0.5*x.^3)+sin(0.5*x.^3));                     % data set(1): equation (1)
y_2 = 4*x.^5 - 5*x.^4 - 20*x.^3 + 10*x.^2 + 40*x + 10 + r;    % data set(2): equation (2)

% ******** making a column vector from y_2 ****************

y1 = y_1 .';
y2 = y_2 .';

%***********************************************************
%      making matrix A which is Vandermonde matrix
%***********************************************************

A = ones(n , m);
for j = 2 : m 
    for i = 1 : n
        
        A(i , j) = x (1, i)^(j-1);
        
    end
end


[Q , R1] = qr(A, 0);  % QR decomposition. Please note that R1 is non-zero part of R such that R = [R1 , 0]

B = A * A';   
[L , D] = cholesky(A' * A);   
bb = A' * y2;
R2  = L * D^(1/2);

 c1 = backsubst(R1 , (Q') * y1);
   c2 = backsubst(R1 , (Q') * y2);

    cc1 = backsubst(R2' , fwrdsubst (R2 , A' * y1));
    cc2 = backsubst(R2' , fwrdsubst (R2 , A' * y2));


  
   
%**************************************************************
%                      Plotting Task (1)
%**************************************************************
%  equation (1)
   figure
   plot (x, A * c1);
   hold on
   plot (x, y_1 ,'o');
   title ('Task1: Approx of Polynomial #(1) via [QR] method')
   legend ('y1: approx via [QR]', 'y1 :data set')
% equation (2)
   figure
   plot (x, A * c2)
   hold on
   plot (x, y_2,'o');
   title ('Task1: Approx of Polynomial #(2) via [QR] method')
   legend ('y2: approx via [QR]', 'y2 :data set')
%**************************************************************
%                      Plotting Task (2)
%**************************************************************
%  equation (1)
    figure
    plot (x, A * cc1)
    hold on
    plot (x, y_1,'o');
    title ('Task2: Approx of Polynomial #(1) via [QR] method')
    legend ('y2: approx via Chol', 'y2 :data set')
%  equation (2)
    figure
    plot (x, A * cc2)
    hold on
    plot (x, y_2,'o');
    title ('Task2: Approx of Polynomial #(2) via [QR] method')
    legend ('y2: approx via Chol', 'y2 :data set')
end




   
%**************************************************************
%                   Cholesky Method:

% B X = A(T) y ----> L * D^(0.5) * D^(0.5) * L(T) * X = A(T) y
%            R1 * R1(T) X = A(T) y ----> R1 * U = A(T) y
%                      R1(T) X = U
%**************************************************************
function [L , D] = cholesky(B)   %R = L * D^(0.5)

    [n , m] = size (B);
    L = zeros (m);
    D = zeros (n);
    
    Bk = B;
    
    if (n == m)
        for k = 1 : n
            lk = Bk(: , k) / Bk(k , k);
            D(k , k) = Bk(k , k);
            Bk = Bk - D(k , k) * lk * lk';
            L(: , k) = lk;            
        end
    else 
        fprint ('Square Matrix is Required for Input of Cholesky Method')
    end
      
end

%**************************************************************
%                      Back Substitution
%**************************************************************
function x = backsubst(U , b)

% finds x where Ux = b for an upper triangular matrix U

[n , m]  = size (U);


% erz = 1 e-12        % required for checking upper-triangular matrix


 x  = zeros(n , 1);
if (n == m)
   
    x (n) = b(n) / U(n , m);
end

for k = (n-1) : -1 : 1
    x(k) = (b(k) - U(k , k+1 : n) * x(k+1 : n)) / U(k , k);
end
end

%**************************************************************
%                      Forward Substitution
%**************************************************************
function x = fwrdsubst (U, b)

[n , m]  = size (U);

if (n == m)
    
    x = zeros (n , 1);
    
    x(1) = b(1)/ U(1 , 1);
    
    for k = 2 : n
    
    %for j = 1 : (k - 1)
        
    x(k) = (b(k) - U(k , 1 : k) * x(1 : k)) / U(k , k);
    
    end
end
end