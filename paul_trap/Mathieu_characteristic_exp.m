%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the characteristic exponent of a Mathieu equation of a, q
% parameters. See "Mathieu Functions and Numerical Solutions of 
% the Mathieu Equation" for details. The n = 15 number can be increased
% for further accuracy, but this seems to be more than enough and
% ultrafast.
%
% Inputs
% a := Mathieu eq. parameter
% q := Mathieu eq. parameter
%
% Outputs
% nu := characteristic exponent (unstable if imaginary part != 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nu] = Mathieu_characteristic_exp(a, q, varargin)
%% Generate matrix of recurrence equation
n = 15; % matrix dimension is odd by symmetry reasons (the matrix is infinite, n is the truncation parameter)
% Calculate upper and lower diagonal 
diags_matrix = zeros(1, n);
n_i = -(n-1); 

if a == 0
    nu_Delta = 1;
else 
    nu_Delta = 0;
end

for i = 1:n
   diags_matrix(i) =  q/((nu_Delta + n_i)^2-a);
   n_i = n_i + 2;
end
upper_diag = diags_matrix(1:end-1);
lower_diag = diags_matrix(2:end);

Delta_matrix = eye(n) + diag(upper_diag, 1) + diag(lower_diag, -1);

%% Calculate characteristic exponent
Delta = det(Delta_matrix);

if a == 0
    nu = 1/pi*acos(2*Delta - 1);
else
    nu = 2/pi*asin(sqrt(Delta*sin(pi/2*sqrt(a))^2));
end
end

