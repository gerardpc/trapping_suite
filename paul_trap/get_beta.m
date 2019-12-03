%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beta / mathieu characteristic exponent from a, q parameters
% Check http://dlmf.nist.gov/28.2#iii for details on the calculation
%
% WARNING: deprecated function. Use instead Mathieu_characteristic_exp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
function[beta] = get_beta(a, q)

%% Approximations
% beta1 = sqrt(a + q^2/2); % Dehmelt approximation
% beta2 = sqrt(a - (a-1)*q^2/(2*(a-1)^2-q^2) - (5*a + 7)*q^4/(32*(a-1)^3*(a-4)) - (9*a^2 + 58*a + 29)*q^6/(64*(a-1)^5*(a-4)*(a-9))); 
%% Get value from Mathieu function
% this should be quite exact as long as q < 0.908. When getting closer the 
% error increases (for 0.908 it's about 2% with 1e-2). For more accuracy
% reduce step size. 
[~, ~, m_cos] = mathieu_functions(0, pi, 1e-3, a, q); 
beta = acos(m_cos(end))/pi;
end