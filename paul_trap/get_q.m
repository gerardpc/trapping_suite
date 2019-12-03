%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the q parameter from Mathieu's equation
% for given a, beta (characteristic exponent). Basically it gets the
% beta(a,q) function and inverts it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[q] = get_q(a, beta)
fun = @(q) Mathieu_characteristic_exp(a, q) - beta;
q = fzero(fun, [0, 0.9080463]);
end