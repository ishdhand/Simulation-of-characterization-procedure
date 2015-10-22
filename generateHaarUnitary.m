%% Generate Haar Unitary 
% Haar-random unitary generator
% to sample $n \times n$  unitary matrix U from Haar measure.
% inputs: integer n, dimension of unitary matrix
% outputs: U nxn unitary matrix sampled from Haar-random distribution.

function U = generateHaarUnitary(n)
[Q,R] = qr((randn(n) + 1i*randn(n))/sqrt(2));
U = Q*diag(diag(R)./abs(diag(R)));
end