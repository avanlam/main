function [phi, P] = deim(uIn)
%DEIM code for generating DEIM interpolation points
% https://github.com/pehersto/adeim

[n,m] = size(uIn);
phi = zeros(m, 1);
[~, phi(1)] = max(abs(uIn(:, 1)));
e = eye(n);
P = e(:, phi(1));
U(:, 1) = uIn(:, 1);
for l = 2:m
    uL = uIn(:, l);
    c = (P'*U)\(P'*uL);
    r = uL - U*c;
    [~, phi(l)] = max(abs(r));
    U(:, l) = uL;
    P(:, l) = e(:, phi(l));
end
phi = sort(phi);
P = e(:, phi);