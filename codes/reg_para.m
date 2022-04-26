function f = reg_para(phi)
% compute the distance regularization term coefficient by equation (4.12)
[phi_x, phi_y] = gradient(phi);
s = sqrt(phi_x.^2 + phi_y.^2);
a = (s >= 1);
b = (s < 1);
dps = a .* (s + 1) .* (s - 1) + b * 2 .* (s - 1) .* (2 * s - 1);
dp2 = a .* (3 * s .* s-1) + b .* (12 * s .* s - 12 * s + 2);
mdp = max(abs(dps), [], 'all');
mp2 = max(abs(dp2), [], 'all');
f = 1 ./ (4 .* max(mdp, mp2) + eps);