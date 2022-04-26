function phi = bolza_edge(phi, g, mu, lambda, alfa, epsilon, timestep, iter)
[gx, gy] = gradient(g);

for k=1:iter
    phi = NeumannBoundCond(phi);        % Neumann boundary condition
    [phix, phiy] = gradient(phi);		% compute the gradient of the level set function phi
    mag = sqrt(phix .^ 2 + phiy .^ 2);	% compute the gradient magtitude of the level set function phi
    Phix = phix ./ (mag + 1e-10);
    Phiy = phiy ./ (mag + 1e-10);
    
    curvature = div(Phix, Phiy);        % curvature of the level set function
    diracPhi = Dirac(phi, epsilon);		% Dirac function
    
    distRegTerm = distReg_p2(phi);      % distance regularization term
    areaTerm = diracPhi .* g;           % area term
    edgeTerm = diracPhi .* (gx .* Phix + gy .* Phiy) + diracPhi .* g .* curvature;      % edge term
    t = (sum(Heaviside(-phi, epsilon), 'all') - lambda * sum(g .* diracPhi .* mag, 'all')) / (sum(g .* Heaviside(-phi, epsilon), 'all')+1e-10);	% adaptive area term coefficient in equation (4.15)
    phi = phi + timestep * (lambda .* edgeTerm + alfa * t .* areaTerm) + mu .* distRegTerm;  % evolution
end


function f = distReg_p2(phi)
[phi_x, phi_y] = gradient(phi);
s = sqrt(phi_x .^ 2 + phi_y .^ 2) + 1e-10;
a = (s > 1);
b = (s <= 1) & (s >= 0);
dps = a .* (s + 1) .* (s -1) + b .* 2 .* (s - 1) .* (2 * s - 1);  % diffusion rate
f = div(dps .* phi_x - phi_x, dps .* phi_y - phi_y) + 4 * del2(phi);


function f = div(nx, ny)
[nxx, ~] = gradient(nx);
[~, nyy] = gradient(ny);
f = nxx + nyy;


function f = Dirac(x, epsilon)
f = (1 / 2 / epsilon) * (1 + cos(pi * x / epsilon));
b = (x <= epsilon) & (x >= -epsilon);
f = f .* b;

function f = Heaviside(x, epsilon)
f = (1 + x / epsilon + sin(pi * x/epsilon) / pi) / 2;
f(x > epsilon) = 1;
f(x < -epsilon) = 0;

function g = NeumannBoundCond(f)
[nrow, ncol] = size(f);
g = f;
g([1 nrow], [1 ncol]) = f([3 nrow-2], [3 ncol-2]);  
g([1 nrow], 2 : end-1) = f([3 nrow-2], 2 : end-1);          
g(2 : end - 1, [1 ncol]) = f(2 : end - 1, [3 ncol-2]);  