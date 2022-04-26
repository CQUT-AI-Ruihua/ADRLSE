close all; clear; clc;


% read image
Img = imread('../images/3_2.bmp');
if size(Img, 3) > 1
    Img = rgb2gray(Img);	  % Grayscale the image when the image is RGB format
end
I = double(Img);

% parameters setting
timestep = 1;                 								% timestep
iter_inner = 20;              								% inner iteration number
iter_outer = 40;			  								% outer iteration number
alfa = 1;					  								% area term: direction of the contour evolution in equation (4.15)
lambda = 5.0;                 								% length term coefficient
epsilon = 1.5;                								% parameter of the Heaviside and Dirac functions
sigma = 1.0;                  								% standard deviation of Gaussian kernel
Img_smooth = imgaussfilt(I, sigma, 'FilterSize', 3);        % Gaussian bulr
[Ix, Iy] = gradient(Img_smooth);                            % compute the gradient of image
f = abs(Ix) + abs(Iy);										% compute the gradient magnitude
g = 2 ./ (1 + exp(f));										% edge indicator function in equation (4.3)

% initialize the level set function as a binary step function
c0 = 2;
[m ,n] = size(I);
initialLSF = c0 * ones(m, n);
initialLSF(5: m - 5, 5: n - 5) = -c0;
phi = initialLSF;

% iteration
for n = 1 : iter_outer
    mu = reg_para(phi);                   					% compute the adaptive distance regularization term coefficient in equation (4.12)
    phi = adrlse_edge(phi, g, mu, lambda, alfa, epsilon, timestep, iter_inner);		% evolution the level set function in equation (4.9)
    if mod(n, 2) == 0										% visiualization of the evolution
        figure(1)
        imagesc(Img); axis off; axis equal; colormap(gray);
		hold on;
		contour(phi, [0 0], 'r', 'linewidth', 2);
        iterNum = [num2str(n), ' iterations'];
        title(iterNum);
    end
end

% visiualization the segmentation results
close all;
figure(2);
imagesc(Img); axis off; axis equal; colormap(gray);
hold on;
contour(phi, [0 0], 'r', 'linewidth', 2);
hold on;contour(initialLSF, [0 0], 'y', 'linewidth', 2);
figure(3);
mesh(phi);