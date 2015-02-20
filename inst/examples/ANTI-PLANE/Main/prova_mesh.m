clear all
close all
clc

% Cerchio!
knots = {[0 0 0 1 1 1],[0 0 0 1 1 1]};
coefs = zeros(4,3,3);
% coefs(:, 1, 1) = [-1, 0, 0, 1]; 
% coefs(:, 2, 1) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 3, 1) = [ 0, 1, 0, 1];
% 
% coefs(:, 1, 2) = [-1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 2, 2) = [0,   0, 0, 1];
% coefs(:, 3, 2) = [1/sqrt(2),   1/sqrt(2), 0, 1/sqrt(2)];
% 
% coefs(:, 1, 3) = [0, -1, 0, 1];
% coefs(:, 2, 3) = [1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 3, 3) = [1,  0, 0, 1];

% QUADRATO
% coefs(:, 1, 1) = [0, 0.0, 0, 1]; 
% coefs(:, 2, 1) = [0, 0.5, 0, 1];
% coefs(:, 3, 1) = [0, 1.0, 0, 1];
% 
% coefs(:, 1, 2) = [0.5, 0.0, 0, 1];
% coefs(:, 2, 2) = [0.5, 0.5, 0, 1];
% coefs(:, 3, 2) = [0.5, 1.0, 0, 1];
% 
% coefs(:, 1, 3) = [1, 0.0, 0, 1];
% coefs(:, 2, 3) = [1, 0.5, 0, 1];
% coefs(:, 3, 3) = [1, 1.0, 0, 1];

%% PACMAN
coefs(:, 1, 1) = [0.0, 0.0, 0.0, 1.0];
coefs(:, 2, 1) = [0.0, (sqrt(2)-1)*sqrt( 2+sqrt(2) )/2, 0.0, sqrt( 2+sqrt(2) )/2];
coefs(:, 3, 1) = [1-sqrt(2)/2, sqrt(2)/2, 0.0, 1.0];

coefs(:, 1, 2) = [1.0, 0.0, 0.0, 1.0];
coefs(:, 2, 2) = [1.0, 0.5, 0.0, 1.0];
coefs(:, 3, 2) = [sqrt(2)/2, 1.0, 0.0, sqrt(2)/2];

coefs(:, 1, 3) = [2.0, 0.0, 0.0, 1.0];
coefs(:, 2, 3) = [2*sqrt( 2+sqrt(2) )/2, (sqrt(2)-1)*sqrt( 2+sqrt(2) )/2, 0.0, sqrt( 2+sqrt(2) )/2];
coefs(:, 3, 3) = [1.0+sqrt(2)/2, sqrt(2)/2, 0.0, 1.0];
%%
srf = nrbmak(coefs, knots);
%srf = nrbtform(srf, vecscale([1 1 0]));

figure
nrbctrlplot(srf)
figure
nrbkntplot(srf)

srf_r = nrbdegelev(srf,[2 2]);

figure
nrbctrlplot(srf_r)
figure
nrbkntplot(srf_r)

%% Soluzione di secondo gradiente!
ls = 0.1;
lt = 0.1;
CIII = 1;
r = linspace(0,1);
s = CIII.*r.^(3/2).*( 3./ (16.*(ls/lt).^2-3) + 1 );
figure, plot(r,s)

d = lt*sqrt(3)/(16*sqrt(2))*sqrt((lt/ls)^4 -32*(lt/ls)^2 +128);
disp(d)