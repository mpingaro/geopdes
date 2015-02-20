clear all; close all; clc
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
KII = 1;
mu = 0;
Ls = 1; Lt = 0.4;

C1 =Ls*Ls;
C2 =Lt*Lt/2;
Radius = 1;
% Definizione della geometria
% crv = nrbline([0 0],[Radius 0]);
% srf = nrbrevolve(crv,[0 0 0],[0 0 Radius],pi);
% % k-raffinament
% new_knots = linspace (0, 1, 20);
% new_knots = new_knots (2:end-1);
% srf_r = nrbkntins(srf, {new_knots, new_knots});
% 
% % Save Geometry
% %nrbexport(srf,'pacman.txt');
% %problem_data.geo_name = 'pacman.txt';
% problem_data.geo_name = srf_r;

%% Secondo Tipo di mesh
%knots = {[0 0 0 1 1 1],[0 0 0 1 1 1]};
%coefs = zeros(4,3,3);
% coefs(:, 1, 1) = [1, 0, 0, 1]; 
% coefs(:, 2, 1) = [0, 0, 0, 1]; %  0.5 prima
% coefs(:, 3, 1) = [0, 0, 0, 1];
% coefs(:, 1, 2) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 2, 2) = [0, 0, 0, 1]; %%
% coefs(:, 3, 2) = [0, 0, 0, 1]; % -0.5 prima %%
% coefs(:, 1, 3) = [0, 1, 0, 1]; 
% coefs(:, 2, 3) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 3, 3) = [-1, 0, 0, 1];

% Quadrato!
% knots = {[0 0 1 1],[0 0 1 1]};
% coefs = zeros(4,2,2);
% coefs(:, 1, 1) = [0, 0, 0, 1]; 
% coefs(:, 2, 1) = [0, 1, 0, 1];
% coefs(:, 1, 2) = [1, 0, 0, 1];
% coefs(:, 2, 2) = [1, 1, 0, 1];

% Cerchio!
% knots = {[0 0 0 1 1 1],[0 0 0 1 1 1]};
% coefs = zeros(4,3,3);
% coefs(:, 1, 1) = [-1, 0, 0, 1]; 
% coefs(:, 2, 1) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 3, 1) = [ 0, 1, 0, 1];
% coefs(:, 1, 2) = [-1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 2, 2) = [0,   0, 0, 1];
% coefs(:, 3, 2) = [1/sqrt(2),   1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 1, 3) = [0, -1, 0, 1];
% coefs(:, 2, 3) = [1/sqrt(2), -1/sqrt(2), 0, 1/sqrt(2)];
% coefs(:, 3, 3) = [1,  0, 0, 1];

%% PacMan 2
coefs = zeros(4,3,3);
knots = {[0 0 0 1 1 1],[0 0 0 1 1 1]};

coefs(:, 1, 1) = [-1, 0, 0, 1]; 
coefs(:, 2, 1) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 3, 1) = [ 0, 1, 0, 1];
%
coefs(:, 1, 2) = [-1, 0, 0, 1]; %
coefs(:, 2, 2) = [0, 0, 0, 1];
coefs(:, 3, 2) = [1/sqrt(2),   1/sqrt(2), 0, 1/sqrt(2)];
%
coefs(:, 1, 3) = [0, 0, 0, 1];
coefs(:, 2, 3) = [1, 0, 0, 1]; %
coefs(:, 3, 3) = [1,  0, 0, 1];
% 
srf = nrbmak(coefs, knots);
srf = nrbtform(srf, vecscale([Radius Radius 0]));
% %%

figure
nrbctrlplot(srf)
% Plottaggio raffinamento
figure
nrbkntplot(srf)

% Elevazione d'ordine
% srf = nrbdegelev(srf,[0 1]);
% k-raffinament

% new_knots = linspace (0, 1, 3);
% new_knots = new_knots (2:end-1);
% srf_r = nrbkntins(srf, {new_knots, new_knots});

%problem_data.geo_name = srf_r;
% figure
% nrbctrlplot(srf_r)
% Plottaggio raffinamento
% figure
% nrbkntplot(srf_r)

% keyboard

problem_data.geo_name = srf;
%%
% ----------------------------------------------------------------------- %
% Type of boundary conditions for each side of the domain
% 1 Lato arco destro            |---> u_d = u_far;
% 2 Lato orizzontale sinistro;  |---> u libero;
% 3 Lato orizzontale destro;    |---> u_d = 0;
% 4 Lato arco sinistro          |---> u_d = u_far;

problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [2 3];
% ----------------------------------------------------------------------- %
% Physical parameters
problem_data.c_diff  = @(x, y) mu.*ones(size(x));          % Laplace
problem_data.d_diff  = @(x, y) C1.*ones(size(x));          % Bilaplace
problem_data.copp = C2;                                    % Bordi
% Source and boundary terms (bulk)
problem_data.f = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt;
 
% Derivate prime delle normali.
problem_data.dnorm_x_x  = @test_circular_plate_couple_dnx_x;
problem_data.dnorm_x_y  = @test_circular_plate_couple_dnx_y;
problem_data.dnorm_y_x  = @test_circular_plate_couple_dny_x;
problem_data.dnorm_y_y  = @test_circular_plate_couple_dny_y;
%
problem_data.dnorm_x_xx = @test_circular_plate_couple_dnx_xx;
problem_data.dnorm_x_xy = @test_circular_plate_couple_dnx_xy;
problem_data.dnorm_x_yy = @test_circular_plate_couple_dnx_yy;
%
problem_data.dnorm_y_xx = @test_circular_plate_couple_dny_xx;
problem_data.dnorm_y_xy = @test_circular_plate_couple_dny_xy;
problem_data.dnorm_y_yy = @test_circular_plate_couple_dny_yy;

% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data

method_data.degree     = [4 4];       % Degree of the splines
method_data.regularity = [3 3];       % Regularity of the splines
method_data.nsub       = [10 10];     % Number of subdivisions
method_data.nquad      = [5 5];       % Points for the Gaussian quadrature rule

% method_data.degree     = [srf_r.dim-1 srf_r.dim-1];       % Degree of the splines
% method_data.regularity = [srf_r.dim-2 srf_r.dim-2];       % Regularity of the splines
% method_data.nsub       = [1 1];                           % Number of subdivisions
% method_data.nquad      = [4 4];                           % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_fracture_antiPlane_NURBS_iso (problem_data, method_data);

% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, 100), linspace(0, 1, 100)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

%PLOT IN MATLAB.
% % Plottagio geometria e control points
% figure
% nrbctrlplot(srf)
% % Plottaggio raffinamento

% figure
% nrbkntplot(srf)
% figure,
% nrbctrlplot(srf_r)
% figure
% nrbkntplot(srf_r)
% Solution
figure
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution')
axis equal
%axis tight

% ----------------------------------------------------------------------- %
% Estrazione valori sulla zona di frattura
leap_SG_COMP(:,1) = [X(1,:)'; X(:,100)];
leap_SG_COMP(:,2) = [eu(1,:)'; eu(:,100)]; 

anal_sol = 2.*sqrt(-leap_SG_COMP(:,1)./(2*pi));
%
figure
plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),leap_SG_COMP(:,1),anal_sol)
axis equal
grid on
print('lip1.eps')
save leap_SG_COMP
% ----------------------------------------------------------------------- %