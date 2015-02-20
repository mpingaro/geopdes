clear all; close all; clc
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
KIII = 1; % Se varia cambiare : boundary_pacman_mixed_bc_h_drchlt  !!!!
mu = 1;
Ls = 0.0; 
Lt = 0.1;

C1 =Ls*Ls;
C2 =Lt*Lt/2;
% Definizione della geometria
%p11 =[-0.5 0]; p12 =[0.5 0]; p22 =[0.5 0.5]; p21 =[-0.5 0.5];
p11 =[-1 0]; p12 =[1 0]; p22 =[1 1]; p21 =[-1 1];
srf = nrb4surf(p11,p12,p21,p22);

% new_knots = linspace (0, 1, 11);  % Mettere dispari!!!!
% new_knots = new_knots (2:end-1);
% srf_r = nrbkntins(srf, {new_knots, new_knots});
% problem_data.geo_name = srf_r;

problem_data.geo_name = srf;
% ----------------------------------------------------------------------- %
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [ 2 4];
% ----------------------------------------------------------------------- %
% Physical parameters
problem_data.c_diff  = @(x, y) mu*ones(size(x));        % Laplace
problem_data.d_diff  = @(x, y) C1*ones(size(x));        % Bilaplace
problem_data.copp = C2;
% Source and boundary terms (bulk)
problem_data.f = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt_new;

% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data

% method_data.degree     = [srf_r.dim-1 srf_r.dim-1];      % Degree of the splines
% method_data.regularity = [srf_r.dim-2 srf_r.dim-2];      % Regularity of the splines
% method_data.nsub       = [1 1];                          % Number of subdivisions
% method_data.nquad      = [4 4];                          % Points for the Gaussian quadrature rule

method_data.degree     = [2 2];                           % Degree of the splines
method_data.regularity = [1 1];                           % Regularity of the splines
method_data.nsub       = [20 20];                         % Number of subdivisions
method_data.nquad      = [4 4];                           % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] = solve_fracture_antiPlane_Quad (problem_data, method_data);

% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, 100), linspace(0, 1, 100)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

%PLOT IN MATLAB.
% Plottagio geometria e control points
figure
nrbctrlplot(srf)
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
leap_SG_COMP(:,1) = [X(:,1);X(100,:)'];
leap_SG_COMP(:,2) = [eu(:,1);eu(100,:)'];
figure
anal_sol = 2.*sqrt(-leap_SG_COMP(:,1)./(2*pi));
plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),leap_SG_COMP(:,1),anal_sol)
save leap_SG_COMP
% ----------------------------------------------------------------------- %
%close all