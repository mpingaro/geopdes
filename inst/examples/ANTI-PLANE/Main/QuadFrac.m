clear all; close all; clc
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
KIII = 1; % Se varia cambiare : boundary_pacman_mixed_bc_h_drchlt  !!!!
mu = 1;
Ls = 0.1; 
Lt = 0.1;

C1 =Ls*Ls;
C2 =Lt*Lt/2;
% Definizione della geometria
p11 =[-1 0]; p12 =[1 0]; p22 =[1 1]; p21 =[-1 1];
srf = nrb4surf(p11,p12,p21,p22);

% Elevazione d'ordine
srf_e = nrbdegelev(srf,[2 2]);
% k-Raffinament
new_knots = linspace (0, 1, 21);  % Mettere dispari!!!!
new_knots = new_knots (2:end-1);
srf_r = nrbkntins(srf_e, {sort([new_knots,0.5,0.5]), new_knots});

problem_data.geo_name = srf_r;
% ----------------------------------------------------------------------- %
problem_data.nmnn_sides   = [1 2 3 4];
problem_data.drchlt_sides = [1 2 4];
% ----------------------------------------------------------------------- %
% Physical parameters
problem_data.c_diff  = @(x, y) mu*ones(size(x));        % Laplace
problem_data.d_diff  = @(x, y) C1*ones(size(x));        % Bilaplace
problem_data.copp = C2;
% Source and boundary terms (bulk)
%problem_data.f = @(x, y) zeros(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt_new;

% CHOICE OF THE DISCRETIZATION PARAMETERS
%clear method_data

method_data.degree     = [srf_r.dim-1 srf_r.dim-1];      % Degree of the splines
method_data.regularity = [srf_r.dim-2 srf_r.dim-2];      % Regularity of the splines
method_data.nsub       = [1 1];                          % Number of subdivisions
method_data.nquad      = [4 4];                          % Points for the Gaussian quadrature rule

% method_data.degree     = [4 4];                          % Degree of the splines
% method_data.regularity = [3 3];                          % Regularity of the splines
% method_data.nsub       = [20 20];                        % Number of subdivisions
% method_data.nquad      = [5 5];                          % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] = solve_fracture_antiPlane_Quad (problem_data, method_data);

% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, 501), linspace(0, 1, 501)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

%PLOT IN MATLAB.
% Plottagio geometria e control points
figure
nrbctrlplot(srf)
% Plottaggio raffinamento
figure
nrbkntplot(srf)
figure,
nrbctrlplot(srf_r)
figure
nrbkntplot(srf_r)
% Solution
figure
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution','fontsize',20)
axis equal

% ----------------------------------------------------------------------- %
% Estrazione valori sulla zona di frattura
leap_SG_COMP(:,1) = [X(:,1);X(501,:)'];
leap_SG_COMP(:,2) = [eu(:,1);eu(501,:)'];

% Soluzione 1 Gradiente
anal_sol = 2.*sqrt(-leap_SG_COMP(:,1)./(2*pi));

%
% Lip
d = Lt*sqrt(3)/(16*sqrt(2)) *sqrt((Lt/Ls)^4 -32*(Lt/Ls)^2 +128);
disp(d)

% Soluzione 2 Gradiente
CIII = 3.2;
%r = linspace(0,1);
r = 0:0.002:1;
s = CIII.*r.^(3/2).*( 3./ (16.*(Ls/Lt).^2-3) + 1 );

% ----------------------------------------------------------------------- %
figure
plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),'- b', 'linewidth',2 );
hold on
plot(-r,s, '- r', 'linewidth',2);
hold on
plot([-d -d],[0 0.8],'- k', 'linewidth',2 );
hold off

gtext( 'Limit zone of influence' , 'fontsize', 20, 'rotation', 0, 'color' , 'k');
gtext( 'of the second gradient (Analytic)' , 'fontsize', 20, 'rotation', 0, 'color' , 'k');
gtext( 'Numerical solution' , 'fontsize', 20, 'rotation', 0, 'color' , 'b');
gtext( 'Analytical solution' , 'fontsize', 20, 'rotation', 0, 'color' , 'r');

title ('Lip Solution (Zoom)','fontsize',20)
xlabel('lip x \in [-1, 1]','fontsize',15)
ylabel('Transverse Displacement','fontsize',15)
axis equal
grid on
save leap_SG_COMP