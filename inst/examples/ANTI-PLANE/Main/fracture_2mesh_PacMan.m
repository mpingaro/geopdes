clear all; close all; clc
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
KII    = 0.5;
G      = 1;
l      = 0.2;
eta    = -0.9;
x1     = 0.1*l;
x2     = 0.4*l;
npt    = 1000;
Ls     = l/sqrt(2);
Lt     = sqrt(eta+1)*l;
Radius = 40*l;
C1     = Ls*Ls;
C2     = Lt*Lt/2 - C1;

%% PACMAN
% knots = {[0 0 0 1 1 1],[0 0 0 1 1 1]};
% %
% coefs = zeros(4,3,3);
% coefs(:, 1, 1) = [-1.0, 0.0, 0.0, 1.0];
% coefs(:, 2, 1) = [-sqrt( 2+sqrt(2) )/2, (sqrt(2)-1)*sqrt( 2+sqrt(2) )/2, 0.0, sqrt( 2+sqrt(2) )/2];
% coefs(:, 3, 1) = [-sqrt(2)/2, sqrt(2)/2, 0.0, 1.0];
% %
% coefs(:, 1, 2) = [0.0, 0.0, 0.0, 1.0];
% coefs(:, 2, 2) = [0.0, 0.5, 0.0, 1.0];
% coefs(:, 3, 2) = [0.0, 1.0, 0.0, sqrt(2)/2];
% %
% coefs(:, 1, 3) = [1.0, 0.0, 0.0, 1.0];
% coefs(:, 2, 3) = [sqrt( 2+sqrt(2) )/2, (sqrt(2)-1)*sqrt( 2+sqrt(2) )/2, 0.0, sqrt( 2+sqrt(2) )/2];
% coefs(:, 3, 3) = [sqrt(2)/2, sqrt(2)/2, 0.0, 1.0];
% %%
% srf = nrbmak(coefs, knots);
% srf = nrbtform(srf, vecscale([Radius Radius 0]));


knots = {[0 0 0 1 1 1],[0 0 0 1 1 1]};
coefs = zeros(4,3,3);
coefs(:, 1, 1) = [1, 0, 0, 1]; 
coefs(:, 2, 1) = [0, 0, 0, 1]; %  0.5 prima
coefs(:, 3, 1) = [0, 0, 0, 1];
coefs(:, 1, 2) = [1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 2, 2) = [0, 0, 0, 1];
coefs(:, 3, 2) = [0, 0, 0, 1]; % -0.5 prima
coefs(:, 1, 3) = [0, 1, 0, 1]; 
coefs(:, 2, 3) = [-1/sqrt(2), 1/sqrt(2), 0, 1/sqrt(2)];
coefs(:, 3, 3) = [-1, 0, 0, 1];

srf = nrbmak(coefs, knots);
srf = nrbtform(srf, vecscale([Radius Radius 0]));

% Elevazione d'ordine
%srf_elev = nrbdegelev(srf,[2 2]);
% k-raffinament
%new_knots = linspace (0, 1, 21);
%new_knots = new_knots (2:end-1);
%srf_raf = nrbkntins(srf_elev, {new_knots,sort([new_knots,0.5,0.5,0.5])});
%srf_raf = nrbkntins(srf_elev, {new_knots,new_knots});

% Definisco la geometria
% problem_data.geo_name = srf_raf;
problem_data.geo_name = srf;

% Plottagio geometria e control points
figure
nrbctrlplot(srf)
figure
nrbkntplot(srf)
% Plottaggio raffinamento
% figure,
% nrbctrlplot(srf_raf)
% figure
% nrbkntplot(srf_raf)

% ----------------------------------------------------------------------- %
problem_data.nmnn_sides     = [];         % --> Lato frattura
problem_data.drchlt_sides_u = [1, 3, 4]; % 1 - 4 Arco, 3 lato incastrato
problem_data.drchlt_sides_r = [1, 3, 4];
% ----------------------------------------------------------------------- %
% Physical parameters
problem_data.c_diff  = @(x, y) G.*ones(size(x));          % Laplace
problem_data.d_diff  = @(x, y) C1.*ones(size(x));         % Bilaplace
problem_data.copp = C2;                                   % Bordi
% Source and boundary terms (bulk)
problem_data.f = @(x, y) zeros(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt;
problem_data.hr = @boundary_pacman_mixed_bc_hr_drchlt;
 
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

% method_data.degree     = [srf_raf.dim-1 srf_raf.dim-1];   % Degree of the splines
% method_data.regularity = [srf_raf.dim-2 srf_raf.dim-2];   % Regularity of the splines
% method_data.nsub       = [1 1];                           % Number of subdivisions
% method_data.nquad      = [srf_raf.dim srf_raf.dim];       % Points for the Gaussian quadrature rule

method_data.degree     = [4 4];   % Degree of the splines
method_data.regularity = [3 3];   % Regularity of the splines
method_data.nsub       = [120 120]; % Number of subdivisions
method_data.nquad      = [5 5];   % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_fracture_antiPlane_NURBS_iso (problem_data, method_data);

% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, npt), linspace(0, 1, npt)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

%PLOT IN MATLAB.
figure
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution')
axis equal
% % Solution
% figure
% plot(leap_SG_COMP(:,1),leap_SG_COMP(:,2),leap_SG_COMP(:,1),anal_sol)
% axis equal
% grid on
% print('lip1.eps')
% save leap_SG_COMP

% Estrazione valori sulla zona di frattura
leap_FG(:,1) = [X(:,1);X(npt,:)'];
leap_FG(:,2) = [eu(:,1);eu(npt,:)'];
figure
plot(leap_FG(:,1),leap_FG(:,2))
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Estrazione valori sulla zona di frattura
leap_SG_COMP(:,1) = X(npt,:);
leap_SG_COMP(:,2) = eu(npt,:); 

% Soluzione 1 Gradiente
anal_sol = (KII/G).*sqrt(-leap_SG_COMP(:,1).*(2/pi)) ;
%
% Soluzione 2 Gradiente
CIII = 2.8;
rr = linspace(0,1);
s = CIII.*rr.^(3/2).*( 3./ (16.*(Ls/Lt).^2-3) + 1 ) ;
%
% Lip
d = Lt*sqrt(3)/(16*sqrt(2)) *sqrt((Lt/Ls)^4 -32*(Lt/Ls)^2 +128) ;
disp(d)
% ----------------------------------------------------------------------- %
figure
plot(leap_SG_COMP(:,1)./l,leap_SG_COMP(:,2).*(G*sqrt(pi/(2*l))/KII),'- r', 'linewidth',2 );
title ('FIGURE 9 a) Variation of crack face sliding displacement w along the crack face','fontsize',15)
xlabel('x_{1}/l','fontsize',20)
ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
axis([-5 0 0 2.5])
grid on

figure
plot(log(-leap_SG_COMP(:,1)./l), log(leap_SG_COMP(:,2).*(G*sqrt(pi/(2*l))/KII)),'- r', 'linewidth',2 );
title ('FIGURE 9 b) Variation of crack face sliding displacement w along the crack face in logarithmic scales','fontsize',15)
xlabel('log_{10}(x_{1}/l)','fontsize',20)
ylabel('log_{10}[wG(\pi/2l)^{1/2}/K_{III}]','fontsize',20)
axis([-3 1 -3 1])
grid on