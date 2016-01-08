clear; 
close all; 
clc;
% DEFINIZIONE GEOMETRIA (PACMAN)
% Costanti di secondo gradiente
KII = 0.5;
G   = 1;
l = 0.1;
eta = 0.0;
x1 = 0.1*l;
x2 = 0.4*l;
npt = 500;
Radius = 10*l;
% ----------------------------------------------------------------------- %
Ls = l/sqrt(2);
Lt = sqrt(eta+1)*l;
C1 = Ls*Ls;
C2 = Lt*Lt/2 - C1;

% Definizione della geometria
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
% srf = nrbdegelev(srf,[0 1]);
% k-raffinament
new_knots = linspace (0, 1, 100);
%new_knots = [0,logspace(-3,0,10)];
%new_knots = [ones(1,10)-logspace(0, -3,10),1];

new_knots = new_knots (2:end-1);
srf_r = nrbkntins(srf, {new_knots, new_knots});

% Secondo modo---
% [~,~,new_knots] = kntrfine(srf.knots,[9 9],[2 2],[1 1]);
% srf_refined = nrbkntins(srf,new_knots);
% Save Geometry
%nrbexport(srf,'pacman.txt');
%problem_data.geo_name = 'pacman.txt';
problem_data.geo_name = srf_r;

% ----------------------------------------------------------------------- %
% Type of boundary conditions for each side of the domain
% 1 Lato arco destro            |---> u_d = u_far;
% 2 Lato orizzontale sinistro;  |---> u libero;
% 3 Lato orizzontale destro;    |---> u_d = 0;
% 4 Lato arco sinistro          |---> u_d = u_far;

problem_data.nmnn_sides  = 2;           % Define Neumann conditions
problem_data.drchlt_sides_u= [1, 3, 4];  % Define Dirichlet conditions u
problem_data.drchlt_sides_r= [1, 3, 4];  % Define Dirichlet condition du/dn
% ----------------------------------------------------------------------- %
% Physical parameters
problem_data.c_diff  = @(x, y) G.*ones(size(x));           % Laplace
problem_data.d_diff  = @(x, y) C1.*ones(size(x));          % Bilaplace
problem_data.copp = C2;                                    % Bordi
% Source and boundary terms (bulk)
problem_data.f = @(x, y) zeros(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @boundary_pacman_mixed_bc_h_drchlt;
problem_data.hr =@boundary_pacman_mixed_bc_g_drchlt;
 
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
method_data.degree     = [srf_r.dim-1 srf_r.dim-1];       % Degree of the splines
method_data.regularity = [srf_r.dim-2 srf_r.dim-2];       % Regularity of the splines
method_data.nsub       = [1 1];                           % Number of subdivisions
method_data.nquad      = [3 3];                           % Points for the Gaussian quadrature rule

% CALL TO THE SOLVER
[geometry, msh, space, u] =...
    solve_fracture_antiPlane_grad (problem_data, method_data);

% POST-PROCESSIN

% EXPORT TO PARAVIEW
output_file = 'FractureModeAntiPlane';
vtk_pts = {linspace(0, 1, npt), linspace(0, 1, npt)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% Solution
figure
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
surf (X, Y, eu)
title ('Numerical solution')
axis equal

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
% Plot solution in any direction
j=1;
np = 100 ;
theta = zeros(1,np+1);
w1 = zeros(1,np+1);
w2 = zeros(1,np+1);
dt = 1/np;
for pn=0:dt:1
    theta(j) = 4*atand(pn);
    
    if pn==0.0
        a1 = pn ;
        a2 = pn ;
        b1 = x1/Radius;
        b2 = x2/Radius;
    elseif pn==1.0
        a1 = 1 ;
        a2 = 1 ;
        b1 = x1/Radius ;
        b2 = x2/Radius ;
    else
        a1 = pn/Radius*x1 ;
        a2 = pn/Radius*x2 ;
        b1 = 1/Radius*x1 ;
        b2 = 1/Radius*x2 ;
    end
    
    [uel1, p1] = sp_eval (u, space, geometry, {a1,b1});
    [uel2, p2] = sp_eval (u, space, geometry, {a2,b2});
    
    w1(j) = uel1 ;
    w2(j) = uel2 ;
    j= j+1;
end

figure, plot(theta, w1*G*sqrt(pi/(2*l))/KII)
title ('FIGURE 11 a) Normalized angular variation of the out-of-plane displacement field w at distance 0.1l','fontsize',15)
xlabel('\theta','fontsize',20)
ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
axis([0 180 0 0.4])
grid on

figure, plot(theta, w2*G*sqrt(pi/(2*l))/KII)
title ('FIGURE 11 b) Normalized angular variation of the out-of-plane displacement field w at distance 0.4l','fontsize',15)
xlabel('\theta','fontsize',20)
ylabel('wG(\pi/2l)^{1/2}/K_{III}','fontsize',20)
axis([0 180 0 0.8])
grid on

% Valutazione flesso sulla lip
p=polyfit(leap_SG_COMP(npt:end,1),leap_SG_COMP(npt:end,2),9);
y=polyval(p,leap_SG_COMP(npt:end,1));
%figure, plot(leap_SG_COMP(npt:end,1),leap_SG_COMP(npt:end,2),leap_SG_COMP(npt:end,1),y)
figure, plot(leap_SG_COMP(npt:end,1),leap_SG_COMP(npt:end,2))

dp=polyder(p) ;
ddp=polyder(dp);
flex=roots(ddp);
disp(abs(max(flex)))
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