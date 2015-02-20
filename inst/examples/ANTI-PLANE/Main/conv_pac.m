% Convergence PacMan Radial

% Point x = -0.5 & y = 0 
x1 = [2, 8, 32, 50, 72, 128, 200, 338, 512];
y1 = [0.508992589486209, 0.505892077632951, 0.504043485287959,...
    0.498780000352767, 0.497400571697213 ,0.496598187367594,...
    0.496138113132429 , 0.495818964165051 ,0.495756887521127];

figure, plot(x1,y1, '-o r', 'linewidth', 2, 'markersize', 6);
title ('Convergence Analysis Of Point x = -0.5 y = 0','fontsize',25)
xlabel('number of elements','fontsize',20)
ylabel('Transverse Displacement','fontsize',20)
axis([0 520 0.49 0.52])
%axis equal
grid on

% Point x = -0.04 & y = 0  -- > Punto interno!!
x2 = [2, 8, 32, 50, 72, 128, 200, 338, 512];
y2 = [0.035719408227171, 0.032100930067478, 0.029783654151943,...
    0.026374740859087, 0.023460448012223, 0.019985416194288,...
    0.018630890464229, 0.018150202265510, 0.018130448854038];
figure, plot(x2,y2, '-o r', 'linewidth', 2, 'markersize', 6);
title ('Convergence Analysis Of Point x = -0.04 y = 0','fontsize',25)
xlabel('number of elements','fontsize',20)
ylabel('Transverse Displacement','fontsize',20)
axis([0 520 0.00 0.04])
%axis equal
grid on