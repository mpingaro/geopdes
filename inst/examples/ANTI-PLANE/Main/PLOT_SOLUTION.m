%% PLOT SOLUTION PACMAN
clear; clc; close all;

load leap_FG;
load leap_SG;
load leap_SG_COMP0_100;

plot(leap_FG(:,1), leap_FG(:,2),'-- k', 'linewidth',2 );
hold on
plot(leap_SG(:,1), leap_SG(:,2),'-- b', 'linewidth',2 );
hold on
plot(leap_SG_COMP(:,1), leap_SG_COMP(:,2),'- r', 'linewidth',2 );
hold off
gtext( '(only) First Gradient' , 'fontsize', 20, 'rotation', 0, 'color' , 'k');
gtext( '(only) Second Gradient' , 'fontsize', 20, 'rotation', 0, 'color' , 'b');
gtext('Complete' , 'fontsize', 20, 'rotation', 0, 'color' , 'r');

%legend '(only) First Gradient' '(only) Second Gradient' 'Complete' 
title('Solution on the Lip','fontsize',30)
xlabel('lip x \in [-1, 1]','fontsize',15)
ylabel('Transverse Displacement','fontsize',15)
grid on