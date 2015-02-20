%% PLOT SOLUTION PACMAN
clear; clc; close all;

ls = 0.1;

% LOAD FILE
load leap_SG_COMP0_000;
load leap_SG_COMP0_025;
load leap_SG_COMP0_050;
load leap_SG_COMP0_075;
load leap_SG_COMP0_100;
load leap_SG_COMP0_125;
load leap_SG_COMP0_150;
load leap_SG_COMP0_175;
load leap_SG_COMP0_200;
load leap_SG_COMP0_225;
load leap_SG_COMP0_250;


y = [leap_SG_COMP0_000(500,2), leap_SG_COMP0_025(500,2), leap_SG_COMP0_050(500,2), leap_SG_COMP0_075(500,2),...
    leap_SG_COMP0_100(500,2), leap_SG_COMP0_125(500,2), leap_SG_COMP0_150(500,2), leap_SG_COMP0_175(500,2),...
    leap_SG_COMP0_200(500,2), leap_SG_COMP0_225(500,2), leap_SG_COMP0_250(500,2)]/ls;


x = [0.000, 0.025, 0.050, 0.075, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250]/ls;

% PLOT SOLUTION
plot( x, y,'- b', 'linewidth',2 );

% plot(leap_SG(:,1), leap_SG(:,2),'-- b', 'linewidth',2 );
% hold on
% 
% plot(leap_SG_COMP(:,1), leap_SG_COMP(:,2),'- r', 'linewidth',2 );
% hold off
% gtext( '(only) First Gradient' , 'fontsize', 20, 'rotation', 0, 'color' , 'k');
% gtext( '(only) Second Gradient' , 'fontsize', 20, 'rotation', 0, 'color' , 'b');
% gtext('Complete' , 'fontsize', 20, 'rotation', 0, 'color' , 'r');
%legend '(only) First Gradient' '(only) Second Gradient' 'Complete' 
%title('Variation of the solution near the apex of the lip','fontsize',25)
xlabel( texlabel('l_{t}/l_{s}'),'fontsize',20)
ylabel( texlabel('w / l_{s} '),'fontsize',20)
grid on