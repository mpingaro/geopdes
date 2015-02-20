%% PLOT SOLUTION PACMAN
clear; clc; close all;

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

% PLOT SOLUTION
plot( leap_SG_COMP0_000(:,1), leap_SG_COMP0_000(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_025(:,1), leap_SG_COMP0_025(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_050(:,1), leap_SG_COMP0_050(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_075(:,1), leap_SG_COMP0_075(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_100(:,1), leap_SG_COMP0_100(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_125(:,1), leap_SG_COMP0_125(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_150(:,1), leap_SG_COMP0_150(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_175(:,1), leap_SG_COMP0_175(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_200(:,1), leap_SG_COMP0_200(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_225(:,1), leap_SG_COMP0_225(:,2),'- r', 'linewidth',1 );
hold on
%
plot( leap_SG_COMP0_250(:,1), leap_SG_COMP0_250(:,2),'- r', 'linewidth',1 );
hold off

% plot(leap_SG(:,1), leap_SG(:,2),'-- b', 'linewidth',2 );
% hold on
% 
% plot(leap_SG_COMP(:,1), leap_SG_COMP(:,2),'- r', 'linewidth',2 );
% hold off
% gtext( '(only) First Gradient' , 'fontsize', 20, 'rotation', 0, 'color' , 'k');
% gtext( '(only) Second Gradient' , 'fontsize', 20, 'rotation', 0, 'color' , 'b');
% gtext('Complete' , 'fontsize', 20, 'rotation', 0, 'color' , 'r');
%legend '(only) First Gradient' '(only) Second Gradient' 'Complete'
gtext( 'l_{t} = 2.5 - 0.0' , 'fontsize', 20, 'rotation', 45,'color' , 'blue' );


%title('Solution on the Lip','fontsize',30)
xlabel('lip x \in [-1, 1]','fontsize',20)
ylabel('Transverse Displacement','fontsize',20)
grid off