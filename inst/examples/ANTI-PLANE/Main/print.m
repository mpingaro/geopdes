%% PLOT SOLUTION PACMAN
clear; clc; close all;

load leap_SG_COMP0_175;

lip = leap_SG_COMP0_175;

for i = 1:size(lip,1)
    fprintf('(%d, %d) \n', lip(i,1),lip(i,2) );
end