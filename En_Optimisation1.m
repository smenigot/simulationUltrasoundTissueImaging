% Setup the optimization
%
% By Sébastien Ménigot (Inserm U930)
% sebastien.menigot@univ-tours.fr

close all;
clear;
clc;
warning off;

%% ---------------------- Optimization parameter -----------------------
optim.Nc        = 2.3;                              %       Cycle number in the excitation
optim.Nitermax  = 5000;                             %       Iteration number
optim.Amp       = [400]*1e3;                        % [Pa]  Pulse amplitude
optim.Npoint    = 40;
%% -------------------------------- Optimization --------------------------
k = 1;
dir = ['Results' num2str(k)];
mkdir(dir);
En_Algorithme(optim, dir);



