% Computation of the cost function
% By Sébastien Ménigot (Inserm U930)
% sebastien.menigot@univ-tours.fr

function [E,medium] = En_FindEnergy(pulse, medium);

%% Medium Parameter
if ~medium.flag
    attn.flag_calcul = true;
    map.flag = true;
else
    attn = medium.attn;
    attn.flag_calcul = false;
    
    map = medium.map;
    map.flag = false;
end

map.c0   = 1584;    % [m/s]    Mean Speed of sound
map.rho0 = 1060;    % [kg/m3]  Mean Density
map.ba0  = 6;       %          Mean Nonlinear Coefficient B/A
map.p0   = 1.013e5; % [Pa]     Ambient pressure

%% Propagation
[time1, echo1, lTS1, map1, attn1] = En_propagation(pulse, map, attn);

medium.flag = true;
medium.map = map1;
medium.attn = attn1;

pulse.A = -pulse.A;

[time2, echo2, lTS2, map2, attn2] = En_propagation(pulse, map, attn1);

%%
m = min([length(echo1) length(echo2)]);
scan_line_PI = echo1(1:m) + echo2(1:m);

%% Energie et CTR
dt = time1(2)-time1(1);

% Zone 1 (to set)
pointA =  length(echo2)*1/4;
pointB =  length(echo2)*3/4;

% Power of zone 1
sum_zone1 = scan_line_PI(pointA:pointB,1);
E.E_zone1 = sum(sum_zone1.^2)/(length(sum_zone1)*dt);
% Power of zone 1
sum_zone2 = scan_line_PI([1:pointA-1,pointB+1:length(scan_line_PI)],1);
E.E_zone2 = sum(sum_zone2.^2)/(length(sum_zone2)*dt);
% Contrast
E.CTR = abs(E.E_zone1 - E.E_zone2)/(E.E_zone1 + E.E_zone2);

%% Valeur a optimiser
E.optim = E.CTR;
fprintf('Contrast = %2.2f\n',E.CTR);
return

