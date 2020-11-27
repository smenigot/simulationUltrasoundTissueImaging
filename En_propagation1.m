% Preparation of the propagation simulation
% By Sébastien Ménigot (Inserm U930)
% sebastien.menigot@univ-tours.fr

function [time, echo, lTS, map, attn] = En_propagation(pulse, map, attn)

%% Parameters defining the excitation:
fc = 4e6;         % Maximal Fundamental frequency	[Hz]
p0 = pulse.A;     % Maximal initial pressure amplitude  [Pa]

%% Parameters defining the grid and medium:
num_axpts = 2^8;  % Number of points in space axially
num_latpts = 2^4; % Number of points in space laterally
c0 = map.c0;	  % Speed of sound in medium    [m/s]
                  % (a c0 is required for each attenuation calculation)
rho0 = map.rho0;% Density of medium 	        [kg/m^3]
                  % (a rho0 is required for each attenuation calculation)
num_pml = 7;      % Number of points in the perfectly matched layer
pml_ampred = 0.025; % Amplitude reduction factor for each point in PML
nodes_lam = 6;    % Nodes per wavelength
stabf = 0.2;      % Courant stability factor
lambda = c0/fc;        % Wavelength 	           [m]
dx = lambda/nodes_lam; % Lateral space increment   [m]
dz = lambda/nodes_lam; % Axial space increment     [m]


%% Map of velocity of sound and density of tissue (to set)

ba = map.ba0+randn(num_axpts,num_latpts); % B/A nonlinearity parameter matrix (scalar if homogeneous)
c =  map.c0+randn(num_axpts,num_latpts); % Speed of sound matrix (scalar if homogeneous)    [m/s]
rho =  map.rho0+randn(num_axpts,num_latpts);  % Density matrix (scalar if homogeneous)           [kg/m^3]

map.ba = ba;
map.c = c;
map.rho = rho;

%% Define time steps over which initial pulse is defined
% and the forcing functions used to stimulate the domain:

dt = stabf*min([dx dz])/max(c(:)); % Timestep is implicitly defined

T2 = 1/(2*1e6) + 1/(2*1e6);
ts = [-4*T2:dt/2:4*T2].';

pulse.fs = 1/(dt/2);
pulse.ts = ts;   
pulse = En_MakePulse(pulse);

% Calculate the desired pressure function:
excit_p = repmat(imag(hilbert(pulse.p)),...
			1,num_latpts);
% excit_p = p0*excit_p/max(abs(hilbert(excit_p(:))));

rho_c = ones(num_axpts,num_latpts)./(rho.*c);
excit_lat = [1:num_latpts];
excit_ax = ones(size(excit_lat))*(num_pml+1);

% Calculate the associated velocity function(s):
excit_vz = excit_p*diag(rho_c(num_pml+1,excit_lat));
excit_vx = zeros(size(excit_p));

% Set velocity excitation to zero in this case:
% (Applying BOTH pressure and velocity produces
%  a superposition, i.e. a doubling in this case.)
excit_vz = zeros(size(excit_p));

%% Parameters defining the attenuation:
% (One section like this required for each unique attenuation value.)
if attn.flag_calcul

    att_flag  = 1;	  % Flag specifying whether to include attenuation
    proc_flag = 1;    % Flag specifying whether to calculate the

    att = 0.45;       % Attenuation out liver                [dB/cm-MHz]
    n = 1.05;          % Frequency dependence        [MHz^n]

    kappa_infty = 1/(rho0*c0^2);    % Steady-state compressibility [1/(kg/s^2-m)]
    n_iters = 5;                    % Number of iterations to fit relaxation
                                    % mechanism constants
    n_pts = 10;                     % Number of points over which to fit relaxation
                                    % mechanism constants

    % Calculate attenuation, if specified by att_flag:
    if att_flag ~= 0
        if proc_flag == 1
            [kappas,taus] = En_fit_kt(0.5*fc,2*fc,att,n,kappa_infty,c0,rho0,n_pts,n_iters);
        end
    end

    if att_flag ~= 0
        kappas1 = kappas(1); % kappa_1 matrix (scalar if homogeneous)
        kappas2 = kappas(2); % kappa_2 matrix (scalar if homogeneous)
        taus1 = taus(1);     % tau_1 matrix (scalar if homogeneous)
        taus2 = taus(2);     % tau_2 matrix (scalar if homogeneous)
    else
        kappas1 = [];
        kappas2 = [];
        taus1 = [];
        taus2 = [];
    end
    
    attn = 	struct('flag',att_flag,...
    'kappas1',kappas1,'kappas2',kappas2,...
    'taus1',taus1,'taus2',taus2, 'flag_calcul', attn.flag_calcul);

    attn.flag_calcul = 0;
end

%% Load up the structures to be passed to the calculation function:
med = struct('c',c,'rho',rho,'ba',ba,'num_pml',num_pml,...
    'pml_ampred',pml_ampred,'stabf',stabf,...
    'num_latpts',num_latpts,'num_axpts',num_axpts,...
    'dx',dx,'dz',dz,'dt',dt,...
    'rho_c',rho_c);

excit = struct('p0',p0,'p',excit_p,...
    'vx',excit_vx,'vz',excit_vz,...
    'lat',excit_lat,'ax',excit_ax);

%% Call the calculation function:

iters = round((num_axpts*dx/c0)/dt + length(ts));
iters = round(iters/2)*2;

[time_pass1,pressure_pass1] = En_solve_sps2d(med,iters,excit,attn);

%% Time
dt = time_pass1(2) - time_pass1(1);

%% Pression

time(:, 1) = (0:length(pressure_pass1(:, 1))-1) .* dt;
echo(:, 1) = En_Transducteur(pressure_pass1(:, 1),1/dt);
echo(1:length(ts)/2,1) = 0; % remove excitation
echo(length(ts)/2:end) = echo(length(ts)/2:end);

lTS = length(ts);
