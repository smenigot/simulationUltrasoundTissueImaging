% solve_sps2d(med,iters,excit,attn,record)
%
% MATLAB function to implement 2-D solution of
% the nonlinear wave equation using a psuedospectral
% / staggered A-B method, including attenuation modeled
% with two relaxation mechanisms.
%
% This version accepts an arbitrary excitation (forcing)
% function(s) at arbitrary points on the 2-D grid.
%
% Anderson, M. E.
% A 2D nonlinear wave propagation solver written in open-source MATLAB code 2000
% Proceeding IEEE Ultrasonic Symposium, 2000, 1351-1354

function [temps,pression] = En_solve_sps2d(med,iters,excit,attn)

%%
nonlin_flag = (sum(med.ba(:)) > 0);

% Calculate time increment [s]
dt = med.dt;

% Calculation of pressure and velocity are "staggered",
% i.e. velocity at integer time steps [0 1 2...] * dt,
% and pressure at half-integer time steps [0.5 1.5 2.5 ...] * dt;
dt2 = dt/2;

% Calculate bulk modulus matrix:
K = med.rho.*med.c.^2.*ones(med.num_axpts,med.num_latpts);

% Define axes for output image:
xax = [-med.num_latpts/2+0.5:med.num_latpts/2-0.5]*med.dx;
zax = [-med.num_pml:med.num_axpts-med.num_pml-1]*med.dz;

% Define Perfectly Matched Boundary Layer (PML)
% Only one taper is used for both x and y axes.
% After Section IV of Yuan, et al., 1997:
alpha_max = -log(med.pml_ampred).*med.c ./ (K.*med.dx);
quadratic = ([med.num_pml:-1:1].'/med.num_pml).^2;

% After equation 4 in Wojcik, et al., 1997:
% Note the factor of -K actually cancels out.
% (Included here for clarity w.r.t. references.)
left_alpha = -K(:,1:med.num_pml) .*...
    alpha_max(:,1:med.num_pml) .*...
    repmat(quadratic.',med.num_axpts,1);
right_alpha = -K(:,med.num_latpts-med.num_pml+1:med.num_latpts) .*...
    alpha_max(:,med.num_latpts-med.num_pml+1:med.num_latpts) .*...
    repmat(fliplr(quadratic.'),med.num_axpts,1);
top_alpha = -K(1:med.num_pml,:) .*...
    alpha_max(1:med.num_pml,:) .*...
    repmat(quadratic,1,med.num_latpts);
bottom_alpha = -K(med.num_axpts-med.num_pml+1:med.num_axpts,:) .*...
    alpha_max(med.num_axpts-med.num_pml+1:med.num_axpts,:).*...
    repmat(flipud(quadratic),1,med.num_latpts);

% Define spatial frequency multipliers used to calculate
% spectral derivatives:
ps_axshifts = [[0:1:med.num_axpts/2-1] [-med.num_axpts/2:1:-1]].'*i*pi;
ps_axshifts = repmat(ps_axshifts,1,med.num_latpts);

ps_latshifts = [[0:1:med.num_latpts/2-1] [-med.num_latpts/2:1:-1]].'*i*pi;
ps_latshifts = repmat(ps_latshifts,1,med.num_axpts);

fft_axscale = 1/(med.num_axpts*med.dz/2);
fft_latscale = 1/(med.num_latpts*med.dx/2);

% Calculate Adams-Bashforth coefficients for ABS4
% integrator in Ghrist, et al. (2000):
ab1 = (dt/24) * 26; ab2 = (dt/24) * -5;
ab3 = (dt/24) * 4 ; ab4 = (dt/24) * -1;

%% Initialize pressure and velocity vectors:
p = zeros(med.num_axpts,med.num_latpts);
px = zeros(med.num_axpts,med.num_latpts);
pz = zeros(med.num_axpts,med.num_latpts);
vx = zeros(med.num_axpts,med.num_latpts);
vz = zeros(med.num_axpts,med.num_latpts);
dpx = zeros(med.num_axpts,med.num_latpts);
dpz = zeros(med.num_axpts,med.num_latpts);
dvx = zeros(med.num_axpts,med.num_latpts);
dvz = zeros(med.num_axpts,med.num_latpts);
div_u = zeros(med.num_axpts,med.num_latpts);

% Initialize relaxation mechanism state variables:
if attn.flag ~= 0
    S = zeros(med.num_axpts,med.num_latpts,4,5);
    dS = zeros(med.num_axpts,med.num_latpts,4,5);
    kt1 = attn.kappas1./attn.taus1;
    kt2 = attn.kappas2./attn.taus2;
    tk_term = kt1 + kt2;
end

% Initialize W and dW matrices:
% (These are 4-D work matrices that are
%  the size of the spatial domain in the
%  first two dimensions stacked as [vx vz px pz]
%  and d[vx vz px pz]/dt
%  across the third dimension for linear conditions,
%  [vx vz px pz div_ux div_uy] and
%  d[vx vz px pz div_ux div_uy]/dt for nonlinear,
%  with the 5 timesteps of these used in the
%  forward time integration stacked across the
%  fourth dimension.)

if nonlin_flag ~= 0
    w = zeros(med.num_axpts,med.num_latpts,6,5);
    dw = zeros(med.num_axpts,med.num_latpts,6,5);
else
    w = zeros(med.num_axpts,med.num_latpts,4,5);
    dw = zeros(med.num_axpts,med.num_latpts,4,5);
end

% Initialize time-step pointer vector:
% (When the work matrices W and dW are updated
%  with each time step, their contents from
%  previous four time-steps are actually left
%  "in place", and a set of pointers to the
%  time-steps are updated. This saves computation
%  time and overhead.)
tsp = [1:5];

%% Main calculation loop
% M=moviein(fix(iters/4));
% compteur = 0;
rec_count = 1;
npourcent = 1;
pourcent_old = 0;

for cr = 1:iters
    
%     if rem(fix(cr/iters*100),25)
%         if (fix(cr/iters*100) - pourcent_old) > 24
%             npourcent = 1;
%         end
%         if npourcent == 1
%             fprintf('%2.0f ',cr/iters*100-1)
%             pourcent_old = cr/iters*100;
%         end
%         npourcent = npourcent + 1;
%     end
    
    % Integer time step: calculation of velocity
    % Calculate spatial partial derivatives:
    % Apply forcing function at points beyond the attenuation taper:
    if cr<=size(excit.p,1)/2
        % The next loop is not efficient, but allows for completely
        % arbitrary excitation locations throughout the medium:
        % (Usually very little time is spent on excitation, anyway.)
        for cpt = 1:length(excit.ax)
            p(excit.ax(cpt),excit.lat(cpt)) = excit.p(cr*2-1,cpt);
        end
    end

    dpx = real(ifft(fft(p.').*ps_latshifts)).' *fft_latscale;
    dpz = real(ifft(fft(p).*ps_axshifts)) *fft_axscale;
    
%% Calculate relaxation mechanism state variables:
    % (Based on equations 17 and 18' in Yuan, et al., 1999)
    if attn.flag ~= 0
        sum1x = S(:,:,1,tsp(4))./attn.taus1 + S(:,:,2,tsp(4))./attn.taus2;
        sum2x = px*tk_term;
        dS(:,:,1,tsp(4)) = -S(:,:,1,tsp(4))./attn.taus1 + px.*kt1;
        dS(:,:,2,tsp(4)) = -S(:,:,2,tsp(4))./attn.taus2 + px.*kt2;

        sum1z = S(:,:,3,tsp(4))./attn.taus1 + S(:,:,4,tsp(4))./attn.taus2;
        sum2z = pz*tk_term;
        dS(:,:,3,tsp(4)) = -S(:,:,3,tsp(4))./attn.taus1 + pz.*kt1;
        dS(:,:,4,tsp(4)) = -S(:,:,4,tsp(4))./attn.taus2 + pz.*kt2;
    end

%% Calculate PML boundary region masks:
    pml_vxleft = left_alpha .* vx(:,1:med.num_pml);
    pml_vxright = right_alpha .* ...
        vx(:,med.num_latpts-med.num_pml+1:med.num_latpts);

    pml_vztop = top_alpha .* vz(1:med.num_pml,:);
    pml_vzbottom = bottom_alpha .* ...
        vz(med.num_axpts-med.num_pml+1:med.num_axpts,:);

    % Calculate ODE coefficients and apply PML tapers:
    dw(:,:,1,tsp(4)) = (-1./med.rho) .* dpx;
    dw(:,1:med.num_pml,1,tsp(4)) =  dw(:,1:med.num_pml,1,tsp(4)) + pml_vxleft;
    dw(:,med.num_latpts-med.num_pml+1:med.num_latpts,1,tsp(4)) = ...
        dw(:,med.num_latpts-med.num_pml+1:med.num_latpts,1,tsp(4)) + pml_vxright;

    dw(:,:,2,tsp(4)) = (-1./med.rho) .* dpz;
    dw(1:med.num_pml,:,2,tsp(4)) = dw(1:med.num_pml,:,2,tsp(4)) + pml_vztop;
    dw(med.num_axpts-med.num_pml+1:med.num_axpts,:,2,tsp(4)) = ...
        dw(med.num_axpts-med.num_pml+1:med.num_axpts,:,2,tsp(4)) + pml_vzbottom;


%% Integrate one timestep forward by 4th order Adams-Bashforth:
    w(:,:,1:2,tsp(5)) = w(:,:,1:2,tsp(4)) + ...
        (ab1*dw(:,:,1:2,tsp(4)) + ab2*dw(:,:,1:2,tsp(3)) + ...
        ab3*dw(:,:,1:2,tsp(2)) + ab4*dw(:,:,1:2,tsp(1)));

    if attn.flag ~= 0
        S(:,:,:,tsp(5)) = S(:,:,:,tsp(4)) + ...
            (ab1*dS(:,:,:,tsp(4))+ab2*dS(:,:,:,tsp(3))+ ...
            ab3*dS(:,:,:,tsp(2))+ab4*dS(:,:,:,tsp(1)));
    end

    % Update velocity vectors:
    vx = w(:,:,1,tsp(5));
    vz = w(:,:,2,tsp(5));

%% Half-integer time step: calculation of pressure

    % Store pressure waveform before timestep:

    % Apply forcing function at points beyond the attenuation taper:
    if cr<=size(excit.p,1)/2
        % The next loop is not efficient, but allows for completely
        % arbitrary excitation locations throughout the medium:
        % (Usually very little time is spent on excitation, anyway.)
        for cpt = 1:length(excit.ax)
            vz(excit.ax(cpt),excit.lat(cpt)) = excit.vz(cr*2,cpt);
            vx(excit.ax(cpt),excit.lat(cpt)) = excit.vx(cr*2,cpt);
        end
    end

    % Calculate spatial partial derivatives:
    dvx = real(ifft(fft(vx.').*ps_latshifts)).' *fft_latscale;
    dvz = real(ifft(fft(vz).*ps_axshifts)) *fft_axscale;

    % Calculate PML boundary region masks:
    pml_pleft = left_alpha .* px(:,1:med.num_pml);
    pml_pright = right_alpha .* ...
        px(:,med.num_latpts-med.num_pml+1:med.num_latpts);

    pml_ptop = top_alpha .* pz(1:med.num_pml,:);
    pml_pbottom = bottom_alpha .* ...
        pz(med.num_axpts-med.num_pml+1:med.num_axpts,:);

%% Calculate ODE coefficients and apply PML tapers:
    % Note: The variables px and py are an artificial decomposition
    % of the pressure p used to calculate PML boundary masks
    % and attenuation state variables along each axis.
    % To decompose p into px and py, as required, implies
    % linear conditions.
    % Under nonlinear conditions, the variable p cannot be
    % calculated by a superposition of px and py.
    % Thus the PML boundary masks and attenuation state variables
    % are calculated using linear forward integration, i.e.
    % the nonlinear contribution, if any, for the current time step
    % is not included.
    % Pressure, under nonlinear conditions, is calculated on the
    % basis of the variable div_u (u = particle displacement),
    % which is in turn calculated by the forward integration of
    % dvx and dvy. The key is that div_u is calculated
    % before the squaring operation in the B/A constitutive
    % equation below.

    if attn.flag ~= 0
        dw(:,:,3,tsp(4)) = -K .* (dvx - sum1x + sum2x);
    else
        dw(:,:,3,tsp(4)) = -K .* dvx;
    end
    dw(:,1:med.num_pml,3,tsp(4)) =  dw(:,1:med.num_pml,3,tsp(4)) + pml_pleft;
    dw(:,med.num_latpts-med.num_pml+1:med.num_latpts,3,tsp(4)) = ...
        dw(:,med.num_latpts-med.num_pml+1:med.num_latpts,3,tsp(4)) + pml_pright;

    if attn.flag ~= 0
        dw(:,:,4,tsp(4)) = -K .* (dvz - sum1z + sum2z);
    else
        dw(:,:,4,tsp(4)) = -K .* dvz;
    end
    dw(1:med.num_pml,:,4,tsp(4)) = dw(1:med.num_pml,:,4,tsp(4)) + pml_ptop;
    dw(med.num_axpts-med.num_pml+1:med.num_axpts,:,4,tsp(4)) = ...
        dw(med.num_axpts-med.num_pml+1:med.num_axpts,:,4,tsp(4)) + pml_pbottom;

    % Under nonlinear conditions, we also forward integrate dvx and
    % dvy to calculate dux and duy (derivatives of particle
    % displacement), which we sum below to form div_u.
    if nonlin_flag ~= 0
        if attn.flag ~= 0
            dw(:,:,5,tsp(4)) = dvx - sum1x + sum2x;
        else
            dw(:,:,5,tsp(4)) = dvx;
        end
        dw(:,1:med.num_pml,5,tsp(4)) =  dw(:,1:med.num_pml,5,tsp(4)) + ...
            pml_pleft./-K(:,1:med.num_pml);
        dw(:,med.num_latpts-med.num_pml+1:med.num_latpts,5,tsp(4)) = ...
            dw(:,med.num_latpts-med.num_pml+1:med.num_latpts,5,tsp(4)) + ...
            pml_pright./-K(:,med.num_latpts-med.num_pml+1:med.num_latpts);
        if attn.flag ~= 0
            dw(:,:,6,tsp(4)) = dvz - sum1z + sum2z;
        else
            dw(:,:,6,tsp(4)) = dvz;
        end
        dw(1:med.num_pml,:,6,tsp(4)) = dw(1:med.num_pml,:,6,tsp(4)) + ...
            pml_ptop./-K(1:med.num_pml,:);
        dw(med.num_axpts-med.num_pml+1:med.num_axpts,:,6,tsp(4)) = ...
            dw(med.num_axpts-med.num_pml+1:med.num_axpts,:,6,tsp(4)) + ...
            pml_pbottom./-K(med.num_axpts-med.num_pml+1:med.num_axpts,:);
    end

    % Integrate one timestep forward by 4th order Adams-Bashforth:
    if nonlin_flag ~= 0
        w(:,:,3:6,tsp(5)) = w(:,:,3:6,tsp(4)) + ...
            (ab1*dw(:,:,3:6,tsp(4)) + ab2*dw(:,:,3:6,tsp(3)) + ...
            ab3*dw(:,:,3:6,tsp(2)) + ab4*dw(:,:,3:6,tsp(1)));
    else
        w(:,:,3:4,tsp(5)) = w(:,:,3:4,tsp(4)) + ...
            (ab1*dw(:,:,3:4,tsp(4)) + ab2*dw(:,:,3:4,tsp(3)) + ...
            ab3*dw(:,:,3:4,tsp(2)) + ab4*dw(:,:,3:4,tsp(1)));
    end

    % Update displacement and pressure matrices:
    px = w(:,:,3,tsp(5));
    pz = w(:,:,4,tsp(5));

    % Under linear conditions, p is the superposition of px and pz.
    % Under nonlinear conditions, p is calculated from div_u.
    if nonlin_flag ~= 0
        div_u = w(:,:,5,tsp(5)) + w(:,:,6,tsp(5));
        p = -K.*(div_u+0.5*med.ba.*(div_u).^2);
    else
        p = px + pz;
    end

    % Roll pointers back one step:
    tsp = [tsp(2:5) tsp(1)];
    
    %% pressure
    temps(cr,1) = cr*dt;

    pression(cr,1) = p(8,8); % transducer
    pression(cr,2) = p(229,8); % tissue
    
    if isnan(p)
        error('NaN pressure matrix')
    end
    
    
end
