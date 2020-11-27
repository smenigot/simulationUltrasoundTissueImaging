function [kappas,taus] = En_fit_kt(f_low,f_high,att,n,...
				kappa_infty,c0,rho0,n_pts,n_iters)
% [kappas,taus] = fit_kt(f_low,f_high,att,n,...
%                          kappa_infty,c0,rho0,n_pts,n_iters)
%
% MATLAB function to fit two-mechanism relaxation terms to 
% a desired power law fit.
% 
% Martin Anderson
% Created:  1-31-00  Last modified:  10-11-00
%
% Accepts the lower frequency limit,    (f_low)   [Hz]
% the upper frequency limit,            (f_high)  [Hz]
% the relative bandwidth,               (rbw)
% the attenuation,                      (att)     [dB/cm-MHz]
% the frequency dependence,             (n)       [MHz^n]  
% the steady-state compressibility,     (kappa_infty) [1/(kg/s^2-m)]
% the steady-state speed of sound,      (c0)      [m/s]
% the steady-state density,             (rho0)    [kg/m^3]
% the number of points at which
% to evaluate the fit,                  (n_pts)
% and the number of iterations.         (n_iters)
%
% The number of relaxation mechanisms is fixed at 2.
% The maximum number of iterations is currently 5.
%
% N.B.: The calculation of wavespeed is currently 
% provided only for completeness (cost is minimal), 
% as the fit is currently based on attenuation only.
% 
% Returns vectors of 
% compressibility values                (kappas)  [1/(kg/s^2-m)]
% and time constants                    (taus)    [s]
% for the optimum fit with uniform weight 
% across the frequency band [f_low:f_high]
% Dep: None
%
% Based on equations (40), (44), and (45) of: 
% Nachman, A. I., J. F. Smith, and R. C. Waag, "An equation for
% acoustic propagation in inhomogeneous media with relaxation 
% losses," JASA 88, pp. 1584-1595, 1990.

%% Convert attenuation units from dB/cm to dB/m:
att = att*100;

%% Define frequency points at which to evaluate fit:
omegas = 2*pi*[f_low:(f_high-f_low)/(n_pts-1):f_high];
ref_att = att*(omegas/(2*pi*1e6)).^n;

min_err = [];

%% Define initial parameter "guesses" and multipliers which 
% define the search points about these on each iteration:
kappas = [kappa_infty; kappa_infty];
taus = [1 1];

mult = [ 10.^[-10:9/8:-1];
	 10.^[-1:0.25:1];
	 2.^[-1:0.25:1];
	 1.4.^[-1:0.25:1];
	 1.1.^[-1:0.25:1]];

%% Main loop:
% Take current "guess" of kappas and taus and generate a 
% set of search points around them defined by a row
% in the multiplier matrix.
% Evaluate fit to desired attenuation curve for every 
% permutation.  Adopt best fit as "guesses" for next
% iteration, which has a more narrow search region.

nepers_dB = 20*log10(exp(1)); % Nepers -> dB/cm conversion factor 

for c_iter = 1:n_iters
    % Generate set of search points:
    kappas_trial1 = mult(c_iter,:)*kappas(1);
    kappas_trial2 = mult(c_iter,:)*kappas(2);
    taus_trial1 = mult(c_iter,:)*taus(1);
    taus_trial2 = mult(c_iter,:)*taus(2);
    lk = length(kappas_trial1);
    lt = length(taus_trial1);
    for ck1 = 1:lk
        for ck2 = 1:lk
            for ct1 = 1:lt
                for ct2 = 1:lt
                    kappas = [kappas_trial1(ck1); kappas_trial2(ck2)];
                    taus = [taus_trial1(ct1); taus_trial2(ct2)];

                    [cs,alphas] = calc_curves(omegas,kappa_infty,rho0,...
                        kappas,taus);

                    fit_att = nepers_dB*alphas;
                    err = sqrt(mean((ref_att-fit_att).^2));

                    if isempty(min_err)
                        min_err = err;
                        min_index = [ck1,ck2,ct1,ct2];
                    end
                    if err < min_err;
                        min_err = err;
                        min_index = [ck1,ck2,ct1,ct2];
                    end

                end
            end
        end
    end

  kappas = [kappas_trial1(min_index(1)); kappas_trial2(min_index(2))];
  taus = [taus_trial1(min_index(3)); taus_trial2(min_index(4))];

  [cs,alphas] = calc_curves(omegas,kappa_infty,rho0,kappas,taus);

end


%% Subroutine to calculate c(omega) and alpha(omega):

function [cs,alphas] = calc_curves(omegas,kappa_infty,rho0,kappas,taus);

% For each frequency, evaluate c and alpha, given kappas and taus:
for cf = 1:length(omegas)
  omega = omegas(cf);
  term1 = sum(kappas./(1+taus.^2*omega^2));
  term2 = kappa_infty + sum(kappas./(1+taus.^2*omega^2));
  term3 = sum(kappas.*taus*omega./(1+taus.^2*omega^2));
  % Eq: (44):
  cs(cf) = sqrt(2/rho0)  / sqrt( kappa_infty + term1 + ...
  		           sqrt( term2^2 + term3^2));
  % Eq: (45):
  alphas(cf) = omega * sqrt(rho0/2) *  sqrt( -kappa_infty ...
     - term1 + sqrt( term2^2 + term3^2));
end
