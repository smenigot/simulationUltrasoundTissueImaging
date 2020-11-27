% function pulse = En_MakePulse( pulse )
% Fonction qui créé l'impulsion
%
% By Sébastien Ménigot (Inserm U930)
% sebastien.menigot@univ-tours.fr

function pulse = En_MakePulse( pulse )

%% Parameters
A0          = pulse.A;
fs          = pulse.fs;
Nc          = pulse.Nc;
Npoints     = length(pulse.bruit);

T_bruit  = 1.15*1e-6;       % length of the excitation
fs_bruit = Npoints/T_bruit; % sampling frequency
%% --- Time ---
dt = 1/fs;                % Sample interval
tp = pulse.ts;

tp_bruit = (tp(1):1/fs_bruit:tp(end));

%% Pulse Ref
fRef = 3e6;
tpRef = (-5000:4999) * 1/fs;
xRef = A0 .* exp(-( (pi*fRef*tpRef/Nc).^2) ) .*sin(2*pi*fRef.*tpRef);
PRef = sum(xRef.*xRef)/(length(xRef)/fs);

Gaussienne = exp(-( (pi*fRef*tp/Nc).^2) );

%% Samples of the signal
signal_bruit = zeros(size(tp_bruit));
l2signal_bruit = fix(length(signal_bruit)/2);
signal_bruit(l2signal_bruit-fix(Npoints/2):l2signal_bruit+fix(Npoints/2)-1) = pulse.bruit;

%% interpolation
signal_bruit_interpole = interp1((0:length(tp_bruit)-1)*(tp_bruit(2)-tp_bruit(1)),signal_bruit,(0:length(tp)-1)*dt,'spline');

%% Apodisation
signal_bruit_apodise = Gaussienne' .* signal_bruit_interpole;

%% Amplitude
A_bruit = sqrt((PRef * (length(signal_bruit_apodise)/fs)) / sum(signal_bruit_apodise.^2)); 
signal_bruit_amplifie = sign(A0) * A_bruit .* signal_bruit_apodise;

%% Transducer
signal_bruit_filtre = En_Transducteur(signal_bruit_amplifie,fs);

%% -- Pulse final ---
pulse.t = (0:length(signal_bruit_interpole)) * dt;
pulse.p = signal_bruit_filtre';

return