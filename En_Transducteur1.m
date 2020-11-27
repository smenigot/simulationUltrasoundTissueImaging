% Transducer
% By Sébastien Ménigot (Inserm U930)
% sebastien.menigot@univ-tours.fr

function pulse = En_Transducteur(pulse,fs_pulse)

      
wb1=2*2.2e6/fs_pulse;
wb2=2*1e6/fs_pulse;
wh1=2*5.8e6/fs_pulse;
wh2=2*7e6/fs_pulse;
Wp = [wb1 wh1];
Ws = [wb2 wh2];
[n2,Wn2] = buttord(Wp,Ws,3,10);

[B2,A2]=butter(n2,Wn2);

%% Filtering
pulse=filtfilt(B2,A2,pulse);

return



