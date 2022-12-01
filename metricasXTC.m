%% C�lculo de m�tricas de XTC (Petrosino - Lizaso, 2022)
clear; clc
% Par�metros de desplazamiento y regularizaci�n
deltaLat=0.00; %desplazamiento lateral [m] (positivo derecha)
deltaFte=0.00; %desplazamiento frontal [m] (positivo hacia altavoces)
beta=0.005; % beta =0 implica matriz sin regularizaci�n
% Par�metros generales
c = 345; fmax=16000; f=10:10:fmax; nf=length(f);
% Par�metros distribuci�n de altavoces
D = 2; % Distancia desde el plano de altavoces al de los o�dos [m]
d0 = 0.22; % Distancia entre o�dos [m]
dp = 1.2; % Distancia entre altavoces [m]
d1 = sqrt(D^2+(dp/2-d0/2)^2); % caso sim�trico ipsilateral
d2 = sqrt(D^2+(dp/2+d0/2)^2); % caso sim�trico contralateral
dLL = sqrt((D-deltaFte)^2+(dp/2-d0/2+deltaLat)^2); % Altavoz L a o�do L
dLR = sqrt((D-deltaFte)^2+(dp/2+d0/2+deltaLat)^2); % Altavoz L a o�do R
dRL = sqrt((D-deltaFte)^2+(dp/2+d0/2-deltaLat)^2); % Altavoz R a o�do L
dRR = sqrt((D-deltaFte)^2+(dp/2-d0/2-deltaLat)^2); % Altavoz R a o�do R
t1=d1/c; t2=d2/c; % Retardos caso sim�trico
tLL=dLL/c; tLR=dLR/c; tRL=dRL/c; tRR=dRR/c; % Retardos con desplazamientos
% C�lculo de las matrices (transferencia)
H=zeros(2,2,nf); % Matriz Hb desplazada seg�n par�metros
Cbeta=zeros(2,2,nf); k=zeros(1,nf); % Matriz C regularizada y k T=zeros(2,2,nf); % Matriz T (transferencia total)
for naux = 1:nf
w=2*pi*f(naux);
    Haux0=[exp(-1i*w*t1)/d1 exp(-1i*w*t2)/d2; ...
         exp(-1i*w*t2)/d2 exp(-1i*w*t1)/d1];
    Haux=[exp(-1i*w*tLL)/dLL exp(-1i*w*tLR)/dLR; ...
   exp(-1i*w*tRL)/dRL exp(-1i*w*tRR)/dRR];
    H(:,:,naux)=Haux; % Matriz H del caso con desplazamiento
    % Cbeta es la matriz XTC con regularizaci�n beta
    % Cbeta se calcula con la matriz Haux0 (sin desplazamiento)
    Cbeta(:,:,naux)=inv((Haux0'*Haux0+beta*eye(2)))*Haux0';
    k(naux)=cond(Cbeta(:,:,naux),2);
    T(:,:,naux)=H(:,:,naux)*Cbeta(:,:,naux); % Caso con desplazamiento
end
% C�lculo de las m�tricas del sistema XTC
% EsiII ipsilateral altavoz L a o�do L
EsiII=abs(squeeze(T(1,1,:)));
% EsiX contralateral altavoz L a o�do R
EsiX=abs(squeeze(T(2,1,:)));
% Crosstalk ipsilateral/contralateral
XTC=EsiII./EsiX;
% Rango din�mico de altavoces
Sfase = squeeze(abs(Cbeta(1,1,:)+Cbeta(2,1,:)));
Scontrafase = squeeze(abs(Cbeta(1,1,:)-Cbeta(2,1,:)));
Smax=max(Sfase,Scontrafase);
%Graficos de rango din�mico de altavoces
figure(1)
referencia=min(Sfase);
semilogx(f,20*log10(Sfase/referencia),'--' ...
    ,f,20*log10(Scontrafase/referencia),'--' ...
    ,f,20*log10(Smax/referencia),'k' ...
    ,'linewidth',1)
axis([20 fmax -10 50]); xticks([20 100 1000 10000]); grid on
legend('S_{fase}','S_{contrafase}','S_{max}','location','best')
title(['Rango din�mico (altavoces) para \beta=' num2str(beta)])
xlabel('frecuencia (Hz)')
ylabel('dB relativos')
%Gr�ficos de Crosstalk, nivel ipsilateral y contralateral
figure(2)
semilogx(f,20*log10(XTC),'k',f,20*log10(EsiII),'--b',f,20*log10(EsiX),'--r')
axis([20 fmax -40 40]); xticks([20 100 1000 10000]); grid on
title(['\chi ideal (sim�trico) para \beta=' num2str(beta) ...
    '; \DeltaLat=' num2str(deltaLat*100) 'cm; ' ...
    '\DeltaFte=' num2str(deltaFte*100) 'cm'])
hold on
nfmin2a=find(20*log10(XTC)>=20,1); 
nfmin2b=nfmin2a+find(20*log10(XTC(nfmin2a+1:end))<=20,1); 
area(f(nfmin2a:nfmin2b),0*f(nfmin2a:nfmin2b)+40,-20 ...
    ,'EdgeColor','none','FaceColor','r','FaceAlpha',0.1)
hold off
legend('\chi (crosstalk)','Esi ipsilateral L\rightarrowL','Esi contralateral L\rightarrowR' ...
    ,'\chi>20 dB','location','best')
xlabel('frecuencia (Hz)')
ylabel('dB relativos')
%Graficos de condici�n de n�mero
figure(3)
semilogx(f,k,'k')
title(['condici�n de n�mero de la matriz C para \beta=' num2str(beta)])
grid on
axis([20 fmax 0 inf])
xticks([20 100 1000 10000]);
xlabel('frecuencia (Hz)')
ylabel('k')
