clear all
close all
clc

%   caricamento dei dati
load('101m.mat')

x=val(1,:);
fc=360;

%   plot del segnale
figure
plot(x)
title('ECG prima derivazione')
xlabel('Campioni')

%   provo a calcolare la fft del segnale per vedere dove tagliare per 
%   rimuovere il trend
XF=fft(x-mean(x));
figure
f=linspace(0,fc/2,10800);
plot(f,abs(XF(1:(length(XF)/2))/length(x)))
title('Trasformata di Fourier del segnale')
xlabel('Frequenze')

%   filtro passa alto
[ord, ft]=buttord([4*2/fc, 0.9], [0.1, 4*2/fc], 0.5, 50); 

[b,a]=cheby2(10,100,4*2/fc, 'high');
figure
freqz(b,a,512,fc)

xf=filtfilt(b,a,x);
figure
subplot(2,1,1)
plot(xf)
title('Filtrato')
subplot(2,1,2)
plot(x)
title('Non filtrato')

%   stima spettrale del segnale originale
x=x-mean(x);
w=rectwin(length(x));
NFFT=length(x);
[Pxx, f]=pwelch(x, w, 0, NFFT, fc);
figure
plot(f,Pxx/max(Pxx))
title('Stime spettrale del segnale non flitrato')

%   stima spettrale del segnale filtrato
xf=xf-mean(xf);
[Pxxf, ff]=pwelch(xf, w, 0, NFFT, fc);
figure
plot(ff,Pxxf/max(Pxxf))
title('Stime spettrale del segnale flitrato')

durata_segnale=length(x)/fc;
ris_teorica=1/durata_segnale;
ris_app=ris_teorica;

y1=zeros(1,length(3:(length(xf)-2)));
y2=zeros(1,length(3:(length(xf)-2)));
y3=zeros(1,length(3:(length(xf)-2)));

for i=3:(length(xf)-2)
    y1(i)=abs(xf(i+1)-xf(i-1)); 
    y2(i)=abs(xf(i+2)-2*xf(i)+xf(i-2));
    y3(i)=1.3*y1(i)+1.1*y2(i);
end
figure
plot(y3(1:2000))

qrs=zeros(1,length(y3)-8);
for i=1:(length(y3)-8)
    %blocco=y3(((i-1)*8+2):(8*i));
    blocco=y3((i+1):(i+8));
    indici=find(blocco>=40);
    if length(indici)>=6
        %qrs((i-1)*8+1)=1;
        qrs(i)=1;
    else
        qrs(i)=0;
    end
end
figure
plot(qrs)
ind=find(qrs==1);
k=2;
QRS=[];
QRS(1)=ind(1);
for i=1:(length(ind)-1)
    diff=ind(i+1)-ind(i);
    if diff>1
        QRS(k)=ind(i+1);
        k=k+1;
    end
end
RR=[];
for i=1:(length(QRS)-1)
    RR(i)=QRS(i+1)-QRS(i);
end
figure
plot(RR/fc)
