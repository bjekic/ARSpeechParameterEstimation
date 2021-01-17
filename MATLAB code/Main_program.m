clear all
close all
clc

% AR model parameters
a1=-2.22; a2=2.89; a3=-3.08; a4=3.27; a5=-2.77; a6=2.35; a7=-1.7; a8=0.75;
A=[1 a1 a2 a3 a4 a5 a6 a7 a8];
B=1;
length_signal=256;
f_pitch=125; % pitch frequency of the synthetized signal


%% Synthetize signal
% Strube's glottal waveform
Tp=8e-3; Ts=3.2e-3; Tn=1.2e-3; Tp_odb=80; T_sample=1e-4; % parameters of the strube wave excitation
Glotal=zeros(Tp_odb,1);
for m=1:Tp_odb
    if m*T_sample<Ts
        Glotal(m)=(sin(pi*m*T_sample/2/Ts))^2;
    elseif m*T_sample<Ts+Tn
        Glotal(m)=cos(pi*(m*T_sample-Ts)/2/Tn);
    else
        Glotal(m)=0;
    end          
end
Glotal=Glotal*ones(1,20); Glotal=Glotal(1:1100);
ug_prim=([Glotal 0]-[0 Glotal]); ug_prim=ug_prim(1:end-1);
ug_sek=([ug_prim 0] - [0 ug_prim]); ug_sek=ug_sek(1:end-1);
Strube=ug_sek(46:1045); Strube=Strube./max(Strube);
s_strube=filter(B,A,Strube); s_pom=s_strube;
s_strube_orig=s_strube; % Synthetized signal wihout noise
eps=0.05; %Probability of outlier occurence
sig1=2e-4; sig2=100*sig1;
pom2=randn(1,1000); pom1=rand(1,1000);

for m=1:length(s_strube)
    if pom1(m)<eps
        pom=sig2*pom2(m);
    else
        pom=sig1*pom2(m);
    end
    s_pom(m)=s_strube(m)+pom;
end
s_strube=s_pom; % Synthetized signal with noise

% Train of Dirac pulses
Ts = 1e-4; Tp = 8e-3; Tp_odb = Tp/Ts; % parameters of the train of dirac pulses excitation

Dirac=zeros(1,1000); visina=1;
Dirac(1)=visina;
for m=1:length(Dirac)
    if ~mod(m,Tp_odb)
        Dirac(m)=visina;
    end
end
s_dirac=filter(B,A,Dirac);
s_dirac_orig=s_dirac; % Synthetized signal wihout noise

for m=1:length(s_dirac)
    if pom1(m)<eps
        pom=sig2*pom2(m);
    else
        pom=sig1*pom2(m);
    end
    s_pom(m)=s_dirac(m)+pom;
end
s_dirac=s_pom; % Synthetized signal with noise

%% Plotting the excitations

figure(1)
set(0,'DefaultAxesTitleFontWeight','normal');
subplot(2,1,1)
plot(Dirac,'k')
xlabel('a)')
ylabel('e[n]','Rotation',0)
yticks([0 1]);

subplot(2,1,2)
plot(Strube,'k'),ylim([1.2*min(Strube) 1])
xlabel({'b)'}), ylabel('e[n]','Rotation',0)
yticks([0 1]);
%% Evaluating the AR model parameters using LSQ, DUTTER and WLSQ algorithms
s=s_dirac_orig;
A_LSQ=zeros(length(s),8); A_DUTTER=zeros(length(s),8); A_WLSQ=zeros(length(s),8);
for m=length_signal:length(s)
   A_LSQ(m,:)=LSQ(s(m-length_signal+1:m)',8); %acquired using windowing
   [A_DUTTER(m,:), d]=Dutter(s(m-length_signal+1:m)',8,A_LSQ(m,:)');
   A_WLSQ(m,:)=WLSQ(s(m-length_signal+1:m)',8,A_DUTTER(m,:)',d);
end
A_LSQ1=A_LSQ; A_DUTTER1=A_DUTTER; A_WLSQ1=A_WLSQ;

s=s_dirac;
A_LSQ=zeros(length(s),8); A_DUTTER=zeros(length(s),8); A_WLSQ=zeros(length(s),8);
for m=length_signal:length(s)
   A_LSQ(m,:)=LSQ(s(m-length_signal+1:m)',8); %acquired using windowing
   [A_DUTTER(m,:), d]=Dutter(s(m-length_signal+1:m)',8,A_LSQ(m,:)');
   A_WLSQ(m,:)=WLSQ(s(m-length_signal+1:m)',8,A_DUTTER(m,:)',d);
end
A_LSQ2=A_LSQ; A_DUTTER2=A_DUTTER; A_WLSQ2=A_WLSQ;

s=s_strube_orig;
A_LSQ=zeros(length(s),8); A_DUTTER=zeros(length(s),8); A_WLSQ=zeros(length(s),8);
for m=length_signal:length(s)
   A_LSQ(m,:)=LSQ(s(m-length_signal+1:m)',8); %acquired using windowing
   [A_DUTTER(m,:), d]=Dutter(s(m-length_signal+1:m)',8,A_LSQ(m,:)');
   A_WLSQ(m,:)=WLSQ(s(m-length_signal+1:m)',8,A_DUTTER(m,:)',d);
end
A_LSQ3=A_LSQ; A_DUTTER3=A_DUTTER; A_WLSQ3=A_WLSQ;

s=s_strube;
A_LSQ=zeros(length(s),8); A_DUTTER=zeros(length(s),8); A_WLSQ=zeros(length(s),8);
for m=length_signal:length(s)
   A_LSQ(m,:)=LSQ(s(m-length_signal+1:m)',8); %acquired using windowing
   [A_DUTTER(m,:), d]=Dutter(s(m-length_signal+1:m)',8,A_LSQ(m,:)');
   A_WLSQ(m,:)=WLSQ(s(m-length_signal+1:m)',8,A_DUTTER(m,:)',d);
end
A_LSQ4=A_LSQ; A_DUTTER4=A_DUTTER; A_WLSQ4=A_WLSQ;

%% Plotting the results
x_osa=length_signal+1:length(A_LSQ);
figure(2)
close(2)
figure(2), hold off
subplot(2,2,1)
plot(x_osa,A_LSQ1(length_signal+1:end,1),'b'), hold all, plot(x_osa,A_DUTTER1(length_signal+1:end,1),'r'), plot(x_osa,A_WLSQ1(length_signal+1:end,1),'k')
ylim([-2.3 -2.16]), xlim([length_signal+1 x_osa(end)]), xlabel('a)')

subplot(2,2,2)
plot(x_osa,A_LSQ2(length_signal+1:end,1),'b'), hold all, plot(x_osa,A_DUTTER2(length_signal+1:end,1),'r'), plot(x_osa,A_WLSQ2(length_signal+1:end,1),'k')
ylim([-2.3 -2.16]), xlim([length_signal+1 x_osa(end)]), xlabel('b)'), legend('LSQ','Dutter','RBLP','Location','Best')
subplot(2,2,3)
plot(x_osa,A_LSQ3(length_signal+1:end,1),'b'), hold all, plot(x_osa,A_DUTTER3(length_signal+1:end,1),'r'), plot(x_osa,A_WLSQ3(length_signal+1:end,1),'k')
ylim([-2.3 -2.16]), xlim([length_signal+1 x_osa(end)]), xlabel('c)')
subplot(2,2,4)
plot(x_osa,A_LSQ4(length_signal+1:end,1),'b'), hold all, plot(x_osa,A_DUTTER4(length_signal+1:end,1),'r'), plot(x_osa,A_WLSQ4(length_signal+1:end,1),'k')
ylim([-2.3 -2.16]), xlim([length_signal+1 x_osa(end)]), xlabel('d)')

sr1=mean(A_WLSQ1(length_signal+1:end,1)); var1=var(A_WLSQ1(length_signal+1:end,1));
disp(['1: mean value: ' num2str(sr1) ', variance: ' num2str(var1)]);
sr2=mean(A_WLSQ2(length_signal+1:end,1)); var2=var(A_WLSQ2(length_signal+1:end,1));
disp(['2: mean value: ' num2str(sr2) ', variance: ' num2str(var2)]);
sr3=mean(A_WLSQ3(length_signal+1:end,1)); var3=var(A_WLSQ3(length_signal+1:end,1));
disp(['3: mean value: ' num2str(sr3) ', variance: ' num2str(var3)]);
sr4=mean(A_WLSQ4(length_signal+1:end,1)); var4=var(A_WLSQ4(length_signal+1:end,1));
disp(['4: mean value: ' num2str(sr4) ', variance: ' num2str(var4)]);

figure(3)
close(3)
figure(3), hold off
subplot(2,1,1)
plot(x_osa,A_LSQ1(length_signal+1:end,1),'b'), hold all, plot(x_osa,A_DUTTER1(length_signal+1:end,1),'r'), plot(x_osa,A_WLSQ1(length_signal+1:end,1),'k--','LineWidth',1.2)
ylim([-2.3 -2.15]), xlim([length_signal+1 x_osa(end)]), xlabel('a)')
ylabel('a_1','Rotation',0)
yticks([-2.3 -2.22 -2.15]);

subplot(2,1,2)
plot(x_osa,A_LSQ2(length_signal+1:end,1),'b'), hold all, plot(x_osa,A_DUTTER2(length_signal+1:end,1),'r'), plot(x_osa,A_WLSQ2(length_signal+1:end,1),'k')
ylim([-2.3 -2.15]), xlim([length_signal+1 x_osa(end)]), xlabel('b)'), legend('LSQ','Dutter','RBLP','Location','Southwest')
ylabel('a_1','Rotation',0)
yticks([-2.3 -2.22 -2.15]);

figure(4)
close(4)
figure(4), hold off

subplot(2,1,1)
plot(x_osa,A_LSQ3(length_signal+1:end,1),'b'), hold all, plot(x_osa,A_DUTTER3(length_signal+1:end,1),'r'), plot(x_osa,A_WLSQ3(length_signal+1:end,1),'k')
ylim([-2.3 -2.15]), xlim([length_signal+1 x_osa(end)]), xlabel('a)')
ylabel('a_1','Rotation',0)
yticks([-2.3 -2.22 -2.15]);

subplot(2,1,2)
plot(x_osa,A_LSQ4(length_signal+1:end,1),'b'), hold all, plot(x_osa,A_DUTTER4(length_signal+1:end,1),'r'), plot(x_osa,A_WLSQ4(length_signal+1:end,1),'k')
ylim([-2.3 -2.15]), xlim([length_signal+1 x_osa(end)]), xlabel('b)'), legend('LSQ','Dutter','RBLP','Location','Southwest')
ylabel('a_1','Rotation',0)
yticks([-2.3 -2.22 -2.15]);



