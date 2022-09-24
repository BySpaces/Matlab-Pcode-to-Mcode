clc;
clear;
closeall;


c=3e8;
fc=3e9;
lambda=c/fc;
Ralpha=1.2;
Aalpha=1.5;
H=8000-50;
Yc=6000;
R0=sqrt(H^2+Yc^2);
Xmin=-200;
Xmax=200;
Yw=300;
La=1.5;
theta=0.886*lambda/La;
Ls=R0*theta;
Vr=150;
Ts=Ls/Vr;
Xwid=Ls+Xmax-Xmin;
Twid=Xwid/Vr;
Ka=2*Vr^2/lambda/R0;
Ba=abs(Ka*Ts);
PRF=Aalpha*Ba;
PRT=1/PRF;
dx=PRT;
N=ceil(Twid/dx);
N=2^nextpow2(N);
x=linspace((Xmin-Ls/2)/Vr,(Xmax+Ls/2)/Vr,N);

X=Vr*x;
PRT=Twid/N;
PRF=1/PRT;
dx=PRT;
Tr=5e-6;
Br=100e6;
Kr=Br/Tr;
Fs=Ralpha*Br;
dt=1/Fs;
Rmin=sqrt((Yc-Yw)^2+H^2);
Rmax=sqrt((Yc+Yw)^2+H^2);
Rm=Rmax-Rmin+c*Tr/2;
M=ceil(2*Rm/c/dt);
M=2^nextpow2(M);
t=linspace(2*Rmin/c-Tr/2,2*Rmax/c+Tr/2,M);

r=c*t/2;
dt=(2*Rmax/c+Tr-2*Rmin/c)/M;
Fs=1/dt;
Ntarget=3;
Ptarget=[0+0.25,R0-1,1;50+1,R0+100,0.8;100,R0+100,0.8];





s_echo=zeros(N,M);
for k=1:3
    R=sqrt(Ptarget(k,2)^2+(X-Ptarget(k,1)).^2);
    delay=2*R/c;
    Delay=ones(N,1)*t-delay'*ones(1,M);
    Phase=1j*pi*Kr*Delay.^2-1j*4*pi*fc*(R'*ones(1,M))/c;
    s_echo=s_echo+Ptarget(k,3)*rectpuls(Delay/Tr).*...
    rectpuls((X-Ptarget(k,1))'*ones(1,M)/Ls).*exp(Phase);
end
figure;
imagesc(real(s_echo));
title('回波数据');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');




fr=(-M/2:M/2-1)*Fs/M;
mf_f=exp(1j*pi*fr.^2/Kr);
s_echo=fftshift(fft(fftshift(s_echo),[],2));
src_f=s_echo.*(ones(N,1)*(mf_f));




Src=fftshift(fft(fftshift(src_f),[],1));




Fa=(-N/2:N/2-1)*PRF/N;
D_Fa_Vr=sqrt(1-lambda^2*Fa.^2/4/Vr^2);
Ksrc=(2*Vr^2*fc^3*D_Fa_Vr.^3/c./Fa).'*(1./r);
for k=1:N
    Hsrc=exp(-1j*pi.*fr.^2./Ksrc(k,:));
    Src(k,:)=Src(k,:).*Hsrc;
end
Src=fftshift(ifft(fftshift(Src),[],2));
figure;
imagesc(abs(Src));
title('二次距离压缩后信号');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');




fa=(-N/2:N/2-1)*PRF/N;
RCM=lambda^2*R0*fa.^2/8/Vr^2;
RCM_point=RCM*2*Fs/c;
N_interp=8;
N_add=N_interp*ceil(max(RCM_point));
Src_RCMC=zeros(N,M);
for k1=1:N
    n_rcm=round(((1:M)+RCM_point(k1)-1)*N_interp+1);
    Src_interp=zeros(1,M*N_interp+N_add);
    Src_interp(1:M*N_interp)=interp(Src(k1,:),N_interp);
    Src_RCMC(k1,:)=Src_interp(n_rcm);
end
figure;
imagesc(abs(Src_RCMC));
title('RCMC后信号');
xlabel('距离时域（采样点）');
ylabel('方位频域（采样点）');

data1=fftshift(fft(fftshift(Src_RCMC(:,510))));
figure;
imagesc(abs(data1));




Src_mfa=zeros(N,M);
for k=1:M
    mfa_f=exp(1j*pi/(-2*Vr^2/lambda/r(k))*fa.^2);
    Src_mfa(:,k)=Src_RCMC(:,k).*mfa_f.';
end




sac=fftshift(ifft(fftshift(Src_mfa),[],1));




figure;
imagesc(abs(sac));
title('点目标成像');
xlabel('距离时域（采样点）');
ylabel('方位时域（采样点）');




[location_a,location_r]=find(abs(sac)==max(max(abs(sac))));
Ns=8;
sac_interp_r=interp(sac(:,location_r).',Ns);
sac_interp_r_abs=abs(sac_interp_r);
sac_interp_r_absmax=sac_interp_r_abs/max(sac_interp_r_abs);
sac_interp_r_log=20*log10(sac_interp_r_absmax);
rr=find(sac_interp_r_log==max(sac_interp_r_log));

sac_interp_a=interp(sac(location_a,:),Ns);
sac_interp_a_abs=abs(sac_interp_a);
sac_interp_a_absmax=sac_interp_a_abs/max(sac_interp_a_abs);
sac_interp_a_log=20*log10(sac_interp_a_absmax);
aa=find(sac_interp_a_log==max(sac_interp_a_log));

pslr_r=pslrfunc(sac_interp_r_abs)
islr_r=islrfunc(sac_interp_r_abs)
pslr_a=pslrfunc(sac_interp_a_abs)
islr_a=islrfunc(sac_interp_a_abs)

figure(7);
subplot(221);
plot(sac_interp_r_log);
axis([rr-150,rr+150,-30,0]);
ylabel('幅度dB');title('(a)距离剖面图幅度');

subplot(222);
plot(sac_interp_a_log);
axis([aa-150,aa+150,-30,0]);
ylabel('幅度dB');title('(b)方位剖面图幅度');

subplot(223);
plot(angle(sac_interp_r));
axis([rr-150,rr+150,-4,4]);
xlabel('距离向（采样点）');ylabel('相位 度');
title('(c)距离剖面图相位');

subplot(224);
plot(angle(sac_interp_a));
axis([aa-150,aa+150,-4,4]);
xlabel('方位向（采样点）');ylabel('相位 度');
title('(d)方位剖面图相位');




NN=ceil(c/2/Br/(Vr/Ba)/2)*2;

Ran_coordinate=1025;
Azi_coordinate=509;
chooseNum=20;
a=sac(Ran_coordinate-ceil(chooseNum/NN/2)*2:Ran_coordinate+ceil(chooseNum/NN/2)*2,...
Azi_coordinate-chooseNum:Azi_coordinate+chooseNum-1);
Interp_num=8;
[range_num,azimuth_num]=size(a);
a_temp=zeros(range_num*(2*Interp_num+1),azimuth_num);
a_final=zeros(range_num*(2*Interp_num+1),azimuth_num*(2*Interp_num+1));
for n=1:azimuth_num

    data_f=fft(a(:,n));
    data=[zeros(1,range_num*Interp_num),data_f.',zeros(1,range_num*Interp_num)];
    a_temp(:,n)=ifft(data).';
end
for n=1:range_num*(2*Interp_num+1)
    data_f=fft(a_temp(n,:));
    data=[zeros(1,azimuth_num*Interp_num),data_f,zeros(1,azimuth_num*Interp_num)];
    a_final(n,:)=ifft(data);
end

figure;imagesc(abs(a_final))
figure;contour(abs(a_final),10);




x=sac_interp_r_abs;
soo=x;
s_number=8*8;
T=dt*8;
N_buling=s_number;
soo_abs=abs(soo);
[C,I]=max(soo_abs);
y=soo_abs.^2;
x1=0;
while(soo_abs(I-x1-1)-soo_abs(I-x1))<0
    M1=x1;
    x1=x1+1;
end
x2=0;
while(soo_abs(I+x2+1)-soo_abs(I+x2))<0
    M2=x2;
    x2=x2+1;
end
P1=I-1-M1;
P2=I+1+M2;

M=(10^(-3/20))*C;

z1=abs(soo_abs(P1)-M);
x1=1;
z1_x1=0;
for k1=P1:I
    cha1=abs(soo_abs(P1+x1)-M);
    if cha1<z1
        z1=cha1;
        z1_x1=x1;
    end
    x1=x1+1;
end
z2=abs(soo_abs(I)-M);
x2=1;
z2_x2=0;
for k2=I:P2
    cha2=abs(soo_abs(I+x2)-M);
    if cha2<z2
        z2=cha2;
        z2_x2=x2;
    end
    x2=x2+1;
end
Th_x1=P1+z1_x1;
Th_x2=I+z2_x2;



if soo_abs(Th_x1)-M<0
    x0_linear=Th_x1;
    x1_linear=Th_x1+1;
else
    x0_linear=Th_x1-1;
    x1_linear=Th_x1;
end
Th_x1_real=...
(M-soo_abs(x1_linear))/(soo_abs(x0_linear)-soo_abs(x1_linear))*x0_linear...
+(M-soo_abs(x0_linear))/(soo_abs(x1_linear)-soo_abs(x0_linear))*x1_linear;

if soo_abs(Th_x2)-M>0
    x0_linear=Th_x2;
    x1_linear=Th_x2+1;
else
    x0_linear=Th_x2-1;
    x1_linear=Th_x2;
end
Th_x2_real=...
(M-soo_abs(x1_linear))/(soo_abs(x0_linear)-soo_abs(x1_linear))*x0_linear...
+(M-soo_abs(x0_linear))/(soo_abs(x1_linear)-soo_abs(x0_linear))*x1_linear;

width=Th_x2_real-Th_x1_real;
IRW=(T/N_buling)*width*c/2;


IRW_juli=IRW



x=sac_interp_a_abs;
soo=x;
s_number=8*8;
T=(dt)*8;
N_buling=s_number;
soo_abs=abs(soo);
[C,I]=max(soo_abs);
y=soo_abs.^2;
x1=0;
while(soo_abs(I-x1-1)-soo_abs(I-x1))<0
    M1=x1;
    x1=x1+1;
end
x2=0;
while(soo_abs(I+x2+1)-soo_abs(I+x2))<0
    M2=x2;
    x2=x2+1;
end

P1=I-1-M1;
P2=I+1+M2;
M=(10^(-3/20))*C;

z1=abs(soo_abs(P1)-M);
x1=1;
z1_x1=0;
for k1=P1:I
    cha1=abs(soo_abs(P1+x1)-M);
    if cha1<z1
        z1=cha1;
        z1_x1=x1;
    end
    x1=x1+1;
end
z2=abs(soo_abs(I)-M);
x2=1;
z2_x2=0;
for k2=I:P2
    cha2=abs(soo_abs(I+x2)-M);
    if cha2<z2
        z2=cha2;
        z2_x2=x2;
    end
    x2=x2+1;
end
Th_x1=P1+z1_x1;
Th_x2=I+z2_x2;



if soo_abs(Th_x1)-M<0
    x0_linear=Th_x1;
    x1_linear=Th_x1+1;
else
    x0_linear=Th_x1-1;
    x1_linear=Th_x1;
end
Th_x1_real=...
(M-soo_abs(x1_linear))/(soo_abs(x0_linear)-soo_abs(x1_linear))*x0_linear...
+(M-soo_abs(x0_linear))/(soo_abs(x1_linear)-soo_abs(x0_linear))*x1_linear;

if soo_abs(Th_x2)-M>0
    x0_linear=Th_x2;
    x1_linear=Th_x2+1;
else
    x0_linear=Th_x2-1;
    x1_linear=Th_x2;
end
Th_x2_real=...
(M-soo_abs(x1_linear))/(soo_abs(x0_linear)-soo_abs(x1_linear))*x0_linear...
+(M-soo_abs(x0_linear))/(soo_abs(x1_linear)-soo_abs(x0_linear))*x1_linear;

width=Th_x2_real-Th_x1_real;
IRW=(T/N_buling)*width*c/2;


IRW_fangwei=IRW


