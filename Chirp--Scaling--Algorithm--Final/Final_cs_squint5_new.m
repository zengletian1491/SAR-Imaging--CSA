%%  Interpretations
%  不考虑高度
%  斜视 5 度,斜平面成像chirp scaling算法
%%  radar parameters
clear all;
close all;
clc;
c = 3e8; 
lambda = 0.03125;
fc = c/lambda;
Tp = 10e-6; 
B = 180e6; 
gama = B/Tp; 
Fs = 200e6;  Fs = 400e6;
H = 0;
Rs = 12.8e3; 
squint_angle = 5;
squint = squint_angle/180*pi ; 
Rb = Rs*cos(squint) ;
Da = 0.4;
L = Rs*tan(lambda/Da/2)*2/cos(squint);
tp_nrn = ceil(Tp*Fs/2)*2;
v = 110;
PRF = 5000/3;
Wg = 1.5e3;
Wa = 4e2;
Pr = c/2/B;
Wr = Wg/cos(squint);
Ts = 2*Wr/c ; 

%%  range parameter
DeltaR = c/2/Fs;
nrn  = ceil((Tp + Ts)*Fs/2)*2;
nrn = 2^ceil(log2(nrn));
Tstart = 2*Rs/c-nrn/2/Fs;
Tnrn = 1/Fs;
Tend = 2*Rs/c+(nrn/2-1)/Fs;
tnrn = [Tstart:Tnrn:Tend].';

%%  azimuth parameter
Rs*lambda/Da/(v/PRF)
DeltaX = v/PRF ;
nan = 8192;  nan = 4096;
tnan = [-nan/2:nan/2-1]/PRF;
delta_theta = 2*atan((nan*v/PRF)/Rs); 
Pa = lambda/2/delta_theta;

%%  target parameter
point_num = 9; 
Distance_a = 100;           
Distance_g = 300;     
squint_g = 0;
point = [...
         -Distance_a*cos(squint_g)-Distance_g*sin(squint_g),-Distance_a*sin(squint_g)+Distance_g*cos(squint_g)   ,   0;
                                -Distance_g*sin(squint_g),                          Distance_g*cos(squint_g) ,   0;
          Distance_a*cos(squint_g)-Distance_g*sin(squint_g), Distance_a*sin(squint_g)+Distance_g*cos(squint_g)   ,   0;
       -Distance_a*cos(squint_g)                         ,-Distance_a*sin(squint_g)                          ,   0;
                                                      0,                                                0,   0;
        Distance_a*cos(squint_g)                         , Distance_a*sin(squint_g)                          ,   0;                       
         -Distance_a*cos(squint_g)+Distance_g*sin(squint_g),-Distance_a*sin(squint_g)-Distance_g*cos(squint_g)   ,   0;
                                 Distance_g*sin(squint_g),                         -Distance_g*cos(squint_g) ,   0;
          Distance_a*cos(squint_g)+Distance_g*sin(squint_g), Distance_a*sin(squint_g)-Distance_g*cos(squint_g)   ,   0];
point = point + [zeros(size(point,1),1),  ones(size(point,1),1)*Rs,  zeros(size(point,1),1)];    
%%  platform parameter
X_ac = 0; 
Y_ac = 0; 
Z_ac = H ; 

X_a = ones(1,nan)*X_ac + v*tnan; 
Y_a = ones(1,nan)*Y_ac;
Z_a = ones(1,nan)*Z_ac;

%%  echo generation
temp = zeros(nrn,1); 
x0 = zeros(nrn,nan) ;
for k = 1 : nan
    if rem(k,400) == 0
        k
    end
    for m = 1 : 1% point_num 
        R_m(k) = sqrt((X_a(k) - point(m,1)).^2 + (Y_a(k) - point(m,2)).^2 + (Z_a(k) - point(m,3)).^2);
        alpha = asin((point(m,1) - X_a(k))./R_m(k));
%         if (abs(alpha - squint) <= lambda/2/Da) 
            win1 = (abs(tnrn - 2*R_m(k)/c)<=Tp/2);
            temp = win1.*exp(j*pi*gama*(tnrn - 2*R_m(k)/c).^2).*exp(-j*4*pi*(R_m(k))/lambda);
            x0(:,k) = x0(:,k) + temp;
%         end
    end
end
% figure,imagesc(abs(x0))
% figure,imagesc(abs(Range_compress(x0,Fs,Tp,gama)))


%%  Chirp Scaling Algorithm
Rb = Rs*cos(squint);
fdc = 2*v*sin(squint)/lambda;
fa = [-nan/2:nan/2-1]/nan*PRF;
fd = fa + fdc;
faM = 2*v*cos(squint)/lambda;
afa = 1./sqrt(1-(fd./faM).^2)-1;
theta = asin(fd/faM);
gamafa = 1./(1/gama-Rb*2*lambda*sin(theta).^2./c.^2./cos(theta).^3);
sref1 = exp(-j*2*pi*fdc*tnan);
x = zeros(nrn,nan);
for n = 1:nrn
    x(n,:) = fftshift(fft(fftshift(x0(n,:).*sref1)));
    H1CS = exp(j*pi*gamafa.*afa.*(tnrn(n)-2*Rs*(1+afa)/c).^2);
    x(n,:) = x(n,:).*H1CS;
end

fr = [-nrn/2:nrn/2-1].'/nrn*Fs;
for m = 1:nan
    x(:,m) = fftshift(fft(fftshift(x(:,m))));
    H2 = exp(j*pi./gamafa(m)./(1+afa(m)).*fr.^2).*exp(j*4*pi*Rs*afa(m)/c.*fr);
    x(:,m) = fftshift(ifft(fftshift(x(:,m).*H2)));
end

rs = [-nrn/2:nrn/2-1].'*(c/2/Fs);
for n = 1:nrn
    H3 = exp(j*2*pi*(Rb+rs(n,1))/v*sqrt(faM.^2-fd.^2))...
        .*exp(-j*4*pi/c^2*gamafa.*afa.*(1+afa).*(rs(n).^2));
    x(n,:) = fftshift(ifft(fftshift(x(n,:).*conj(sref1).*H3)));
end
figure,imagesc(abs(x))
x_f0 = fftshift(fft(fftshift(x,1),[],1),1);
x_f0 = fftshift(fft(fftshift(x_f0,2),[],2),2);  %  二维频谱不混叠
figure,imagesc(abs(x_f0))  


%%
Image_f = x;
Interp_num = 8;
pos_row = 4897;   (pos_row - nrn/2 - 1)*c/2/Fs
pos_col = 3122;   (pos_col - nan/2 - 1)*v/PRF
num_row = 64;
num_col = 256;
data_0 = Interpolate2D(Image_f(pos_row-num_row/2+1:pos_row+num_row/2,pos_col-num_col/2+1:pos_col+num_col/2),Interp_num);
figure,contour(abs(data_0),40)
% fig_mesh(Image_f(pos_row+(-8*2:8*2-1),pos_col+(-8*2:8*2-1)),40);

if 1
x_r_ori = Image_f(pos_row-num_row/2+1:pos_row+num_row/2,pos_col);
x_r_irf = zeros(num_row*Interp_num,1);
x_r_irf(num_row*Interp_num/2-num_row/2+1:num_row*Interp_num/2+num_row/2,1) = fftshift(fft(fftshift(x_r_ori)));
x_r_irf = fftshift(ifft(fftshift(x_r_irf)));
figure,plot(20*log10(abs(x_r_irf)./max(abs(x_r_irf))))
PSLR(x_r_irf)                        %  -13.278801650672033
ISLR(x_r_irf)                        %  -10.025710459125843
disp('距离向分辨率理论值： （m）')
c/2/B                                %  0.833333333333333
disp('距离向分辨率实际值： （m）')
IRW(x_r_irf)/Interp_num*(c/2/Fs)     %  0.703125000000000

x_a_ori = Image_f(pos_row,pos_col-num_col/2+1:pos_col+num_col/2);
x_a_irf = zeros(1,num_col*Interp_num);
x_a_irf(1,num_col*Interp_num/2-num_col/2+1:num_col*Interp_num/2+num_col/2) = fftshift(fft(fftshift(x_a_ori)));
x_a_irf = fftshift(ifft(fftshift(x_a_irf)));
figure,plot(20*log10(abs(x_a_irf)./max(abs(x_a_irf))))
PSLR(x_a_irf)                        %  -13.255949728326140
ISLR(x_a_irf)                        %  -10.242338896862037
disp('方位向分辨率理论值： （m）')
lambda/2/(v*nan/PRF/Rs)              %  0.739820075757576
disp('方位向分辨率实际值： （m）')
IRW(x_a_irf)/Interp_num*(v/PRF)      %  0.668250000000000
end