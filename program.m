%---------- Exercise 1 ------------


%An intro to signal analysis 

%---- 1Α ----
Fs = 4000; Channels = 1; bits = 16;
r = audiorecorder(Fs, bits, Channels);
duration = 3; disp('Recording Started');
recordblocking(r, duration);
disp('Recording Stopped');


%---- 1Β ----
X = getaudiodata(r);
filename = 'myname.wav';
audiowrite(filename, X, Fs)
%sound(X,Fs,bits);
t = 0:1/Fs:(length(X)-1)/Fs;
figure(1);
plot(t, X, 'LineWidth',1.5);
start_seg = 8000*0.45;
end_seg = 8000*0.55; 
window = X(start_seg:end_seg);
figure(2);
plot(linspace(0, 0.05, length(window)), window);
sound(window, Fs, bits);

%---- 1Γ ----
s = normalize(X, 'range', [-1, 1]);
figure(3);
plot(t, s);

n = linspace(0,200,200);
M = 200;
N = 100;
w = 0.54 - 0.46*cos(2*pi*n/M);
u = conv(s.*s, w);
s = normalize(s, 'range', [min(u), max(u)]);
figure(4);
plot(t,s);
hold on;
plot(linspace(0, 3, length(u)), u);

%---- 1Δ ----
fw = fft(window, 1024);
f = linspace(0, 8000, length(fw));
figure(6);
plot(f, abs(fw));
figure(7);
plot(f, 20*log10(abs(fw)));


%---------- Exercise 2 ------------

%Using the Fourier Transform to reconstruct the outline of an image

leaf = imread('leaf3.png');
[row, col]=find(leaf);
figure(200);
imshow(leaf);
boundaries = bwtraceboundary(leaf,[row(1),col(1)],'N');
x_cord = boundaries(:,1);
y_cord=boundaries(:,2);
figure(1);
plot(x_cord);
hold on; 
plot(y_cord);
hold off;

z = complex(x_cord,y_cord);

Z = fft(z);
figure(201);
plot(abs(Z), 'LineWidth', 2.0);
N = length(z);
z_10 = zeros(N:1);
z_50 = zeros(N:1);
z_200 = zeros(N:1);
for n = 1:N
    z_10(n) = z_m(N, 10, n, Z);
end
x_10 = abs(round(real(z_10)));
y_10 = abs(round(imag(z_10)));
l_10 = l_m(x_10, y_10);
figure(2);
imshow(l_10);

for n = 1:N
    z_50(n) = z_m(N, 50, n, Z);
end
x_50 = abs(round(real(z_50)));
y_50 = abs(round(imag(z_50)));
l_50 = l_m(x_50, y_50);
figure(3);
imshow(l_50);

for n = 1:N
    z_200(n) = z_m(N, 200, n, Z);
end
x_200 = abs(round(real(z_200)));
y_200 = abs(round(imag(z_200)));
l_200 = l_m(x_200, y_200);
figure(4);
imshow(l_200);

z_10_i = zeros(N:1);
z_50_i = zeros(N:1);
z_200_i = zeros(N:1);
for n = 1:N
    z_10_i(n) = z_m(N, 10/2, n, Z)+z_m_s(N, 10, n, Z);
end
x_10_i = abs(round(real(z_10_i)));
y_10_i = abs(round(imag(z_10_i)));
l_10_i = l_m(x_10_i, y_10_i);
figure(5);
imshow(l_10_i);

for n = 1:N
    z_50_i(n) = z_m(N, 50/2, n, Z)+z_m_s(N, 50, n, Z);
end
x_50_i = abs(round(real(z_50_i)));
y_50_i = abs(round(imag(z_50_i)));
l_50_i = l_m(x_50_i, y_50_i);
figure(6);
imshow(l_50_i);

for n = 1:N
    z_200_i(n) = z_m(N, 200/2, n-1, Z)+z_m_s(N, 200, n-1, Z);
end
x_200_i = abs(round(real(z_200_i)));
y_200_i = abs(round(imag(z_200_i)));
l_200_i = l_m(x_200_i, y_200_i);
figure(7);
imshow(l_200_i);

figure(8);
imshow(leaf);


%-------------------------------------------

%Testing it again for a different image


heart = imread('41WWBQjvPqL._SR600,315_PIWhiteStrip,BottomLeft,0,35_SCLZZZZZZZ_FMpng_BG255,255,255.png');
heart = im2bw(heart, 0.5);
figure(9);
imshow(heart);

[row, col] = find(heart);
boundary = bwtraceboundary(heart, [row(1), col(1)], 'S');
x_cord = boundary(:,1);
y_cord = boundary(:,2);
z = complex(x_cord,y_cord);
Z = fft(z); 
N = length(z);
z_200 = zeros(N,1);
for n=1:N
    z_200(n)=z_m(N, 2000/2, n-1, Z) + z_m_s(N, 2000, n-1, Z);
end
x_200 = abs(round(real(z_200)));
y_200 = abs(round(imag(z_200)));
l_200 = l_m(x_200, y_200);
figure(10);
imshow(l_200);





%---------------------------------------------

function z_m = z_m(N, M, n, z)
    z_m = 0;
    for k=1:M+1
        z_m = z_m + (1/N)*(z(k)*exp(1i*2*pi*(k-1)*n/N));
    end
end

function z_m_s = z_m_s(N, M, n, z)
    z_m_s = 0;
    K = N - (M/2);
    for k=K:N-1
        z_m_s = z_m_s + (1/N)*(z(k+1)*exp(1i*2*pi*k*n/N));
    end
end

function l_m = l_m(x, y)
    for i=1:length(x)
        l_m(x(i)+1, y(i)+1)=1;
    end
end

%---------- Exercise 3 ------------

%Bandpass filter creation from scratch using poles and zeros

%--------- 3.1 ---------

poles = [complex(0.51, 0.68); complex(0.51,-0.68)];
zeros = [0.8; -1];
figure(1);
zplane(zeros,poles);
[b,a] = zp2tf(zeros,poles,0.15);

figure(3);
freqz(b,a,2001);

figure(4);
impz(b,a,30);

figure(5);
stepz(b,a,30);
close all;

poles2 = [complex(0.57, 0.76); complex(0.57,-0.76)];
poles3 = [complex(0.6, 0.8); complex(0.6,-0.8)];
poles4 = [complex(0.63, 0.84); complex(0.63,-0.84)];

[b2,a2] = zp2tf(zeros,poles2,0.15);
[b3,a3] = zp2tf(zeros,poles3,0.15);
[b4,a4] = zp2tf(zeros,poles4,0.15);

figure(201);
zplane(zeros,poles2);

figure(8);
stepz(b2,a2,70);

figure(301);
freqz(b2,a2,2001);

figure(202);
zplane(zeros,poles3); 

figure(9);
stepz(b3,a3,70);

figure(302);
freqz(b3,a3,2001);

figure(203);
zplane(zeros,poles4);

figure(10);

stepz(b4,a4,70);



n = linspace(-pi, pi, 101);
x = gensig('sin',6.6,100,1) + gensig('sin',2.85,100,1);
y=filter(b,a,x);
figure(11);
plot(n, x);
hold on;
plot(n, y);


poles5 = [complex(0.68, 0.51); complex(0.68, -0.51); complex(0.68, 0.51); complex(0.68, -0.51)];
figure(12);
zplane(zeros,poles5);
[b5,a5] = zp2tf(zeros,poles5,0.15);

figure(13);
freqz(b5,a5,2001);

y1=filter(b5,a5,x);
figure(14);
plot(n,x);
hold on;
plot(n,y1);
close all;
clear all; clc;

%----- 3.2 -------

[viola, Fs] = audioread('viola_series.wav');
t = linspace(0, 400000/Fs, 400000);
figure(15);
plot(t, viola);
sound(viola, Fs);

dF = Fs/400000;
f = -Fs/2:dF:Fs/2-dF;

fourier = fft(viola);

figure(16);
plot(linspace(-Fs/2,Fs/2,length(viola)), abs(fftshift(fourier)));
xlim([100 1500])

[viola_note,Fs2] = audioread('viola_note.wav');

figure(401);
plot(viola_note);

fourier2 = fft(viola_note);

figure(18);
plot(linspace(-Fs/2,Fs/2,length(viola_note)), abs(fftshift(fourier2)));
xlim([100, 1500])

z = [0.9075; 0.9075i; -0.9075; -0.9075i];
p = [0.995+0.092*1i; 0.995+0.092*1i; 0.995-0.092*1i; 0.995-0.092*1i];
figure(20);
zplane(z, p);
[b, a] = zp2tf(z, p, 10^-7);
figure(21);
freqz(b,a);
[h, ~] = impz(b, a);
har2 = conv(h, viola_note);
f = linspace(-Fs/2,Fs/2,length(har2));
figure(19);
plot(f, abs(fftshift(fft(har2))));
xlim([0, 700])
figure(23);
plot(linspace(0, length(har2)/Fs,length(har2)), har2);
xlim([0, 0.4]);

close all;


u = 0.89;
z = [u; u*1i; -u; -u*1i];
x= 0.99;
y= 0.14;
p2 = [x+y*1i; x+y*1i; x-y*1i; x-y*1i;];
[b2, a2] = zp2tf(z, p2, 10^-7);
figure(24);
zplane(z, p2);
figure(25);
freqz(b2, a2);
[h, ~] = impz(b2, a2);
sig = conv(h, viola_note);
f = linspace(-Fs/2,Fs/2,length(sig));
figure(26);
plot(f, abs(fftshift(fft(sig))));
xlim([0,1000])
figure(27);
plot(linspace(0, length(sig)/Fs,length(sig)), sig);




%---------- Exercise 4 ------------

%Creating multiple bandpass filters from scratch in order
%to split apart a mixed signal of two signals. 

%------ 4 --------


[mix, Fs] = audioread('mixture.wav');
%sound(mix, Fs);


fmix = fft(mix);
Fs_n = Fs/(2*pi);
figure(1);
f =linspace(-Fs/2, Fs/2, length(mix));
f_n = linspace(-Fs_n/2, Fs_n/2, length(mix));
plot(f, abs(fftshift(fmix)));
xlim([0, 3000]);

fmix_n = normalize(abs(fmix), 'scale', [0, 2*pi]);

figure(2);
plot(f_n, abs(fftshift(fmix_n)));
xlim([0, 200]);

close all;



%{
figure(200);
plot(f, fftshift(fft(mix)));
xlim([0 3000])
%}

th1 = bpf(mix, 0.98, 0.18, 0.9, 10^-7);
th1(end+1:80000)=0;
harm12 = bpf(mix, 0.93, 0.33, 0.93, 10^-7);
harm12(end+1:80000)=0;
harm13 = bpf(mix, 0.88, 0.47, 0.9, 10^-7);
harm13(end+1:80000)=0;
harm14 = bpf(mix, 0.75, 0.65, 0.9, 7*10^-5);
harm14(end+1:80000)=0;
harm15 = bpf(mix, 0.65, 0.75, 0.9, 7*10^-5);
harm15(end+1:80000)=0;

th2 = bpf(mix, 0.98, 0.14, 0.9, 10^-7);
th2(end+1:96000)=0;
harm22 = bpf(mix, 0.945, 0.25, 0.9, 10^-7);
harm22(end+1:96000)=0;
harm23 = bpf(mix, 0.91, 0.39, 0.9, 10^-7);
harm23(end+1:96000)=0;
harm24 = bpf(mix, 0.8535, 0.5199, 0.9, 12*10^-7);
harm24(end+1:96000)=0;
harm25 = bpf(mix, 0.794, 0.5806, 0.9, 7*10^-5);
harm25(end+1:96000)=0;

signal1 = th1 + harm12 + harm13 + harm14 + harm15;
signal2 = th2 +harm22 + harm23 + harm24 + harm25;

fourier1 = fft(signal1);
fourier2 = fft(signal2);

[mix2, Fs] = audioread('mixture2.wav');
[reed2, Fsr2] = audioread('reed_acoustic_037-057-127.wav');

th3 = bpf(mix2, 0.99, 0.0795, 0.9, 10^-7);
th3(end+1:70000)=0;
harm32 = bpf(mix2, 0.945, 0.15, 0.9, 10^-7);
harm32(end+1:70000)=0;
harm33 = bpf(mix2, 0.9342, 0.273, 0.91, 10^-7);
harm33(end+1:70000)=0;
harm34 = bpf(mix2, 0.91, 0.345, 0.91, 10^-7);
harm34(end+1:70000)=0;
harm35 = bpf(mix2, 0.897, 0.42, 0.91, 10^-7);
harm35(end+1:70000)=0;

signal3 = th3 + harm32 + harm33 + harm34 + harm35;

figure(3);
plot(signal3);

figure(4);
plot(reed2);

%{
[flute, Fsf] = audioread('flute_acoustic_002-069-025.wav');
[reed, Fsr] = audioread('reed_acoustic_037-065-075.wav');

figure(3);
plot(f, abs(fftshift(fft(flute))));
xlim([0, 3000]);

figure(4);
plot(f, abs(fftshift(fft(reed))));
xlim([0, 3000]);

figure(201);
plot(linspace(-Fs/2, Fs/2, length(fourier2)), abs(fftshift(fourier2)));
xlim([0 3000])

figure(7);
plot(flute);

figure(8);
plot(signal1);

figure(9);
plot(reed);

figure(10);
plot(signal2);
%}

function bpf = bpf(mix, x, y, u, K)
    z = [u*1i; -u*1i; u; -u;];
    p = [x+y*1i; x+y*1i; x-y*1i; x-y*1i;];
    [b, a] = zp2tf(z, p, K);
    [h, ~] = impz(b, a);
    bpf = conv(h, mix);
end
