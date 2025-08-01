%  SSSP Homework - 12/06/2025
%  Matteo Di Giovanni & Alessandro Mancuso

clear; close all; clc

%% Simscape File for Ground Truth

load('ssc_output.mat')

% To evaluate time performances 
tic

%% Sampling Frequency
fs = 192e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 2;  % [seconds]

%% Input Signal

% Time Axis
t = 0:Ts:stop_time;
t = t(1:end-1);

% Signal Amplitude
A = 10;
vin = A * sweeptone(1.99, 0.01, fs, 'SweepFrequencyRange', [500 20000]);

%% Circuit Parameters
% Piezoelectric MEMS loudspeaker in Free-field conditions

% Transformer Ratio
alpha = 3.7e-4;
Seff = 2.0e-5;

% Electrical Domain
Re = 4;
Cp = 2.4e-8;

% Mechanical Domain
Rm = 9.7e-3;
Mm = 1e-6;
Cm = 2.2e-3;

% Acoustic Domain
Cbc = 3.6e-13; 
Ltube1 = 1e2;
Ltube2 = 1e2;
Ctube = 6.5e-13;
Rac = 5e6;

%% Removing Ideal Transformers (Mechanical Domain only)

gamma1 = alpha^-1;
gamma2 = Seff;

% Resistive Elements
R1 = Re/gamma1^2;
R2 = Rm;
R3 = Rac*gamma2^2;

% Dynamic Elements
% Capacitors
C1 = Cp*gamma1^2;
C2 = Cm;
C3 = Cbc/gamma2^2;
C4 = Ctube/gamma2^2;

% Inductors
L1 = Mm;
L2 = Ltube1*gamma2^2;
L3 = Ltube2*gamma2^2;

% Source
Fin = alpha * vin;

%% Setting of Free Parameters (Adaptation Conditions)
Z2 = R1;
Z6 = Ts/(2*C1);
Z9 = R2;
Z10 = (2*L1)/Ts;
Z11 = Ts/(2*C2);
Z14 = Ts/(2*C3);
Z17 = (2*L2)/Ts;
Z20 = Ts/(2*C4);
Z22 = R3; % Output
Z23 = (2*L3)/Ts;

% Make Ports of Adaptors Reflection-Free 
Z21 = Z22 + Z23;
Z19 = Z21;
Z18 = (Z19*Z20)/(Z19+Z20);
Z16 = Z18;
Z15 = Z16 + Z17;
Z13 = Z15;
Z12 = (Z14*Z13)/(Z14+Z13);
Z8 = Z12;
Z7 = Z8 + Z9 + Z10 + Z11;
Z5 = Z7;
Z4 = (Z5*Z6)/(Z5+Z6);
Z1 = Z4;
Z3 = Z1 + Z2; % Facing ideal voltage source, can't be adapted


%% Computing Scattering Matrices
% Series Junctions
gammaSer1 = Z2/(Z1+Z2);
Sser1 = [gammaSer1, (gammaSer1-1), (gammaSer1-1);
        -gammaSer1, (1-gammaSer1), -gammaSer1   ;
          -1      , -1           ,  0           ;];

% 5-Port Series junction obtained through B matrix method
Z_vec = [Z7, Z8, Z9, Z10, Z11];
Z = diag(Z_vec);

% Since it's a series junction all elements share the same current
B = [1, 1, 1, 1, 1];

Sser2 = eye(5) - 2*Z*B'*inv(B*Z*B')*B;

gammaSer3 = Z16/(Z16+Z17);
Sser3 = [ 0       , -1           , -1         ;
       -gammaSer3 , (1-gammaSer3), -gammaSer3 ;
     (gammaSer3-1), (gammaSer3-1), gammaSer3  ;];

gammaSer4 = Z22/(Z22+Z23);
Sser4 = [ 0       , -1           , -1         ;
       -gammaSer4 , (1-gammaSer4), -gammaSer4 ;
     (gammaSer4-1), (gammaSer4-1), gammaSer4  ;];

% Parallel Junctions
gammaPar1 = Z5/(Z5+Z6);
Spar1 = [ 0       , (1-gammaPar1), gammaPar1  ;
          1       ,   -gammaPar1 , gammaPar1  ;
          1       , (1-gammaPar1),(gammaPar1-1);];

gammaPar2 = Z13/(Z13+Z14);
Spar2 = [ 0       , (1-gammaPar2), gammaPar2  ;
          1       ,   -gammaPar2 , gammaPar2  ;
          1       , (1-gammaPar2),(gammaPar2-1);];

gammaPar3 = Z19/(Z19+Z20);
Spar3 = [ 0       , (1-gammaPar3), gammaPar3  ;
          1       ,   -gammaPar3 , gammaPar3  ;
          1       , (1-gammaPar3),(gammaPar3-1);];

%% Initialization of Waves

% Manual initialization of all variables a1–a23 and b1–b23 to zero.
% An alternative (and more scalable) approach would have been: 
% a = zeros(1, 23); b = zeros(1, 23); 
% followed by index-based access (e.g., a(i), b(i)) in the for loop.
a1 = 0; a2 = 0; a3 = 0; b1 = 0; b2 = 0; b3 = 0;
a4 = 0; a5 = 0; a6 = 0; b4 = 0; b5 = 0; b6 = 0;
a7 = 0; a8 = 0; a9 = 0; b7 = 0; b8 = 0; b9 = 0;
a10 = 0; a11 = 0; a12 = 0; b10 = 0; b11 = 0; b12 = 0;
a13 = 0; a14 = 0; a15 = 0; b13 = 0; b14 = 0; b15 = 0;
a16 = 0; a17 = 0; a18 = 0; b16 = 0; b17 = 0; b18 = 0;
a19 = 0; a20 = 0; a21 = 0; b19 = 0; b20 = 0; b21 = 0;
a22 = 0; b22 = 0; a23 = 0; b23 = 0;


%% Initialization of Output Signals
Fout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(Fin)
    
    % Managing Dynamic Elements (Inductors, Capacitors);
    a6 = b6; % C1
    a10 = -b10; % L1
    a11 = b11; % C2
    a14 = b14; % C3
    a17 = -b17; % L2
    a20 = b20; % C4
    a23 = -b23; % L3

    % Forward Scan
    b21 = Sser4(1,:)*[0;a22;a23];
    a19 = b21;

    b18 = Spar3(1,:)*[0;a19;a20];
    a16 = b18;

    b15 = Sser3(1,:)*[0;a16;a17];
    a13 = b15;

    b12 = Spar2(1,:)*[0;a13;a14];
    a8 = b12;

    b7 = Sser2(1,:)*[0;a8;a9;a10;a11];
    a5 = b7;

    b4 = Spar1(1,:)*[0;a5;a6];
    a1 = b4;
    
    b3 = Sser1(3,:)*[a1;a2;0];

    % Local Root Scattering
    a3 = 2*Fin(n) - b3;

    % Backward Scan
    b1 = Sser1(1,:)*[a1;a2;a3];
    a4 = b1;

    b2 = Sser1(2,:)*[a1;a2;a3];

    b5 = Spar1(2,:)*[a4;a5;a6];
    a7 = b5;

    b6 = Spar1(3,:)*[a4;a5;a6];

    b8 = Sser2(2,:)*[a7;a8;a9;a10;a11];
    a12 = b8;

    b9 = Sser2(3,:)*[a7;a8;a9;a10;a11];

    b10 = Sser2(4,:)*[a7;a8;a9;a10;a11];

    b11 = Sser2(5,:)*[a7;a8;a9;a10;a11];
    
    b13 = Spar2(2,:)*[a12;a13;a14];
    a15 = b13;

    b14 = Spar2(3,:)*[a12;a13;a14];

    b16 = Sser3(2,:)*[a15;a16;a17];
    a18 = b16;

    b17 = Sser3(3,:)*[a15;a16;a17];

    b19 = Spar3(2,:)*[a18;a19;a20];
    a21 = b19;

    b20 = Spar3(3,:)*[a18;a19;a20];

    b23 = Sser4(3,:)*[a21;a22;a23];

    b22 = Sser4(2,:)*[a21;a22;a23];

    % Read Output
    Fout(n) = (a22+b22)/2;

end

%% Output Plots

% Computing acoustic pressure
pout = Fout ./ Seff;

% Time Domain Plots
figure
set(gcf, 'Color', 'w');
plot(t, pout, 'b', 'Linewidth', 2);
hold on
plot(gt(1, :), gt(2, :), 'r--', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$p_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Pressure - Time Domain','Fontsize',18,'interpreter','latex');

% Frequency domain Plots
nfft = 2^20;
res = fs/nfft;
f = (0:nfft/2-1) * res;

ir = impzest(vin/A, pout');
ir_gt = impzest(vin/A, gt(2, 1:end-1)');

tf = fft(ir, nfft);
tf_gt = fft(ir_gt, nfft);

abs_tf = abs(tf(1:nfft/2));
abs_tf_gt = abs(tf_gt(1:nfft/2));

figure
set(gcf, 'Color', 'w');
semilogx(f, 20*log10(abs_tf/2e-5), 'b', 'Linewidth', 2);
hold on
semilogx(f, 20*log10(abs_tf_gt/2e-5), 'r--', 'Linewidth', 2);
grid on;
xlim([500, 20000])
xlabel('Frequency [Hz]','Fontsize',16,'interpreter','latex');
ylabel('$\mathrm{SPL}\,[\mathrm{dB}_\mathrm{SPL}]$','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Sound Pressure Level - Frequency Domain','Fontsize',16,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, pout - gt(2, 1:end-1), 'k', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$\mathcal{E}_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',16,'interpreter','latex');

%% Compute Mean Squared Error (MSE)

mse = mean((pout - gt(2, 1:end-1)).^2);
disp('MSE = ')
disp(mse)

% To evaluate time performances 
toc
