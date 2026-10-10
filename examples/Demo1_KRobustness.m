%----------------------------------------------------------------------------------------------------
% "If you want to find the secrets of the universe, think in terms of energy, frequency & vibration."
%                                                            — Nikola Tesla
%----------------------------------------------------------------------------------------------------
% *[Overview]
%  This script evaluates the robustness of the SRSA  with respect to the specified K_guess value. 
%  The same complex-valued time series is analyzed using several K_guess values, and the resulting
%  diagnostic cross-sections at Q_L = 1 and L = N are compared in Figure (2).
%
% *[Caption for Figures]
%  Figure (1): Singular-value spectrum of the Hankel matrix constructed from the simulated time series
%
%  Figure (2):
%      Panel (a): diagnostic cross-sections with K_guess = 2.
%      Panel (b): diagnostic cross-sections with K_guess = 7.
%      Panel (c): diagnostic cross-sections with K_guess = 10.
%      Panel (d): diagnostic cross-sections with K_guess = 50.
%      Panel (e): diagnostic cross-sections with K_guess = 100.
%      Panel (f): diagnostic cross-sections with K_guess = 200.
%
%  Figure (3):
%      Panel (a): Normalized bilateral Fourier amplitude spectrum.
%      Panel (b): Full 2-D SRSA coherence response computed using K_true = 4, corresponding to the
%                 known number of complex exponential components in the simulated signal.
%
% *[Computational Environment]
%  The calculations were performed using MATLAB R2024b Update 3 on a workstation equipped with an
%  Intel Core i9-14900K processor (24 cores and 32 threads), 128 GB of RAM, and an NVIDIA GeForce
%  RTX 4090 D GPU with 24 GB of memory. The process-based parallel pool used eight workers. The
%  present implementation does not invoke GPU computing; therefore, the GPU did not contribute to
%  the reported runtimes.
%
% *[Representative Runtime]
%  The total wall-clock time for this demo was less than 1 minute on the stated workstation,
%  including numerical computation, parallel-pool initialization, and figure generation.
%----------------------------------------------------------------------------------------------------
% Author       :   Zhifeng Chen
% E-mail       :   zfchen@whu.edu.cn
% Created      :   2026-01-20
% Last updated :   2026-06-14
% Last updated :   2026-10-10
%====================================================================================================

%% Complex-valued sequences are used to display bilateral positive- and negative-frequency spectra.
clc;
clear;
close all

rng(20260025,'twister');                 % Fixed random seed for reproducibility

DeltaT = 1;                              % Sampling interval (s)
N      = 1000;                           % Total number of samples
time   = (0:N-1)*DeltaT;                 % Time (s)

% Each component has unit amplitude.
% The real and imaginary parts of the noise both have a standard deviation of 4.
yn = exp(1i*2*pi* 0.02*time) + ...
     exp(1i*2*pi*-0.05*time) + ...
     exp(1i*2*pi* 0.10*time) + ...
     exp(1i*2*pi*-0.45*time) + ...
     4*randn(1,N) + 1i*4*randn(1,N);

%% Singular-value spectrum
P = fix(N/2);                            % Number of Hankel columns
Q = N-P+1;                               % Number of Hankel rows

H = zeros(Q,P);                          % Initialize the Hankel matrix

for ii = 0:Q-1
    H(ii+1,:) = yn(ii+1:ii+P);           % Construct the Hankel matrix
end

singularV  = svd(H,'econ');              % Singular values in descending order
singularSN = (1:numel(singularV)).';      % Singular-value indices

%% Figure (1): singular-value spectrum
fig1 = figure(1);
clf(fig1);

plot(singularSN,singularV,'d');

xlabel('Singular-Value Index','FontName','Times New Roman','FontSize',12);

ylabel('Singular Value','FontName','Times New Roman','FontSize',12);

%% Fourier spectrum
Nfft = 3*N;
win  = hann(N);

% Bilateral frequency vector
FFTf = (-Nfft/2:Nfft/2-1).'/(Nfft*DeltaT);

% The coherent gain of the Hann window is corrected using sum(win).
FFTp = abs(fftshift(fft(yn(:).*win,Nfft)))/sum(win);

% Min-max normalization for comparison with the normalized SRSA spectrum
NormFFTp = (FFTp-min(FFTp))/(max(FFTp)-min(FFTp));

%% SRSA spectra
tic;
R  = 1;                                 % Frequency-grid oversampling factor
f0 = (0:R*N-1)/(DeltaT*N*R);

% Keep the input frequency grid unchanged for all SRSA calculations.
SRSbins = f0-f0(fix(R*N/2)+1);

K_guess1 = 2;
[SRSf,LVec1,SRSp1] = DoSRSA(yn,DeltaT,SRSbins,K_guess1,'AlgorithmMode','svd');

K_guess2 = 7;
[~,LVec2,SRSp2] = DoSRSA(yn,DeltaT,SRSbins,K_guess2,'AlgorithmMode','svd');

K_guess3 = 10;
[~,LVec3,SRSp3] = DoSRSA(yn,DeltaT,SRSbins,K_guess3,'AlgorithmMode','svd');

K_guess4 = 50;
[~,LVec4,SRSp4] = DoSRSA(yn,DeltaT,SRSbins,K_guess4,'AlgorithmMode','svd');

K_guess5 = 100;
[~,LVec5,SRSp5] = DoSRSA(yn,DeltaT,SRSbins,K_guess5,'AlgorithmMode','svd');

K_guess6 = 200;
[~,LVec6,SRSp6] = DoSRSA(yn,DeltaT,SRSbins,K_guess6,'AlgorithmMode','svd');

K_true = 4;
[~,Q_LVecTrue,SRSpTrue] = DoSRSA(yn,DeltaT,SRSbins,K_true,'AlgorithmMode','svd');
toc;
%%
fig2 = figure(2);clf(fig2);

set(fig2, ...
    'Units','centimeters', ...
    'Position',[3,3,24,18], ...
    'Color','w');

tiledlayout(fig2,6,1, ...
    'TileSpacing','compact', ...  
    'Padding','compact'); 

% subplot(a)
nexttile
plot(SRSf,SRSp1(1,:),'k-','LineWidth',1);
text(0.01,0.80,'(a)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);
text(0.16,0.80, ...
    sprintf(['Diagnostic cross-section: ','$Q_L=1$, $L=N=1000$, $K_{guess}=%d$'],2), ...
    'Interpreter','latex', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','bold');

% subplot(b)
nexttile
plot(SRSf,SRSp2(1,:),'k-','LineWidth',1);
text(0.01,0.80,'(b)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);
text(0.16,0.80, ...
    sprintf(['Diagnostic cross-section: ','$Q_L=1$, $L=N=1000$, $K_{guess}=%d$'],7), ...
    'Interpreter','latex', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','bold');

% subplot(c)
nexttile
plot(SRSf,SRSp3(1,:),'k-','LineWidth',1);
ylabel('Normalized Amplitude','FontName','Times New Roman','FontSize',12);
text(0.01,0.80,'(c)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);
text(0.16,0.80, ...
    sprintf(['Diagnostic cross-section: ','$Q_L=1$, $L=N=1000$, $K_{guess}=%d$'],10), ...
    'Interpreter','latex', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','bold');

% subplot(d)
nexttile
plot(SRSf,SRSp4(1,:),'k-','LineWidth',1);
text(0.01,0.80,'(d)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);
text(0.16,0.80, ...
    sprintf(['Diagnostic cross-section: ','$Q_L=1$, $L=N=1000$, $K_{guess}=%d$'],50), ...
    'Interpreter','latex', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','bold');

% subplot(e)
nexttile
plot(SRSf,SRSp5(1,:),'k-','LineWidth',1);
text(0.01,0.80,'(e)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);
text(0.16,0.80, ...
    sprintf(['Diagnostic cross-section: ','$Q_L=1$, $L=N=1000$, $K_{guess}=%d$'],100), ...
    'Interpreter','latex', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','bold');

% subplot(f)
nexttile
plot(SRSf,SRSp6(1,:),'k-','LineWidth',1);
xlabel('Frequency (Hz)','FontName','Times New Roman','FontSize',12)
text(0.01,0.80,'(f)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);
text(0.16,0.80, ...
    sprintf(['Diagnostic cross-section: ','$Q_L=1$, $L=N=1000$, $K_{guess}=%d$'],200), ...
    'Interpreter','latex', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','bold');

%% Figure (3)

fig3 = figure(3);clf(fig3);

set(fig3, ...
    'Units','centimeters', ...
    'Position',[3,3,36,12], ...
    'Color','w');

tiledlayout(fig3,1,2, ...
    'TileSpacing','compact', ...  
    'Padding','compact'); 

% subplot(a): normalized Fourier spectrum
nexttile;
plot(FFTf,NormFFTp,'k-','LineWidth',1);
ylim([0,1.2]);
xlabel('Frequency (Hz)','FontName','Times New Roman','FontSize',12);
ylabel('Normalized Amplitude','FontName','Times New Roman','FontSize',12);

text(0.04,0.90,'(a)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontWeight','bold', ...
    'FontSize',12);

xline([-0.45,-0.05,0.02,0.10],'r--');

legend('Fourier Spectrum', ...
    'Box','off', ...
    'FontName','Times New Roman', ...
    'FontSize',12);

% subplot(b): full 2D SRSA spectrum
nexttile;

pcolor(SRSf,Q_LVecTrue(:,2),SRSpTrue);
shading('interp');
colormap('pink');

xlabel('Frequency (Hz)','FontName','Times New Roman','FontSize',12);

ylabel('Subwindows $Q_L$','Interpreter','latex','FontName','Times New Roman','FontSize',12);

text(0.04,0.90,'(b)', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontWeight','bold', ...
    'FontSize',12, ...
    'Color','w');

text(0.45,0.90,'Full 2D Coherence Response', ...
    'Units','normalized', ...
    'FontName','Times New Roman', ...
    'FontWeight','bold', ...
    'FontSize',12, ...
    'Color','w');

text(0.70,0.20, ...
    sprintf('$K_{true}=%d$',K_true), ...
    'Interpreter','latex', ...
    'Units','normalized', ...
    'FontSize',12, ...
    'FontWeight','bold', ...
    'Color','w');