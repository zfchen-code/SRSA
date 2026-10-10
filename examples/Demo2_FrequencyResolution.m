%----------------------------------------------------------------------------------------------------
% "What we observe is not nature itself, but nature exposed to our method of questioning."
%                                                        — Werner Heisenberg
%----------------------------------------------------------------------------------------------------
% *[Overview]
%  This script illustrates the high frequency-resolution capability of SRSA using two real-valued
%  multitone signals composed of five and nine closely spaced sinusoidal components, respectively.
%  Because each real-valued sinusoid corresponds to a conjugate pair of complex exponentials, the
%  ground-truth model orders are K1_true = 10 for y1 and K2_true = 18 for y2.
%
%  This demo reproduces Figure 2(e, g) of the manuscript. Its analysis parameters are set for the 
%  specified synthetic series. If the sampling interval or synthetic series is changed, reassess 
%  the input model order and spectral stacking interval, and update the corresponding row indices.
%
% *[Caption for Figure (1)]
%      Panel (a): Simulated multitone signal y1.
%      Panel (b): Singular-value spectrum of the Hankel matrix constructed from y1
%      Panel (c): Bilateral Fourier amplitude spectrum of y1.
%      Panel (d): Full 2-D SRSA coherence response of y1 computed with K_input = 10.
%      Panel (e): Selected Q_L range of the 2-D SRSA coherence response shown in panel (d).
%      Panel (f): Logarithmic stacked 1-D SRSA spectrum of y1.

% *[Caption for Figure (2)]
%      Panel (a): Simulated multitone signal y2.
%      Panel (b): Singular-value spectrum of the Hankel matrix constructed from y2
%      Panel (c): Bilateral Fourier amplitude spectrum of y2.
%      Panel (d): Full 2-D SRSA coherence response of y2 computed with K_input = 18.
%      Panel (e): Selected Q_L range of the 2-D SRSA coherence response shown in panel (d).
%      Panel (f): Logarithmic stacked 1-D SRSA spectrum of y2.
%
% *[Computational Environment]
%  The calculations were performed using MATLAB R2024b Update 3 on a workstation equipped with an
%  Intel Core i9-14900K processor (24 cores and 32 threads), 128 GB of RAM, and an NVIDIA GeForce
%  RTX 4090 D GPU with 24 GB of memory. The process-based parallel pool used eight workers. The
%  present implementation does not invoke GPU computing; therefore, the GPU did not contribute to
%  the reported runtimes.
%
% *[Representative Runtime]
%  This demo runs in a matter of seconds on the stated workstation. Runtime may vary with hardware
%  and parallel configuration.
%----------------------------------------------------------------------------------------------------
% Author       :   Zhifeng Chen
% E-mail       :   zfchen@whu.edu.cn
% Created      :   2026-01-20
% Last updated :   2026-06-20
% Last updated :   2026-10-10
%====================================================================================================
%% Data simulation with integer phase reduction
clc;
clear;

sampleDays = 30;                         % Positive integer sampling interval

validateattributes(sampleDays, {'double'}, ...
    {'scalar','real','finite','integer','positive','<=',flintmax});

DeltaT = sampleDays/365.25;

% Exact sample-count relation:
% n*sampleDays <= 41*365.25 = 59901/4 days.
% uint64 division avoids floating-point endpoint decisions.
nLast = double(idivide(uint64(59901), ...
    uint64(4)*uint64(sampleDays), 'floor'));

n = 0:nLast;
N = numel(n);

% t is used only for plotting.
t = 1982.5 + n*DeltaT;

% Exact integer arithmetic for all samples in this fixed time span.
% Multiplying sampleDays*n first keeps intermediates bounded.
timeNumerator    = 3965*1461 + 8*(sampleDays*n);
cycleDenominator = 2000*1461;

frequencyMilli1 = [1180,1185,1190,1195,1200];
phase1          = [1.1,0.0,0.0,0.1,1.1];

frequencyMilli2 = [2712,2723,2734,2745,2756,2767,2778,2789,2800];

y1 = zeros(1,N);
for k = 1:numel(frequencyMilli1)
    [c,s] = rationalPhaseTrig( ...
        frequencyMilli1(k),timeNumerator,cycleDenominator);

    if phase1(k) == 0
        component = c;
    else
        component = cos(phase1(k))*c - sin(phase1(k))*s;
    end

    y1 = y1 + component;
end

y2 = zeros(1,N);
for k = 1:numel(frequencyMilli2)
    c = rationalPhaseTrig( ...
        frequencyMilli2(k),timeNumerator,cycleDenominator);

    y2 = y2 + c;
end
%% Singular-value spectra

P = fix(N/2);                           % Number of Hankel columns
Q = N-P+1;                              % Number of Hankel rows

% Construct and decompose the Hankel matrix of y1.
H1 = zeros(Q,P);

for ii = 0:Q-1
    H1(ii+1,:) = y1(ii+1:ii+P);
end

singularV1  = svd(H1,'econ');           % Complete set of singular values
singularSN1 = (1:numel(singularV1)).';

clear H1

% Construct and decompose the Hankel matrix of y2.
H2 = zeros(Q,P);

for ii = 0:Q-1
    H2(ii+1,:) = y2(ii+1:ii+P);
end

singularV2  = svd(H2,'econ');           % Complete set of singular values
singularSN2 = (1:numel(singularV2)).';

clear H2

%% Fourier amplitude spectra

R1    = 200;                             
R2    = 200; 
Nfft = R1*N;
win  = hann(N);                         % Spectral window

% Bilateral frequency vector
FFTf = (-Nfft/2:Nfft/2-1).'/(Nfft*DeltaT);

% Bilateral spectra with exact window coherent-gain correction
FFTp1 = 2*abs(fftshift(fft(y1(:).*win,Nfft)))/sum(win);
FFTp2 = 2*abs(fftshift(fft(y2(:).*win,Nfft)))/sum(win);

%% SRSA spectra

tic;

Deltaf1 = 1/(DeltaT*N*R1);
Deltaf2 = 1/(DeltaT*N*R2);

bins1 = 1.17:Deltaf1:1.215;
bins2 = 2.60:Deltaf2:3.000;

K1_input = 10;%K_elbow
K2_input = 18;%K_elbow

[Bins1,LVEC1,SRSp1] = DoSRSA(y1,DeltaT,bins1,K1_input,'AlgorithmMode','svd');

[Bins2,LVEC2,SRSp2] = DoSRSA(y2,DeltaT,bins2,K2_input,'AlgorithmMode','svd');

toc;

displayQL = [200,300];
stackRows = 21:31;

stackQL1 = LVEC1(stackRows([1,end]),2);
stackQL2 = LVEC2(stackRows([1,end]),2);

%% Figure (1) 

fig1 = figure(1);clf(fig1);

set(fig1, ...
    'Units','centimeters', ...
    'Position',[3,3,36,20], ...
    'Color','w');

tiledlayout(fig1,2,3, ...
    'TileSpacing','compact', ...  
    'Padding','compact');        
colormap('pink');

% Subplot (a): simulated sequence y1

nexttile;

plot(t,y1,'k-');

xlabel('Time (yr)','FontName','Times New Roman','FontSize',10);

ylabel('Amplitude','FontName','Times New Roman','FontSize',10);

legend('Simulated Sequence 1 ($y_1$)', ...
    'Box','off', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

text(0.04,0.92,'(a)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

% Subplot (b): singular-value spectrum of y1

nexttile;

ld11 = plot(singularSN1,singularV1,'d-');
hold on;

xlim([8,30]);
ylim([-1,5]*1e-11);

xlabel('Singular-Value Index','FontName','Times New Roman','FontSize',10);

ylabel('Singular Value','FontName','Times New Roman','FontSize',10);

text(0.04,0.92,'(b)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

idx = singularSN1 == 11;

ld12 = plot(singularSN1(idx),singularV1(idx),'d', ...
    'MarkerEdgeColor','r', ...
    'MarkerFaceColor','r', ...
    'MarkerSize',4, ...
    'LineWidth',1.2);

hold off;

legend([ld11,ld12], ...
    'Singular-Value Spectrum','Spectrum Elbow', ...
    'Box','off', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

set(gca,'Layer','top');

% Subplot (c): Fourier spectrum of y1

nexttile;

plot(FFTf,FFTp1,'k-');
xline(1.180:0.005:1.2,'c--');
xlim([1.06,1.31]);
ylim([0,5]);

ylabel('Amplitude','FontName','Times New Roman','FontSize',10);

text(0.20,0.90,'$y_1$ Fourier Spectrum', ...
    'Units','normalized', ...
    'Color','k', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

text(0.06,0.90,'(c)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

% Subplot (d): full 2D SRSA spectrum of y1

nexttile;

pcolor(Bins1,LVEC1(:,2),log(SRSp1));
shading interp;
hold on;

rectangle('Position',[1.17,displayQL(1),1.21-1.17,displayQL(2)-displayQL(1)], ...
          'EdgeColor','c', ...
          'LineWidth',1.5, ...
          'LineStyle','-');

hold off;

colorbar('northoutside','Color','k');

xlim([1.17,1.21]);

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',10);

ylabel('$Q_L$','Interpreter','latex','FontName','Times New Roman','FontSize',10);

text(0.04,0.92,'(d)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10, ...
    'Color','w');

text(0.10,0.87,'$y_1$ Coherence Response', ...
    'Units','normalized', ...
    'Color','m', ...
    'Interpreter','latex', ...
    'FontSize',10);

text(0.28,0.80,sprintf('($K_{input}=%d$)',K1_input), ...
    'Units','normalized', ...
    'Color','m', ...
    'Interpreter','latex', ...
    'FontSize',10);

text(0.02,0.56,'Selected trade-off range', ...
    'Units','normalized', ...
    'Color','c', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

% Subplot (e): selected trade-off range of y1

nexttile;

pcolor(Bins1,LVEC1(:,2),log(SRSp1));
shading interp;

colorbar('northoutside','Color','k');

xlim([1.17,1.21]);
ylim(displayQL);

ylabel('$Q_L$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

text(0.06,0.90,'(e)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10, ...
    'Color','w');

text(0.17,0.90,'$y_1$ Coherence Response', ...
    'Units','normalized', ...
    'Color','m', ...
    'Interpreter','latex', ...
    'FontSize',10);

text(0.10,0.10, ...
    sprintf('$Q_L \\in [%d,%d]$', ...
    displayQL(1),displayQL(2)), ...
    'Units','normalized', ...
    'Color','c', ...
    'Interpreter','latex', ...
    'FontWeight','bold');

% Subplot (f): stacked 1D SRSA spectrum of y1

nexttile;

stackedSRSp1 = mean(SRSp1(stackRows,:),1);
plot(Bins1,log(stackedSRSp1),'k-');
xline(1.180:0.005:1.2,'c--');
xlim([1.17,1.21]);
ylim([-30,5]);

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',10);

ylabel('Logarithmic Coherence Response','FontName','Times New Roman','FontSize',10);

text(0.15,0.90,'$y_1$ 1D SRSA Spectrum', ...
    'Units','normalized', ...
    'Color','k', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

text(0.06,0.90,'(f)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

text(0.15,0.10, ...
    sprintf('Stacking $Q_L \\in [%d,%d]$', ...
    stackQL1(1),stackQL1(2)), ...
    'Units','normalized', ...
    'Color','k', ...
    'Interpreter','latex', ...
    'FontWeight','bold');

%% Figure (2) y2

fig2 = figure(2);clf(fig2);

set(fig2, ...
    'Units','centimeters', ...
    'Position',[3,3,36,20], ...
    'Color','w');

tiledlayout(fig2,2,3, ...
    'TileSpacing','compact', ...  
    'Padding','compact');        
colormap('pink');

% Subplot (a): simulated sequence y2

nexttile;

plot(t,y2,'k-');

xlabel('Time (yr)','FontName','Times New Roman','FontSize',10);

ylabel('Amplitude','FontName','Times New Roman','FontSize',10);

legend('Simulated Sequence 2 ($y_2$)', ...
    'Box','off', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

text(0.04,0.92,'(a)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

% Subplot (b): singular-value spectrum of y2

nexttile;

ld21 = plot(singularSN2,singularV2,'d-');
hold on;

xlim([16,38]);
ylim([-1,5]*1e-11);

xlabel('Singular-Value Index','FontName','Times New Roman','FontSize',10);

ylabel('Singular Value','FontName','Times New Roman','FontSize',10);

text(0.04,0.92,'(b)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

idx = singularSN2 == 19;

ld22 = plot(singularSN2(idx),singularV2(idx),'d', ...
    'MarkerEdgeColor','r', ...
    'MarkerFaceColor','r', ...
    'MarkerSize',4, ...
    'LineWidth',1.2);

hold off;

legend([ld21,ld22], ...
    'Singular-Value Spectrum','Spectrum Elbow', ...
    'Box','off', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

set(gca,'Layer','top');

% Subplot (c): Fourier spectrum of y2

nexttile;

plot(FFTf,FFTp2,'k-');
xline(2.712:0.011:2.800,'c--');
xlim([2.5,3.0]);
ylim([0,5]);

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',10);

ylabel('Amplitude','FontName','Times New Roman','FontSize',10);

text(0.20,0.90,'$y_2$ Fourier Spectrum', ...
    'Units','normalized', ...
    'Color','k', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

text(0.06,0.90,'(c)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

% Subplot (d): full 2D SRSA spectrum of y2

nexttile;

pcolor(Bins2,LVEC2(:,2),log(SRSp2));
shading interp;
hold on;

rectangle('Position',[2.66,displayQL(1),2.85-2.66,displayQL(2)-displayQL(1)], ...
          'EdgeColor','c', ...
          'LineWidth',1.5, ...
          'LineStyle','-');

hold off;

colorbar('northoutside','Color','k');

xlim([2.66,2.85]);

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',10);

ylabel('$Q_L$','Interpreter','latex','FontName','Times New Roman','FontSize',10);

text(0.04,0.92,'(d)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10, ...
    'Color','w');

text(0.10,0.87,'$y_2$ Coherence Response', ...
    'Units','normalized', ...
    'Color','m', ...
    'Interpreter','latex', ...
    'FontSize',10);

text(0.28,0.80,sprintf('($K_{input}=%d$)',K2_input), ...
    'Units','normalized', ...
    'Color','m', ...
    'Interpreter','latex', ...
    'FontSize',10);

text(0.02,0.56,'Selected trade-off range', ...
    'Units','normalized', ...
    'Color','c', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

% Subplot (e): selected trade-off range of y2

nexttile;

pcolor(Bins2,LVEC2(:,2),log(SRSp2));
shading interp;

colorbar('northoutside','Color','k');

xlim([2.70,2.81]);
ylim(displayQL);

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',10);

ylabel('$Q_L$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',10);

text(0.06,0.90,'(e)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10, ...
    'Color','w');

text(0.17,0.90,'$y_2$ Coherence Response', ...
    'Units','normalized', ...
    'Color','m', ...
    'Interpreter','latex', ...
    'FontSize',10);

text(0.10,0.10, ...
    sprintf('$Q_L \\in [%d,%d]$', ...
    displayQL(1),displayQL(2)), ...
    'Units','normalized', ...
    'Color','c', ...
    'Interpreter','latex', ...
    'FontWeight','bold');

% Subplot (f): stacked 1D SRSA spectrum of y2

nexttile;

stackedSRSp2 = mean(SRSp2(stackRows,:),1);
plot(Bins2,log(stackedSRSp2),'k-');
xline(2.712:0.011:2.800,'c--');
xlim([2.70,2.81]);
ylim([-22,5]);

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',10);

ylabel('Logarithmic Coherence Response','FontName','Times New Roman','FontSize',10);

text(0.15,0.90,'$y_2$ 1D SRSA Spectrum', ...
    'Units','normalized', ...
    'Color','k', ...
    'FontName','Times New Roman', ...
    'FontSize',10, ...
    'Interpreter','latex');

text(0.06,0.90,'(f)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

text(0.15,0.10, ...
    sprintf('Stacking $Q_L \\in [%d,%d]$', ...
    stackQL2(1),stackQL2(2)), ...
    'Units','normalized', ...
    'Color','k', ...
    'Interpreter','latex', ...
    'FontWeight','bold');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Local function
function [c,s] = rationalPhaseTrig(F,timeNumerator,D)
% Compute cos(2*pi*F*timeNumerator/D) and the corresponding sine.
% All phase reduction is performed with exact integer arithmetic.

    r = mod(F*timeNumerator,D); % Integer remainder: 0 <= r < D

    % Nearest quarter-cycle: q = 0,1,2,3,4.
    % D is even, so D/2 is an exact integer.
    q = floor((4*r + D/2)/D);

    % Exact integer residual; floating-point division occurs only here.
    residual = 4*r - q*D;
    x = residual/(2*D);         % -0.25 <= x <= 0.25

    cx = cospi(x);
    sx = sinpi(x);

    c = zeros(size(r));
    s = zeros(size(r));
    quadrant = mod(q,4);

    idx = quadrant == 0;
    c(idx) =  cx(idx); s(idx) =  sx(idx);

    idx = quadrant == 1;
    c(idx) = -sx(idx); s(idx) =  cx(idx);

    idx = quadrant == 2;
    c(idx) = -cx(idx); s(idx) = -sx(idx);

    idx = quadrant == 3;
    c(idx) =  sx(idx); s(idx) = -cx(idx);
end