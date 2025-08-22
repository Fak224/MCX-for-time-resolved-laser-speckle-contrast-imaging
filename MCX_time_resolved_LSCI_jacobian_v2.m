%% MCXLAB Jacobian Simulation and Gated Analysis
% Author:Faezeh Akbari
% Date:08/22/2025
% Description:
%   - Runs MCXLAB photon simulations for multiple pulses
%   - Accumulates Jacobian data over pulses
%   - Computes LSCI and gated SPAD-like signals
%   - Saves results and generates visualization/video

clc; clear; close all;

%% ================== USER PARAMETERS ==================
nPulses     = 1000;          % Number of pulses to integrate
nPhotons    = 1e8;           % Photons per pulse
volSize     = [50, 50, 25];  % Volume dimensions
gateWidth   = 1277;          % Gate width (time bins)
outputDir   = 'C:\Users\fak224\OneDrive - University of Kentucky\Desktop\mcxlab\Run #19 Wide field\SPAD Proposal';

%% ================== INITIALIZE ==================
jac_data_integ = zeros([volSize, 750]);   % Pre-allocate Jacobian storage
jac_gate_spad  = [];

%% ================== SIMULATION LOOP ==================
for r = 1:nPulses
    fprintf('Running pulse %d/%d...\n', r, nPulses);

    % ----- Configure MCXLAB -----
    cfg.nphoton   = nPhotons;
    cfg.vol       = uint8(ones(volSize));
    cfg.prop      = [0 0 1 1; 0.003 10 0.9 1.37];  % [mua mus g n]
    cfg.issrcfrom0= 1;
    cfg.srctype   = 'cone';
    cfg.srcparam1 = [0.41 0 0 0];  % 20° in radians ≈ 0.3491
    cfg.srcpos    = [0 25 -40];

    theta         = 28 * pi / 180; % Tilt angle
    cfg.srcdir    = [sin(theta) 0 cos(theta)];

    % Detector positions (grid in x-y, z=1)
    [x, y]        = meshgrid(2:48, 2:48);
    cfg.detpos    = [x(:), y(:), ones(numel(x),1), ones(numel(x),1)];

    cfg.vol(:,:,1)= 0;   % Absorbing boundary for reflectance
    cfg.issaveref = 1;
    cfg.gpuid     = 1;
    cfg.autopilot = 1;

    % Time gating
    cfg.tstart    = 0;
    cfg.tend      = 13.5e-9;
    cfg.tstep     = 18e-12;

    % ----- Run MCXLAB -----
    [fluence,detpt,~,seeds,~] = mcxlab(cfg);

    % ----- Replay for Jacobian -----
    newcfg              = cfg;
    newcfg.seed         = seeds.data;
    newcfg.outputtype   = 'jacobian';
    newcfg.detphotons   = detpt.data;
    flux2               = mcxlab(newcfg);

    jac_data            = flux2.data;
    jac_data_integ      = jac_data_integ + jac_data;
end

%% ================== SAVE RESULTS ==================
save(fullfile(outputDir, 'Integration_1000pulses_jacobian.mat'), 'jac_data_integ');

%% ================== LSCI COMPUTATION ==================
jac_lsci = sum(jac_data_integ, 4);

figure;
img = squeeze(log(jac_lsci(:,20,2:end)));
imagesc(imrotate(img, -90));
caxis([-10 1]); axis square; colorbar;
title('LSCI', 'FontSize', 12);
ylabel(colorbar, 'Jac', 'FontSize', 12);

%% ================== GATED SPAD-LIKE COMPUTATION ==================
fprintf('Computing gated signals...\n');
nFrames = size(jac_data_integ, 4) - gateWidth;

for i = 1:nFrames
    jac_gate_spad(:,:,:,i) = sum(jac_data_integ(:,:,:,i:i+gateWidth), 4);
end

%% ================== GENERATE VIDEO ==================
videoPath = fullfile(outputDir, sprintf('GateWidth_%d_Jacobian.mp4', gateWidth));
v = VideoWriter(videoPath, 'MPEG-4');
v.FrameRate = 10;
open(v);

fig = figure('Color', 'w', 'Position', [100, 100, 600, 400]);

for i = 1:nFrames
    imagesc(imrotate(squeeze(jac_gate_spad(:,20,:,i)), 90));
    caxis([0 0.5]); axis square; colorbar;
    title(sprintf('Gate %d', i), 'FontSize', 12);
    ylabel(colorbar, 'Jac', 'FontSize', 10);

    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(v);
fprintf('Video saved to: %s\n', videoPath);
