% two_oscillator_sym_d.m
%
% This code produces those dynamical regimes possible with two units,
% with a firing-rate model of a neural circuit.
% Two units (representing neural
% assemblies) are simulated with firing rates that increase linearly with
% input current above a system-dependent threshold (Ithresh) and that saturate at a value
% rmax = 100Hz.
% Units can excite or inhibit each other to different degrees according to
% the connectivity matrix W in which W(i,j) represents connection strength
% from unit j to unit i.
% The flag "attractor-flag" is used to determine the type of simulation.
%
% This code is a solution of Tutorial 7.2 in the textbook,
% An Introductory Course in Computational Neuroscience
% by Paul Miller, Brandeis University, 2017
%%

clear

tmax = 2;  % default value of maximum time to simulate
Ntrials = 4;

% case 7 produces two distinct oscillators of different
% frequencies. Changes in initial conditions -- or a current pulse
% -- can switch between oscillators.
Ithresh = [-15; -12; 1; 2];
W = [1.2 -0.5 0.5 0; -1.5 1.9 0 1.2; -1.1 -0.6 0 -0.8; -1.5 -1.3 -0.8 0];
rinit1 = [0; 15; 0; 5];

Ntrials = 4;
Nunits = 4;

rng(1)

%% Figure positions for subplots
set(0,'DefaultLineLineWidth',2,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultAxesFontSize',14,...
    'DefaultAxesFontWeight','Bold');
figure(1)
clf

%% Set up the time vector
dt = 0.0001;
tvec = 0:dt:tmax;
Nt = length(tvec);

r = zeros(Nunits,Nt);   % array of rate of each cell as a function of time
tau = 0.010;            % base time constant for changes of rate

sigma_noise = 0.6/sqrt(dt);

%% Set up axes and label them for the figures

for trial = 1:Ntrials
    Inoise = sigma_noise*randn(size(r));
    %     if ( trial == 1 )
    r(:,1) = rinit1;                % Initialize firing rate
    %     else
    %         r(:,1) = rinit2;                % Initialize firing rate
    %     end
    %% Now simulate through time
    for i=2:Nt
        I = W'*r(:,i-1) + Inoise(:,i-1);                   % total current to each unit
        newr = r(:,i-1) + dt/tau*(I-Ithresh-r(:,i-1));  % Euler-Mayamara update of rates
        r(:,i) = max(newr,0);                           % rates are not negative        
    end
    
    
    subplot(Ntrials,1,trial)
    plot(tvec,r(1,:),'m')
    hold on
    plot(tvec,r(4,:),'b:')
    plot(tvec,r(2,:),'r')
    plot(tvec,r(3,:),'c:')
    axis([0 1.5 0 75])

    
end
