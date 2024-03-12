
% David Li | Yale University | Feb 2024
% Regenerates synthetic biaxial data from fitted 4FF parameters and fits
% with bilayer pre-stretch 4FF model for G&R simulations.

clear; clc; %close all;
set(0,'DefaultAxesLineWidth',1);
set(0,'DefaultLineLineWidth',1);
set(0,'DefaultAxesFontSize',12);

global parT0 parC0 PLexp0
global Gmmax Gcmax Gmmin Gcmin c2mmax c2cmax

% upper and lower bounds for SMC & collagen deposition stretch & moduli
Gmmin = 1.00; Gcmin = 1.00;
Gmmax = 1.20; Gcmax = 1.25;
c2mmax = 2e1; c2cmax = 2e1;

mmHg_to_kPa = 0.133322;		% convert mmHg to kPa

%% Set protocols/pressures and specimen
protocols = 3;								% number of protocols for both Pd and Fl tests: 1, 2, 3
pres      = 3;								% experimental pressure (1 dias, 2 MAP, 3 sys)
num       = 'all';							% specimen number to fit

% generate and load synthetic biaxial data
path = pwd;
scriptSpecMats;
load([path,'\biax_C1041G_',num,'.mat']);
disp(['Specimen: C1041G_',num]);
disp(['Protocols: ',num2str(protocols)]);
disp(['Pressure: ',num2str(P(pres)),' mmHg']);

%% ------------------------------  ATA  ------------------------------
Po = P(pres)*mmHg_to_kPa;					% experimental pressure [kPa]
rotf0 = dotf/2;								% outer radius at traction-free configuration (tf) [mm]
lzivo = lziv;								% in-vivo axial stretch (from tf)
riio = rii;									% in-vivo radii at o

Pd0 = Pd;									% pressure-diameter pairs
Ll0 = Ll;									% force-stretch pairs
dlexp0 = dlexp;								% diameter-stretch pairs
PLexp0 = PLexp;								% pressure-force pairs

emcMo = emcM;								% local mass fractions of media
ecAo = ecA;									% local mass fractions of adventitia
emco = [emcMo ecAo];						% local mass fractions [eM mM cM eA cA] at o

%% ------------------------------  REGRESSION  ------------------------------
% % ATA fit (WT)
% c   = 90;									% c elastin
% Get = 1.9;									% elastin circumferential deposition stretch in media and adventitia
% Gez = 1.6;									% elastin axial deposition stretch in media and adventitia
% Bt  = 1/4/4;								% fraction of circumferential collagen within each layer
% Bz  = 1/4/4;								% fraction of axial collagen within the adventitia
% alp = pi/2/3;                               % orientation of diagonal collagen wrt axial direction [rad]

% % ATA fit (WT 2023)
% c   = 100;									% c elastin
% Get = 2.15;									% elastin circumferential deposition stretch in media and adventitia
% Gez = 1.87;									% elastin axial deposition stretch in media and adventitia
% Bt  = 0.0114;								% fraction of circumferential collagen within each layer
% Bz  = 0.0202;                               % fraction of axial collagen within the adventitia
% alp = deg2rad(40);                        % orientation of diagonal collagen wrt axial direction [rad]

% % ATA fit (C1041G->mgRall, no elastin evol 2023)
% c   = 55;									% c elastin
% Get = 1.80;									% elastin circumferential deposition stretch in media and adventitia
% Gez = 2.23;									% elastin axial deposition stretch in media and adventitia
% Bt  = 0.0130;								% fraction of circumferential collagen within each layer
% Bz  = 0.0050;								% fraction of axial collagen within the adventitia
% alp = deg2rad(50.5);							% orientation of diagonal collagen wrt axial direction [rad]

% % ATA fit (C1041G->C1041G 2023)
% c   = 80;									% c elastin
% Get = 2.20;									% elastin circumferential deposition stretch in media and adventitia
% Gez = 2.14;									% elastin axial deposition stretch in media and adventitia
% Bt  = 0.1500;								% fraction of circumferential collagen within each layer
% Bz  = 0.0160;								% fraction of axial collagen within the adventitia
% alp = deg2rad(48.1);						% orientation of diagonal collagen wrt axial direction [rad]

% % ATA fit (C1041G->mgRall, close to Marfan model paper)
% c   = 42;									% c elastin
% Get = 1.74;									% elastin circumferential deposition stretch in media and adventitia
% Gez = 2.24;									% elastin axial deposition stretch in media and adventitia
% Bt  = 0.0120;								% fraction of circumferential collagen within each layer
% Bz  = 0.0120;								% fraction of axial collagen within the adventitia
% alp = deg2rad(51);						% orientation of diagonal collagen wrt axial direction [rad]

% ATA parameters (R Soc paper)
c   = 42.03;									% c elastin
Get = 1.735;									% elastin circumferential deposition stretch in media and adventitia
Gez = 2.244;									% elastin axial deposition stretch in media and adventitia
Bt  = 0.0121;								% fraction of circumferential collagen within each layer
Bz  = 0.0117;								% fraction of axial collagen within the adventitia
alp = deg2rad(51.12);						% orientation of diagonal collagen wrt axial direction [rad]

PAR = [Get Gez Bt Bz alp];				% initial guesses for PAR
lb  = [1.00 1.00 0.0 0.0 0.00];			% lower bounds  /  [0.0 1.00 1.00 0.0 0.0 0.00]
ub  = [2.25 2.25 1/2 1/2 pi/2];			% upper bounds  /  [inf 2.25 2.25 1/2 1/2 pi/2]

% parT0 = [150 1 150 1 1.05 1.05 c];		% initial guesses for parT0 / [100 10 100 10 1.10 1.10] [150 10 150 10 1.255 1.25]
% parT0 = [100 1 100 1 1.20 1.05 c];		% initial guesses for parT0 / [100 10 100 10 1.10 1.10] [150 10 150 10 1.255 1.25]
% parT0 = [1.26 30 665 2.1 1.19 1.24 c];		% initial guesses for ATA USE THIS
parT0 = [1.89 30 1572 1.5 1.01 1.01 c];		% initial guesses for both ATA
% parT0 = [200 eps 3000 10 1.19 1.05 c];		% initial guesses for C1041G DTA
% parT0 = [360 1.5 1000 10 1.19 1.11 c];		% initial guesses for mgR DTA
% parT0 = [300 2 500 2 1.19 1.24 c];		% initial guesses for WT ATA
% parT0 = [480 0.15 815 1.15 1.1 1.2 c];		% initial guesses for WT ATA 2023
% parT0 = [1.26 29.99 665.6 2.14 1.19 1.2499 c];		% ATA R Soc paper

% perform least squares curve fitting to determine PAR, parT, and parC
dlexp = dlexp0;					% diameter-stretch pairs from ALL tests Pd and Ll
PLexp = PLexp0;					% Pressure-Force pairs from ALL tests Pd and Ll

N0 = size(dlexp0,1);						% P-d and L-l measurements at day 0
setpar = 2;									% sets of parameters [parT parC] to determine: 1, 2

options = optimoptions(@lsqcurvefit,'Display','iter',...
                                    'FunctionTolerance',1e-6,...
                                    'StepTolerance',1e-6,...
                                    'OptimalityTolerance',1e-6,...
                                    'FiniteDifferenceType','central',...
                                    'FiniteDifferenceStepSize',6e-4);
[PAR,rnorm,~,eflag,~] = lsqcurvefit(@(PAR,dl) PreForGR(PAR,dl,emco,riio,N0,setpar),PAR,dlexp,PLexp,lb,ub,options);
clear options

%                [Get    Gez    Bt     Bz     alp]
fprintf('PAR =\n %10.4f %10.4f %10.4f %10.4f %10.4f \n',[PAR(1:end-1) rad2deg(PAR(end))])
%                  [c1m_t  c2m_t  c1c_t  c2c_t  Gm     Gc     c      c1m_c  c2m_c  c1c_c  c2c_c]
fprintf('parTC =\n %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f \n',[parT0,parC0])

%% -------------------------------  PLOTS  --------------------------------

PL0 = PreForGRTC(parC0,dlexp0,PAR,parT0,emco,riio,N0);	% values C1041G

figure; set(gcf,'units','normalized','outerposition',[.3 .25 .4 .73])

subplot(221)
hold on; grid on
plot(Pd0(:,2),PLexp0(1:size(Pd0,1),1)/mmHg_to_kPa,'db','markersize',3,'linew',0.75)
plot(Pd0(:,2),PL0(1:size(Pd0,1),1)/mmHg_to_kPa,'xr','markersize',3)
xlabel 'Outer Diameter [mm]'
ylabel 'Pressure [mmHg]'
set(gca,'xlim',[0.2 3],'ylim',[-10 160]) % ATA
set(gca,'fontsize',11)
legend('Data','Fit','location','northwest')

subplot(222)
hold on; grid on
plot(PLexp0(1:size(Pd0,1),1)/mmHg_to_kPa,PLexp0(1:size(Pd0,1),2),'db','markersize',3,'linew',0.75)
plot(PL0(1:size(Pd0,1),1)/mmHg_to_kPa,PL0(1:size(Pd0,1),2),'xr','markersize',3)
xlabel 'Pressure [mmHg]'
ylabel 'Transducer axial force [mN]'
set(gca,'xlim',[-10 160],'ylim',[-20 50]) % ATA
set(gca,'fontsize',11)

subplot(223)
hold on; grid on
plot(lzivo*Ll0(:,4),PLexp0(size(Pd0,1)+1:end,2),'db','markersize',3,'linew',0.75)
plot(lzivo*Ll0(:,4),PL0(size(Pd0,1)+1:end,2),'xr','markersize',3)
xlabel 'Axial stretch'
ylabel 'Transducer axial force [mN]'
set(gca,'xlim',[1.2 2.0],'ylim',[0 50]) % ATA
set(gca,'fontsize',11)

subplot(224)
hold on; grid on
plot(lzivo*Ll0(:,4),PLexp0(size(Pd0,1)+1:end,1)/mmHg_to_kPa,'db','markersize',3,'linew',0.75)
plot(lzivo*Ll0(:,4),PL0(size(Pd0,1)+1:end,1)/mmHg_to_kPa,'xr','markersize',3)
xlabel 'Axial stretch'
ylabel 'Pressure [mmHg]'
set(gca,'xlim',[1.2 2.0],'ylim',[0 160]) % ATA
set(gca,'fontsize',11)

%% -----------------  Axial Stretches / Bilayer Stresses  -----------------

% compute rotf and lziv o
options = optimoptions(@fsolve,'Display','None','FunctionTolerance',eps,'StepTolerance',eps,'OptimalityTolerance',eps);

lzoo = 1;   % axial stretch wrt in vivo (i.e., 1)
rolz0 = fsolve(@(rl) [0;0]-LaplaceAxial(rl,PAR,parT0,parC0,riio,emco,riio,lzoo),[rotf0,1/lzivo],options);
lzlz0 = [lzivo 1/rolz0(2)];   % exp. lziv, pred. lziv (from tfo to o)

clear options

% compute layer-specific stresses
PreStrSc
% fprintf('od h =\n%10.0f %10.0f\n',riio(3)*2*1e3,(riio(3)-riio(1))*1e3);

% save computed values for further computations in evolution form
save([path,'\ATA_C1041G_',num,'_params.mat'],...
	 'PAR','parT0','parC0',...
	 'emco','riio','so','PPfo','rotf0','lzivo','rolz0');
