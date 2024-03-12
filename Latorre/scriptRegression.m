%
%**** This script determines the parameters of a CMM-based bilayered
%     4-fiber thin-wall G&R model of the arterial wall from biaxial data
%     collected at two different timepoints: original (o) and evolved (h).
%     The original homeostatic state serves as reference for the evolution.
%     Parameters in PAR (e.g., for elastin) remain fixed during evolution.
%     Parameters in parT/parC are allowed to evolve (can remain fixed too).
%
%**** The script can also be used to determine parameters of a CMM-based
%     hyperelastic model of the arterial wall if both data sets coincide,
%     that is, from data collected at a single timepoint (e.g., o or h).
%     In this case, the referential homeostatic state is case-specific.
%     All parameters (in particular, those in PAR) are also case-specific.
%
%     (note: some parameters, in particular deposition stretches, may reach
%            bound contraints, which could then be extended further.
%            If bound constraints are already realistic, potential
%            convergence-related warnings of the finite-difference squeme
%            may be turned off with the command <warning('off','last')>.
%            Constraining bounds for these cases is equivalent to prescribe
%            the parameter itself (e.g., Gc), if known from experiments.
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
clear all
%
global parT0 parC0 parT4 parC4
global PLexp0 PLexp4
global Gmmax Gcmax Gmmin Gcmin
%
Gmmin = 1.00; Gcmin = 1.00;
Gmmax = 1.25; Gcmax = 1.25;
%
mmHg_to_kPa = 0.13332;
%
%#################################  DTA  ##################################
%
protocols = 3;								% number of protocols for both Pd and Fl tests: 1, 2, 3
%
%------------------------------  DAY 0 (o) --------------------------------
%
rotf0 = 0.4440;								% outer radius at traction-free configuration (tf) [mm]
htf   = 0.1120;								% total thickness at tf [mm]
ritf  = rotf0-htf;							% inner radius at tf [mm]
lzivo = 1.62;								% in-vivo axial stretch (from tf)
%
SM = 0.7;									% medial cross-sectional area fraction from histology
SA = 0.3;									% adventitial cross-sectional area fraction from histology
%
rMAtf = sqrt(SA*ritf^2+SM*rotf0^2);     	% M-A interface radius at tf [mm]
%
meanpar0 = [18.536,16.593,0.108,25.37,0.036,0.078,1.719,0.5024];	% [c,c1t,c2t,c1z,c2z,c1d,c2d,alp]
Po = 106*mmHg_to_kPa;						% original homeostatic pressure [kPa]
%
dLPd = MeandL(meanpar0,[Po lzivo],[ritf rotf0]);	% in-vivo outer diameter and force
%
roo  = dLPd(1)/2;									% outer radius at o [mm] (consistent with mean material parameters, likely different to EXPER(R,5,T)/2*1e-3)
rMAo = sqrt(roo^2+1/lzivo*(rMAtf^2-rotf0^2));		% M-A radius at o [mm]
rio  = sqrt(roo^2+1/lzivo*(ritf^2-rotf0^2));		% inner radius at o [mm]
%
riio = [rio rMAo roo];								% in-vivo radii at o
%
[Pd0,Ll0] = BiaxialDataMean(meanpar0,[ritf rotf0 lzivo],protocols);		% generate biaxial testing data
%
Pd0(:,4) = Pd0(:,4)/lzivo;					% reset axial stretch from o
Ll0(:,4) = Ll0(:,4)/lzivo;					% reset axial stretch from o
%
dlexp0 = [Pd0(:,[2 4]);Ll0(:,[2 4])];		% experimental diameter-stretch pairs from tests Pd and Ll
PLexp0 = [Pd0(:,[1 3]);Ll0(:,[1 3])];		% experimental Pressure-Force pairs from tests Pd and Ll
%
emcMo = [0.4714    0.4714    0.0572];		% local mass fractions of medial [elastin muscle collagen]
ecAo  = [0.0333    0.9667];					% local mass fractions of adventitial [elastin collagen]
%
emco = [emcMo ecAo];						% local mass fractions [eM mM cM eA cA] at o
%
%------------------------------  DAY 28 (h) -------------------------------
%
rotf4 = 0.5745;								% outer radius at traction-free configuration (tf) [mm]
htf   = 0.2100;								% total thickness at tf [mm]
ritf  = rotf4-htf;							% inner radius at tf [mm]
lzivh = 1.34;								% in-vivo axial stretch (from tf)
%
SM = 0.4301;								% medial cross-sectional area fraction from histology
SA = 0.5699;								% adventitial cross-sectional area fraction from histology
%
rMAtf = sqrt(SA*ritf^2+SM*rotf4^2);     	% M-A interface radius at tf [mm]
%
meanpar4 = [20.652,7.593,1.285,1.074,3.821,0.109,7.154,0.6058];		% [c,c1t,c2t,c1z,c2z,c1d,c2d,alp]
Ph = 144*mmHg_to_kPa;						% evolved homeostatic pressure [kPa]
%
dLPd = MeandL(meanpar4,[Ph lzivh],[ritf rotf4]);	% in-vivo outer diameter and force
%
roh  = dLPd(1)/2;									% outer radius at h [mm] (consistent with mean material parameters, likely different to EXPER(R,5,T)/2*1e-3)
rMAh = sqrt(roh^2+1/lzivh*(rMAtf^2-rotf4^2));		% M-A radius at h [mm]
rih  = sqrt(roh^2+1/lzivh*(ritf^2-rotf4^2));		% inner radius at h [mm]
%
riih = [rih rMAh roh];								% in-vivo radii at h
%
[Pd4,Ll4] = BiaxialDataMean(meanpar4,[ritf rotf4 lzivh],protocols);		% generate biaxial testing data
%
Pd4(:,4) = Pd4(:,4)/lzivh;					% reset axial stretch from h
Ll4(:,4) = Ll4(:,4)/lzivh;					% reset axial stretch from h
%
dlexp4 = [Pd4(:,[2 4]);Ll4(:,[2 4])];		% experimental diameter-stretch pairs from tests Pd and Ll
PLexp4 = [Pd4(:,[1 3]);Ll4(:,[1 3])];		% experimental Pressure-Force pairs from tests Pd and Ll
%
emcMh = [0.2799    0.6271    0.0930];		% local mass fractions of medial [elastin muscle collagen]
ecAh  = [0.0063    0.9937];					% local mass fractions of adventitial [elastin collagen]
%
emch = [emcMh ecAh];						% local mass fractions [eM mM cM eA cA] at h
%
%------------------------------  REGRESSION  ------------------------------
%
c   = 90;									% c elastin
Get = 1.9;									% elastin circumferential deposition stretch in media and adventitia
Gez = 1.6;									% elastin axial deposition stretch in media and adventitia
Bt  = 1/4/4;								% fraction of circumferential collagen within each layer
Bz  = 1/4/4;								% fraction of axial collagen within the adventitia
alp = pi/2/3;								% orientation of diagonal collagen wrt axial direction [rad]
%
PAR = [c Get Gez Bt Bz alp];				% initial guesses for PAR
%
lb = [0.0 1.00 1.00 0.0 0.0 0.00];			% lower bounds  /  [0.0 1.00 1.00 0.0 0.0 0.00]  /  0.9*PAR
ub = [inf 2.25 2.25 1/2 1/2 pi/2];			% upper bounds  /  [inf 2.25 2.25 1/2 1/2 pi/2]  /  1.1*PAR
%
parT0 = [100 10 100 10 1.10 1.10];			% initial guesses for parT0
parT4 = parT0;								% initial guesses for parT4
%
%** perform least squares curve fitting to determine PAR, parT, and parC
%
dlexp = [dlexp0;dlexp4];					% diameter-stretch pairs from ALL tests Pd and Ll
PLexp = [PLexp0;PLexp4];					% Pressure-Force pairs from ALL tests Pd and Ll
%
N0 = size(dlexp0,1);						% P-d and L-l measurements at day 0
setpar = 2;									% sets of parameters [parT parC] to determine: 1, 2
%
options = optimoptions(@lsqcurvefit,'Display','None','FunctionTolerance',1e-6,'StepTolerance',1e-6,'OptimalityTolerance',1e-6);
[PAR,rnorm,~,eflag,~] = lsqcurvefit(@(PAR,dl) PreForGR(PAR,dl,emco,riio,emch,riih,N0,setpar),PAR,dlexp,PLexp,lb,ub,options);
clear options
%
fprintf('PAR =\n%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',PAR)
fprintf('parTC0 =\n%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',[parT0,parC0])
fprintf('parTC4 =\n%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',[parT4,parC4])
%
%% -------------------------------  PLOTS  --------------------------------
%
N4 = size(dlexp4,1);						% P-d and L-l measurements at week 4
%
PL0 = PreForGRTC(parC0,dlexp0,PAR,parT0,emco,riio,emco,riio,N0);	% values at day 0
PL4 = PreForGRTC(parC4,dlexp4,PAR,parT4,emch,riih,emco,riio,N4);	% values at day 28
%
figure
%
subplot(221)
hold on
grid on
plot(Pd0(:,2),PLexp0(1:size(Pd0,1),1)/mmHg_to_kPa,'ok','markersize',5,'linew',0.75)
plot(Pd0(:,2),PL0(1:size(Pd0,1),1)/mmHg_to_kPa,'ok','markersize',3,'markerfacecolor','k','markeredgecolor','none')
plot(Pd4(:,2),PLexp4(1:size(Pd4,1),1)/mmHg_to_kPa,'sk','markersize',5,'linew',0.75)
plot(Pd4(:,2),PL4(1:size(Pd4,1),1)/mmHg_to_kPa,'sk','markersize',3,'markerfacecolor','k','markeredgecolor','none')
xlabel 'Outer Diameter [mm]'
ylabel 'Pressure [mmHg]'
set(gca,'xlim',[0.6 1.5],'ylim',[0 160])
set(gca,'fontsize',11)
%
subplot(222)
hold on
grid on
plot(PLexp0(1:size(Pd0,1),1)/mmHg_to_kPa,PLexp0(1:size(Pd0,1),2),'ok','markersize',5,'linew',0.75)
plot(PL0(1:size(Pd0,1),1)/mmHg_to_kPa,PL0(1:size(Pd0,1),2),'ok','markersize',3,'markerfacecolor','k','markeredgecolor','none')
plot(PLexp4(1:size(Pd4,1),1)/mmHg_to_kPa,PLexp4(1:size(Pd4,1),2),'sk','markersize',5,'linew',0.75)
plot(PL4(1:size(Pd4,1),1)/mmHg_to_kPa,PL4(1:size(Pd4,1),2),'sk','markersize',3,'markerfacecolor','k','markeredgecolor','none')
xlabel 'Pressure [mmHg]'
ylabel 'Transducer axial force [mN]'
set(gca,'xlim',[0 160],'ylim',[0 40])
set(gca,'fontsize',11)
%
subplot(223)
hold on
grid on
plot(lzivo*Ll0(:,4),PLexp0(size(Pd0,1)+1:end,2),'ok','markersize',5,'linew',0.75)
plot(lzivo*Ll0(:,4),PL0(size(Pd0,1)+1:end,2),'ok','markersize',4,'markersize',3,'markerfacecolor','k','markeredgecolor','none')
plot(lzivh*Ll4(:,4),PLexp4(size(Pd4,1)+1:end,2),'sk','markersize',5,'linew',0.75)
plot(lzivh*Ll4(:,4),PL4(size(Pd4,1)+1:end,2),'sk','markersize',4,'markersize',3,'markerfacecolor','k','markeredgecolor','none')
xlabel 'Axial stretch'
ylabel 'Transducer axial force [mN]'
set(gca,'xlim',[1.2 1.8],'ylim',[0 40])
set(gca,'fontsize',11)
%
subplot(224)
hold on
grid on
plot(lzivo*Ll0(:,4),PLexp0(size(Pd0,1)+1:end,1)/mmHg_to_kPa,'ok','markersize',5,'linew',0.75)
plot(lzivo*Ll0(:,4),PL0(size(Pd0,1)+1:end,1)/mmHg_to_kPa,'ok','markersize',3,'markerfacecolor','k','markeredgecolor','none')
plot(lzivh*Ll4(:,4),PLexp4(size(Pd4,1)+1:end,1)/mmHg_to_kPa,'sk','markersize',5,'linew',0.75)
plot(lzivh*Ll4(:,4),PL4(size(Pd4,1)+1:end,1)/mmHg_to_kPa,'sk','markersize',3,'markerfacecolor','k','markeredgecolor','none')
xlabel 'Axial stretch'
ylabel 'Pressure [mmHg]'
set(gca,'xlim',[1.2 1.8],'ylim',[0 160])
set(gca,'fontsize',11)
%
%% -----------------  Axial Stretches / Bilayer stresses  -----------------
%
%** compute rotf and lziv at day 0
%
options = optimoptions(@fsolve,'Display','none','FunctionTolerance',eps,'StepTolerance',eps,'OptimalityTolerance',eps);
%
lzoo = 1;   % axial stretch from o to o (i.e., 1)
rolz0 = fsolve(@(rl) [0;0]-LaplaceAxial(rl,PAR,parT0,parC0,riio,emco,riio,lzoo),[rotf0,1/lzivo],options);
lzlz0 = [lzivo 1/rolz0(2)]   % exp. lziv, pred. lziv (from tfo to o)
%
%** compute rotf and lziv at day 28
%
lzoh = 1;   % axial stretch from o to h (assumed 1)
rolz4 = fsolve(@(rl) [0;0]-LaplaceAxial(rl,PAR,parT4,parC4,riio,emch,riih,lzoh),[rotf4,1/lzivh],options);
lzlz4 = [lzivh 1/rolz4(2)]   % exp. lziv, pred. lziv (from tfh to h)
%
clear options
%
%** compute layer specific stresses
%
PreStrSc
%
%** save computed values for further computations in evolution form
%
save('DTA_nlr.mat',...
	 'PAR','parT0','parC0','parT4','parC4',...
	 'emco','riio','so','PPfo','rotf0','lzivo','rolz0',...
	 'emch','riih','sh','PPfh','rotf4','lzivh','rolz4')
%