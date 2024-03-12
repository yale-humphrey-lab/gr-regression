%
%**** This script computes the mechanobiologically equilibrated evolution
%     of a CMM-based bilayered 4-fiber thin-wall G&R model of an artery.
%
%     It also computes pressure-diameter and axial force-length tests
%     performed in silico at selected timepoints, for which G&R is frozen.
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
clear all
%
load('DTA_nlr.mat')				% load computed values from nonlinear regression
%
mmHg_to_kPa = 0.13332;
%
scrsz = get(0,'ScreenSize');
% 
%** material data
%
fctTC = 1/3;					%* potential factor for parTC evolution
%
Dtau = 1/(riih(1)/riio(1))^3-1;	%* Dtau = (Qh/Qo)/(ah(1)/ao(1))^3 - 1 (>0)
Dsgm = sh(7)/so(7)-1;			%* Dsgm = sgmIh/sgmIo - 1             (<0)
%
KicM = 2.00;					%* assumed
KscM = 2.50;					%* assumed (if KscM=KicM*Dsgm/Dtau, as in IAA, then KscM<0, so KfcM needed if Qh/Qo=1!)
KfcM = KscM*Dtau-KicM*Dsgm;		%* determined as a result, it must be KfcM > 0
%
etaq = 1;						%* assume kom = koc, i.e. etaq := etaqm = kom/koc = 1 = koc/koc = etaqc 
%
etaUm = 0.8/etaq;				%* KimM/KicM = KsmM/KscM = KfmM/KfcM = etam/etaqm
etaUc = 5/3/etaq;				%* KicA/KicM = KscA/KscM = KfcA/KfcM = etac/etaqc
%
Keta = [KicM KscM KfcM etaUm etaUc etaq];	%* gain parameters and eta's
%
rho = 1050;						%* density of soft tissue
%
rio  = riio(1);					% inner radius at o
rMAo = riio(2);					% M-A radius at o
roo  = riio(3);					% outer radius at o
%
hMo = rMAo-rio;					% medial thickness at o
hAo = roo-rMAo;					% adventitial thickness at o
%
Tmax = 0;						%* WITHOUT active tone
CB = 0.8326; CS = 0.5*CB;
lM = 1.1; l0 = 0.4;
flact = 1-((lM-1)/(lM-l0))^2;	% with lM-lact = lM-rio/rio = lM-1
%
actpar = [Tmax flact CB CS];	% active stress parameters
%
days = 28*1;					% total G&R time [days]
sts = 5;                        % number of equil. states to compute
Ds = days/(sts-1);				% time step length (time is just a parameter)
s = (0:Ds:days)';				% G&R times
%
Po    = PPfo(2);				% inner pressure at o
sgmIo = so(7);					% first principal invariant of mean stress at o
fo    = PPfo(3);				% vessel axial force at o
%
Ph   = PPfh(2);					% inner pressure at h
fctr = Ph/Po;					% fold increase in pressure from o to h
%
P = 1 + (fctr-1)*s/days;		% P/Po, linear increase (only interested in initial & final timepoints)
%
P = P*Po;						% inner pressure evolution
Q = ones(size(P));				% constant cardiac output
lz = 1;							% constant axial stretch wrt state o
%
inf = s/days;                   % inflammation, linear increase (only interested in initial & final timepoints)
%
parT04 = parT0 + (inf.^fctTC)*(parT4-parT0);
parC04 = parC0 + (inf.^fctTC)*(parC4-parC0);
%
options = optimoptions(@fsolve,'Display','None','FunctionTolerance',eps,'StepTolerance',eps,'OptimalityTolerance',eps);
%
xsol = [rio hMo hAo emco(3) fo];	% initial guess for geometry, med. coll. mass fraction, axial force
%
figure
set(gcf,'position',[0.1*scrsz(3) 0.25*scrsz(4) 0.6*scrsz(3) 0.5*scrsz(4)])
%
for i = 1:1:length(s)
	%
	parT = parT04(i,:);
	parC = parC04(i,:);
	%
	load = [P(i) Q(i) lz inf(i)];
	%
	%** compute geometry, med. coll. mass fraction, axial force for given loads
	%
	[xsol,fval,eF] = fsolve(@(x) Thin2LEquilEvol(x,load,PAR,parT,emco,riio,sgmIo,Keta,actpar),xsol,options);
	%
	ri(i)    = xsol(1);         % inner radius
	hM(i)    = xsol(2);         % medial thickness
	hA(i)    = xsol(3);         % adv. thickness
	emc(i,3) = xsol(4);         % med. coll. mass fraction
	f(i)     = xsol(5);         % vessel axial force
	%
	rii = [ri(i) ri(i)+hM(i) ri(i)+hM(i)+hA(i)]; % inner, MA, outer radii
	%
	ltM = (2*ri(i)+hM(i))/(2*rio+hMo);               % circ. stretch media
	ltA = (2*ri(i)+2*hM(i)+hA(i))/(2*rio+2*hMo+hAo); % circ. stretch adv.
	lrM = hM(i)/hMo;                                 % radial stretch media
	lrA = hA(i)/hAo;                                 % radial stretch adv.
	%
	JM = ltM*lrM*lz;            % volume ratio media
	JA = ltA*lrA*lz;            % volume ratio adv.
	%
	emc(i,2) = emco(2)/JM*(JM*emc(i,3)/emco(3))^(etaUm*etaq);	% medial smc mass fraction (mass growth)
	emc(i,5) = emco(5)/JA*(JM*emc(i,3)/emco(3))^(etaUc*  1 );	% adventitial collagen mass fraction (mass growth)
	%
	emc(i,1) = emco(1)/JM;		% medial elastin mass fraction (no growth)
	emc(i,4) = emco(4)/JA;		% adventitial elastin mass fraction (no growth)
	%
	%** compute/plot tf configuration and in-vivo stretches
	%
	rolz = rolz0 + (P(i)-Po)/(Ph-Po)*(rolz4-rolz0);		%* interpolated initial guess
	%
	[rolz,fval,eF] = fsolve(@(rl) [0;0]-LaplaceAxial(rl,PAR,parT,parC,riio,emc(i,:),rii,1),rolz,options);
	rotf = rolz(1);				%* outer radii at tf (computed)
	lziv = lz/rolz(2);			%* in-vivo axial stretch (computed)
    %
    %   ... or
	%
	% 	rotfi = rolz0(1); lzivi = lz/rolz0(2);						%* rotf and lziv at 0 days
	% 	rotff = rolz4(1); lzivf = lz/rolz4(2);						%* rotf and lziv at 4 weeks
	% 	rotf  = rotf0 + (rotf-rotfi)/(rotff-rotfi)*(rotf4-rotf0);	%* outer radii at tf (experimental - interpolated)
	% 	lziv  = lzivo + (lziv-lzivi)/(lzivf-lzivi)*(lzivh-lzivo);	%* in-vivo axial stretch (experimental - interpolated)
	%
	ritf = sqrt(rotf^2+lziv*(ri(i)^2-(ri(i)+hM(i)+hA(i))^2));	%* inner radii at tf
	ltiv = (2*ri(i)+hM(i)+hA(i))/(ritf+rotf);					%* in-vivo circumferential stretch
	%
	subplot(333)
	hold on
	bar((i-1)/(sts-1)*100-2,ltiv,3,'FaceColor',[0.2 0.2 0.2])
	bar((i-1)/(sts-1)*100+2,lziv,3,'FaceColor',[0.7 0.7 0.7])
	ylabel '\lambda^{iv}'
	set(gca,'fontsize',11,'xlim',[-7 107],'ylim',[1 2],'xtick',s/s(end)*100,'xticklabel',num2str(s))	% 
	legend('\lambda_{\theta}^{iv}','\lambda_{z}^{iv}','location','northwest')
	%
	%** simulate/plot Pd and Ll tests at i-th day
	%
	j = 0;
	%
	for Pres = (10:5:140)*mmHg_to_kPa
		%
		j = j+1;
		%
		[roPd(j),fval,eF] = fsolve(@(r) Pres-Laplace(r,lz,PAR,parT,parC,riio,emc(i,:),rii,1),ri(i)+hM(i)+hA(i),options);	%* outer radius at P-d test
		%
		PLPd(:,j) = LaplaceAxial([roPd(j) lz],PAR,parT,parC,riio,emc(i,:),rii,1);    % Pressure and transducer force at P-d test
		%
	end
	%
	riPd = sqrt(roPd.^2+1/lz*(ri(i)^2-(ri(i)+hM(i)+hA(i))^2));	% inner radii at P-d test
	%
	subplot(231)
	hold on
	plot(2*roPd,PLPd(1,:)/mmHg_to_kPa,'linew',1.5)
	xlabel 'Outer Diameter [mm]'
	ylabel 'Pressure [mmHg]'
	set(gca,'fontsize',11)
	%
	subplot(232)
	hold on
	plot((riPd+roPd)/(ritf+rotf),PLPd(1,:).*riPd./(roPd-riPd),'linew',1.5)
	xlabel 'Circ. stretch [-]'
	ylabel 'Circ. stress [kPa]'
	set(gca,'fontsize',11)
	%
	j = 0;
	%
	Pres = P(i);				%* associated homeostatic pressure P(i), or 100*mmHg_to_kPa
	%
	for Trf = 2:2:38			%* transducer force
		%
		j = j+1;
		%
		LLl(j) = Trf;			%* transducer force at L-l test
		%
		[rolz,fval,eF] = fsolve(@(rl) [Pres;Trf]-LaplaceAxial(rl,PAR,parT,parC,riio,emc(i,:),rii,1),[ri(i)+hM(i)+hA(i),lz],options);
		%
		roLl(j) = rolz(1);		%* outer radius at L-l test
		lLl(j) = rolz(2);		%* axial stretch at L-l test
		%
	end
	%
	riLl = sqrt(roLl.^2+1/lz*(ri(i)^2-(ri(i)+hM(i)+hA(i))^2));	% inner radii at L-l test
	%
	subplot(234)
	hold on
	plot(lziv*lLl,LLl,'linew',1.5)
	xlabel 'Axial stretch [-]'
	ylabel 'Transducer axial force [mN]'
	set(gca,'fontsize',11)
	%
	subplot(235)
	hold on
	plot(lziv*lLl,(LLl+pi*riLl.^2*P(i))/pi/((ri(i)+hM(i)+hA(i))^2-ri(i)^2),'linew',1.5)
	xlabel 'Axial stretch [-]'
	ylabel 'Axial stress [kPa]'
	set(gca,'fontsize',11)
	%
	%** compute/plot stiffness and stored energy
	%
	ce  = PAR(1);  Get = PAR(2);  Gez = PAR(3);
	Bt  = PAR(4);  Bz  = PAR(5);  alp = PAR(6);
	c1m = parT(1); c2m = parT(2); c1c = parT(3); 
	c2c = parT(4); Gm  = parT(5); Gc  = parT(6);
	%
	betaM = [Bz 1-Bz];						% medial betas [BzM 2*BdM]
	betaA = [Bt Bz 1-Bt-Bz];				% adventitial betas [BtA BzA 2*BdA]
	%
	phiM = [emc(i,1:2) emc(i,3)*betaM];		% local mass fractions of medial [e mt cz 2*cd]
	phiA = [emc(i,4)   emc(i,5)*betaA];		% local mass fractions of adventitial [e ct cz 2*cd]
	%
	Geh = [1/Get/Gez*lrM 1/Get/Gez*lrA Get*ltM Get*ltA Gez*lz Gez*lz];	    %* Foh*Ge = [GerMh GerAh GetMh GetAh GezMh GezAh]
	%
	stM = phiM(1)*ce*(Geh(3)^2-Geh(1)^2) + ...
		  phiM(2)*c1m*(Gm^2-1)*exp(c2m*(Gm^2-1)^2)*Gm^2 + ...
		  phiM(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
	%
	stA = phiA(1)*ce*(Geh(4)^2-Geh(2)^2) + ...
		  phiA(2)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
		  phiA(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
	%
	st = (stM*hM(i)+stA*hA(i))/(hM(i)+hA(i));	% mean circum. stress
	%
	szM = phiM(1)*ce*(Geh(5)^2-Geh(1)^2) + ...
		  phiM(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
		  phiM(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
	%
	szA = phiA(1)*ce*(Geh(6)^2-Geh(2)^2) + ...
		  phiA(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
		  phiA(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
	%
	sz = (szM*hM(i)*(2*ri(i)+hM(i))+szA*hA(i)*(2*ri(i)+2*hM(i)+hA(i)))/((hM(i)+hA(i))*(2*ri(i)+hM(i)+hA(i)));	% mean axial stress
	%
	ctM = phiM(1)*2*ce*Geh(3)^2 + ...
		  phiM(2)*2*c1m*Gm^2*exp(c2m*(Gm^2-1)^2)*((Gm^2-1)+Gm^2*(1+2*c2m*(Gm^2-1)^2)) + ...
	      phiM(4)*2*c1c*Gc^2*exp(c2c*(Gc^2-1)^2)*((Gc^2-1)+Gc^2*(1+2*c2c*(Gc^2-1)^2)*sin(alp)^2)*sin(alp)^2;
	%
	ctA = phiA(1)*2*ce*Geh(4)^2 + ...
		  phiA(2)*2*c1c*Gc^2*exp(c2c*(Gc^2-1)^2)*((Gc^2-1)+Gc^2*(1+2*c2c*(Gc^2-1)^2)) + ...
	      phiA(4)*2*c1c*Gc^2*exp(c2c*(Gc^2-1)^2)*((Gc^2-1)+Gc^2*(1+2*c2c*(Gc^2-1)^2)*sin(alp)^2)*sin(alp)^2;
	%
	ct = (ctM*hM(i)+ctA*hA(i))/(hM(i)+hA(i));	% mean circum. stiffness
	%
	czM = phiM(1)*2*ce*Geh(5)^2 + ...
		  phiM(3)*2*c1c*Gc^2*exp(c2c*(Gc^2-1)^2)*((Gc^2-1)+Gc^2*(1+2*c2c*(Gc^2-1)^2)) + ...
	      phiM(4)*2*c1c*Gc^2*exp(c2c*(Gc^2-1)^2)*((Gc^2-1)+Gc^2*(1+2*c2c*(Gc^2-1)^2)*cos(alp)^2)*cos(alp)^2;
	%
	czA = phiA(1)*2*ce*Geh(6)^2 + ...
		  phiA(3)*2*c1c*Gc^2*exp(c2c*(Gc^2-1)^2)*((Gc^2-1)+Gc^2*(1+2*c2c*(Gc^2-1)^2)) + ...
	      phiA(4)*2*c1c*Gc^2*exp(c2c*(Gc^2-1)^2)*((Gc^2-1)+Gc^2*(1+2*c2c*(Gc^2-1)^2)*cos(alp)^2)*cos(alp)^2;
	%
	cz = (czM*hM(i)*(2*ri(i)+hM(i))+czA*hA(i)*(2*ri(i)+2*hM(i)+hA(i)))/((hM(i)+hA(i))*(2*ri(i)+hM(i)+hA(i)));	% mean axial stiffness
	%
	subplot(336)
	hold on
	bar((i-1)/(sts-1)*100-2,ct/1000,3,'FaceColor',[0.2 0.2 0.2])
	bar((i-1)/(sts-1)*100+2,cz/1000,3,'FaceColor',[0.7 0.7 0.7])
	ylabel 'c_{ijkl} [MPa]'
	set(gca,'fontsize',11,'xlim',[-7 107],'xtick',s/s(end)*100,'xticklabel',num2str(s))
	legend('c_{\theta\theta\theta\theta}','c_{zzzz}','location','northwest')
	%
	WM = phiM(1)*ce/2*(Geh(1)^2+Geh(3)^2+Geh(5)^2-3) + ...
		 phiM(2)*c1m/4/c2m*(exp(c2m*(Gm^2-1)^2)-1) + ...
		 phiM(3)*c1c/4/c2c*(exp(c2c*(Gc^2-1)^2)-1) + ...
	     phiM(4)*c1c/4/c2c*(exp(c2c*(Gc^2-1)^2)-1); % medial energy per unit current medial volume
	%
	WA = phiA(1)*ce/2*(Geh(2)^2+Geh(4)^2+Geh(6)^2-3) + ...
		 phiA(2)*c1c/4/c2c*(exp(c2c*(Gc^2-1)^2)-1) + ...
		 phiA(3)*c1c/4/c2c*(exp(c2c*(Gc^2-1)^2)-1) + ...
	     phiA(4)*c1c/4/c2c*(exp(c2c*(Gc^2-1)^2)-1); % adv. energy per unit current adv. volume
	%
    %** mean total stored energy per current (W) or reference (WR) unit volume
    %
	W  = (JM*(2*rio+hMo)*hMo*WM + JA*(2*rio+2*hMo+hAo)*hAo*WA)/(JM*(2*rio+hMo)*hMo + JA*(2*rio+2*hMo+hAo)*hAo);
	WR = (JM*(2*rio+hMo)*hMo*WM + JA*(2*rio+2*hMo+hAo)*hAo*WA)/((2*rio+hMo+hAo)*(hMo+hAo));
	%
	subplot(339)
	hold on
	bar((i-1)/(sts-1)*100-2,W,3,'FaceColor',[0.2 0.2 0.2])
	bar((i-1)/(sts-1)*100+2,WR,3,'FaceColor',[0.7 0.7 0.7])
	xlabel 'G&R time [days]'
	ylabel 'Stored energy [kJ/m^3]'
	set(gca,'fontsize',11,'xlim',[-7 107],'xtick',s/s(end)*100,'xticklabel',num2str(s))
	legend('W','W_R','location','northwest')
	%
	drawnow
	%
end
%
subplot(231)
legend(num2str(s),'location','best')
%
%% Plot results
%
ltM = (2*ri+hM)/(2*rio+hMo);				%* computed medial circumferential stretch
ltA = (2*ri+2*hM+hA)/(2*rio+2*hMo+hAo);		%* computed adventitial circumferential stretch
lrM = hM/hMo;								%* computed medial radial stretch
lrA = hA/hAo;								%* computed adventitial radial stretch
JM  = ltM'.*lrM'*lz;						%* computed medial Jacobian (change of medial mass)
JA  = ltA'.*lrA'*lz;						%* computed adventitial Jacobian (change of adventitial mass)
%
figure
set(gcf,'position',[0.3*scrsz(3) 0.25*scrsz(4) 0.6*scrsz(3) 0.5*scrsz(4)])
%
%** luminal radius
%
subplot(241)
hold on
plot(s,ri/rio,'linew',1.5)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('$a/a_o$','interpreter','latex')
set(gca,'fontsize',11)
%
%** medial wall thickness
%
subplot(242)
hold on
plot(s,hM/hMo,'linew',1.5)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('$h_M/h_{Mo}$','interpreter','latex')
set(gca,'fontsize',11)
%
%** adventitial wall thickness
%
subplot(243)
hold on
plot(s,hA/hAo,'linew',1.5)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('$h_A/h_{Ao}$','interpreter','latex')
set(gca,'fontsize',11)
%
%** arterial wall thickness
%
subplot(244)
hold on
plot(s,(hM+hA)/(hMo+hAo),'linew',1.5)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('$h/h_o$','interpreter','latex')
set(gca,'fontsize',11)
%
%** referential mass density of medial smc relative to homeostatic
%
subplot(245)
hold on
plot(s,JM.*emc(:,2)/emco(2),'linew',1.5)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('G\&R time [days]','interpreter','latex')
ylabel('$\rho^m_{MR}/\rho^m_{Mo}$','interpreter','latex')
set(gca,'fontsize',11)
%
%** referential mass density of adventitial collagen relative to homeostatic
%
subplot(246)
hold on
plot(s,JA.*emc(:,5)/emco(5),'linew',1.5)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('G\&R time [days]','interpreter','latex')
ylabel('$\rho^c_{AR}/\rho^c_{Ao}$','interpreter','latex')
set(gca,'fontsize',11)
%
%** stimulus function Upsilon
%
subplot(247)
hold on
plot(s,ones(size(s)),'linew',1.5)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('G\&R time [days]','interpreter','latex')
ylabel('$\Upsilon^c$','interpreter','latex')
set(gca,'fontsize',11)
%
%** inner pressure & inflammation
%
subplot(248)
hold on
plot(s,P/Po,'k','linew',1.5)
plot(s,inf,'r','linew',1.5)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('G\&R time [days]','interpreter','latex')
hl = legend('$P/P_o$','$\varrho_f/\varrho_{fm}$');
set(hl,'interpreter','latex','location','southeast')
set(gca,'fontsize',11)
%
%** save computed values for further computations with full model
%
save('DTA_cmm.mat',...
	 'PAR','parT0','parC0','parT4','parC4',...
	 'emco','riio','so','PPfo','rotf0','lzivo','rolz0',...
	 'emch','riih','sh','PPfh','rotf4','lzivh','rolz4',...
	 'KicM','KscM','KfcM','etaUm','etaUc','etaq','fctTC')
%