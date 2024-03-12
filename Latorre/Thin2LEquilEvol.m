function reseqs = Thin2LEquilEvol(x,load,PAR,parT,emco,riio,sgmIo,Keta,actpar)
%
%**** solves mechanobiollogically equilibrated equations for a
%     bilayered uniform (per layer) model of the arterial wall
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
%** unknowns
%
ri    = x(1);					% inner radius
hM    = x(2);					% medial thickness
hA    = x(3);					% adventitial thickness
phicM = x(4);					% medial collagen mass fraction
f     = x(5);					% vessel axial force
%
%** load
%
P   = load(1);					% inner pressure
Q   = load(2);					% normalized flow rate Q/Qo
lz  = load(3);					% axial stretch
inf = load(4);					% inflammation
%
%** PAR
%
c   = PAR(1);					% c elastin
Get = PAR(2);					% circumferential deposition stretch elastin
Gez = PAR(3);					% axial deposition stretch elastin
Bt  = PAR(4);					% fraction of circumferential collagen within the adventitia
Bz  = PAR(5);					% fraction of axial collagen within the adventitia
alp = PAR(6);					% orientation of diagonal collagen wrt axial direction
%
betaM = [Bz 1-Bz];				% medial betas [bzM 2*bdM]
betaA = [Bt Bz 1-Bt-Bz];		% adventitial betas [btA bzA 2*bdA]
%
%** parT
%
c1m = parT(1);					% c1 muscle (tension)
c2m = parT(2);					% c2 muscle (tension)
c1c = parT(3);					% c1 collagen (tension)
c2c = parT(4);					% c2 collagen (tension)
Gm  = parT(5);					% medial circumferential deposition stretch (combined medial collagen and smc)
Gc  = parT(6);					% deposition stretch (collagen)
%
%** emco and riio
%
phieMo = emco(1);				% medial elastin mass fraction
phimMo = emco(2);				% medial smc mass fraction
phicMo = emco(3);				% medial collagen mass fraction
phieAo = emco(4);				% adventitial elastin mass fraction
phicAo = emco(5);				% adventitial collagen mass fraction
%
rio  = riio(1);					% inner radius at o
rMAo = riio(2);					% M-A radius at o
roo  = riio(3);					% outer radius at o
%
hMo = rMAo-rio;					% medial thickness
hAo = roo-rMAo;					% adventitial thickness
%
%** K and eta's
%
KicM  = Keta(1); KscM  = Keta(2); KfcM  = Keta(3);	% gain parameters
etaUm = Keta(4); etaUc = Keta(5); etaq  = Keta(6);	% eta's
%
%** active parameters
%
Tmax = actpar(1); flact = actpar(2); CB = actpar(3); CS = actpar(4);
%
C = CB - CS*(Q/(ri/rio)^3-1);				% vaso-constrictor/dilator ratio
%
ltMoh = (2*ri+hM)/(2*rio+hMo);				% incremental circum. stretch M
ltAoh = (2*ri+2*hM+hA)/(2*rio+2*hMo+hAo);	% incremental circum. stretch A
lrMoh = hM/hMo;								% incremental radial stretch M
lrAoh = hA/hAo;								% incremental radial stretch A
lzoh  = lz;									% incremental axial stretch M&A
%
JM = ltMoh*lrMoh*lzoh;						% volume ratio media
JA = ltAoh*lrAoh*lzoh;						% volume ratio adventitia
%
phimM = phimMo/JM*(JM*phicM/phicMo)^(etaUm*etaq);	% medial smc mass fraction (mass growth)
phicA = phicAo/JA*(JM*phicM/phicMo)^(etaUc*  1 );	% adventitial collagen mass fraction (mass growth)
%
phieM = phieMo/JM;							% medial elastin mass fraction (no growth)
phieA = phieAo/JA;							% adventitial elastin mass fraction (no growth)
%
Geh = [1/Get/Gez*lrMoh 1/Get/Gez*lrAoh Get*ltMoh Get*ltAoh Gez*lzoh Gez*lzoh];	%* Ge*Foh = [GerMh GerAh GetMh GetAh GezMh GezAh]
%
phiM = [phieM phimM phicM*betaM];			% local mass fractions of medial [e mt cz 2*cd]
phiA = [phieA phicA*betaA];					% local mass fractions of adventitial [e ct cz 2*cd]
%
%** circumferential stress in media and adventitia
%
stM = phiM(1)*c*(Geh(3)^2-Geh(1)^2) + ...
	  phiM(2)*c1m*(Gm^2-1)*exp(c2m*(Gm^2-1)^2)*Gm^2 + ...
	  phiM(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2 + ...
	  phiM(2)*Tmax*1*(1-exp(-C^2))*flact;			%* 1 = a/a;
%
stA = phiA(1)*c*(Geh(4)^2-Geh(2)^2) + ...
	  phiA(2)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	  phiA(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
%
%** axial stress in media and adventitia
%
szM = phiM(1)*c*(Geh(5)^2-Geh(1)^2) + ...
	  phiM(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	  phiM(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
szA = phiA(1)*c*(Geh(6)^2-Geh(2)^2) + ...
	  phiA(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	  phiA(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
%** mechanobiollogically equilibrated equations in residual form R -> 0
%
reseqs = [ KicM*((P*ri/(hM+hA)+f/(pi*(hM+hA)*(2*ri+hM+hA)))/sgmIo-1) - KscM*(Q/(ri/rio)^3-1) + KfcM*inf	% Ups = 1
	  
		   phieM+phimM+phicM-1											% phiM = 1
	       phieA+phicA-1												% phiA = 1

	       stM*hM+stA*hA - P*ri											% circumferential equilibrium equation
	       szM*pi*hM*(2*ri+hM) + szA*pi*hA*(2*ri+2*hM+hA) - f ];		% axial equilibrium equation
%
end