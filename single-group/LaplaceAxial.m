function PF = LaplaceAxial(rolz,PAR,parT,parC,riio,emch,riih,lzoh)
%
%**	it computes the distensional Pressure (P) and transducer
%	Force (F) from both Laplace and axial equilibrium
%	equations for given:
%
%	- current outer radius (rolz(1))
%	- axial stretch from h configuration (rolz(2))
%	- mechanical properties and h geometry (input)
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
%** retrieve outer radius and axial stretch
%
ro = rolz(1);
lz = rolz(2);
%
%** PAR
%
% c   = PAR(1);					% c elastin
%
Get = PAR(1);					% circumferential deposition stretch elastin
Gez = PAR(2);					% axial deposition stretch elastin
Bt  = PAR(3);					% fraction of circumferential collagen within the adventitia
Bz  = PAR(4);					% fraction of axial collagen within the adventitia
alp = PAR(5);					% orientation of diagonal collagen wrt axial direction
%
betaM = [Bz 1-Bz];				% medial betas [bzM 2*bdM]
betaA = [Bt Bz 1-Bt-Bz];		% adventitial betas [btA bzA 2*bdA]
%
%** parT
%
c1mt = parT(1);					% c1t muscle
c2mt = parT(2);					% c2t muscle
c1ct = parT(3);					% c1t collagen
c2ct = parT(4);					% c2t collagen
Gm   = parT(5);					% circumferential deposition stretch (combined medial collagen and smc)
Gc   = parT(6);					% deposition stretch (collagen)
%
c    = parT(7);					% c elastin
%
%** parC
%
c1mc = parC(1);					% c1c muscle
c2mc = parC(2);					% c2c muscle
c1cc = parC(3);					% c1c collagen
c2cc = parC(4);					% c1c collagen
%
%** riio
%
rio  = riio(1);						% inner radius at homeostatic configuration o
rMAo = riio(2);						% M-A radius at homeostatic configuration o
roo  = riio(3);						% outer radius at homeostatic configuration o
%
hMo = rMAo-rio;						% medial thickness at o
hAo = roo-rMAo;						% adventitial thickness at o
%
%** emch and riih
%
phiM = [emch(1:2) emch(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
phiA = [emch(4)   emch(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
rih  = riih(1);						% inner radius at homeostatic configuration h
rMAh = riih(2);						% M-A radius at homeostatic configuration h
roh  = riih(3);						% outer radius at homeostatic configuration h
%
hMh = rMAh-rih;						% medial thickness at h
hAh = roh-rMAh;						% adventitial thickness at h
%
ri  = sqrt(ro.^2+1./lz*(rih^2-roh^2));	% inner radii at P-d test
rMA = sqrt(ro.^2+1./lz*(rMAh^2-roh^2));	% M-A radii at P-d test
%
hM = rMA-ri;				% medial thicknesses at P-d test
hA = ro-rMA;				% adventitial thicknesses at P-d test
%
SM = pi./lz*(rMAh^2-rih^2);	% medial cross-sectional area at P-d test
SA = pi./lz*(roh^2-rMAh^2);	% adventitial cross-sectional area at P-d test
%
ltM = (2*ri+hM)/(2*rih+hMh);	% circumferential stretch media
ltA = (2*ro-hA)/(2*roh-hAh);	% circumferential stretch adventitia
%
lrM = 1./ltM./lz;			% radial stretch media (approx. hM/hMh)
lrA = 1./ltA./lz;			% radial stretch adventitia (approx. hA/hAh)
%
ldM = sqrt(ltM.^2*sin(alp)^2+lz.^2*cos(alp)^2);	% diagonal stretch media
ldA = sqrt(ltA.^2*sin(alp)^2+lz.^2*cos(alp)^2);	% diagonal stretch adventitia
%
% lzoh;								% incremental axial stretch (input)
%
lrMoh = hMh/hMo;					% incremental radial stretch media
lrAoh = hAh/hAo;					% incremental radial stretch adventitia
%
ltMoh = (2*rih+hMh)/(2*rio+hMo);	% incremental circumferential stretch media
ltAoh = (2*roh-hAh)/(2*roo-hAo);	% incremental circumferential stretch adventitia
%
Ge = [1/Get/Gez*lrMoh 1/Get/Gez*lrAoh Get*ltMoh Get*ltAoh Gez*lzoh Gez*lzoh];	%* Foh*Ge
%
%** circ. stress media
%
c1m = c1mt*ones(size(ro));
c2m = c2mt*ones(size(ro));
c1m(ltM<1/Gm) = c1mc;
c2m(ltM<1/Gm) = c2mc;
c1cd = c1ct*ones(size(ro));
c2cd = c2ct*ones(size(ro));
c1cd(ldM<1/Gc) = c1cc;
c2cd(ldM<1/Gc) = c2cc;
%
stM = phiM(1)*c*(Ge(3)^2*ltM.^2-Ge(1)^2*lrM.^2) + ...
	  phiM(2)*c1m.*(Gm^2*ltM.^2-1).*exp(c2m.*(Gm^2*ltM.^2-1).^2)*Gm^2.*ltM.^2 + ...
	  phiM(4)*c1cd.*(Gc^2*ldM.^2-1).*exp(c2cd.*(Gc^2*ldM.^2-1).^2)*Gc^2.*ltM.^2*sin(alp)^2;
%
%** circ. stress adventitia
%
c1c = c1ct*ones(size(ro));
c2c = c2ct*ones(size(ro));
c1c(ltA<1/Gc) = c1cc;
c2c(ltA<1/Gc) = c2cc;
c1cd = c1ct*ones(size(ro));
c2cd = c2ct*ones(size(ro));
c1cd(ldA<1/Gc) = c1cc;
c2cd(ldA<1/Gc) = c2cc;
%
stA = phiA(1)*c*(Ge(4)^2*ltA.^2-Ge(2)^2*lrA.^2) + ...
	  phiA(2)*c1c.*(Gc^2*ltA.^2-1).*exp(c2c.*(Gc^2*ltA.^2-1).^2)*Gc^2.*ltA.^2 + ...
	  phiA(4)*c1cd.*(Gc^2*ldA.^2-1).*exp(c2cd.*(Gc^2*ldA.^2-1).^2)*Gc^2.*ltA.^2*sin(alp)^2;
%
%** axial stress media
%
lzM = lz;
%
c1c = c1ct*ones(size(ro));
c2c = c2ct*ones(size(ro));
c1c(lzM<1/Gc) = c1cc;
c2c(lzM<1/Gc) = c2cc;
c1cd = c1ct*ones(size(ro));
c2cd = c2ct*ones(size(ro));
c1cd(ldM<1/Gc) = c1cc;
c2cd(ldM<1/Gc) = c2cc;
%
szM = phiM(1)*c*(Ge(5)^2*lzM.^2-Ge(1)^2*lrM.^2) + ...
	  phiM(3)*c1c.*(Gc^2*lzM.^2-1).*exp(c2c.*(Gc^2*lzM.^2-1).^2)*Gc^2.*lzM.^2 + ...
	  phiM(4)*c1cd.*(Gc^2*ldM.^2-1).*exp(c2cd.*(Gc^2*ldM.^2-1).^2)*Gc^2.*lzM.^2*cos(alp)^2;
%
%** axial stress adventitia
%
lzA = lz;
%
c1c = c1ct*ones(size(ro));
c2c = c2ct*ones(size(ro));
c1c(lzA<1/Gc) = c1cc;
c2c(lzA<1/Gc) = c2cc;
c1cd = c1ct*ones(size(ro));
c2cd = c2ct*ones(size(ro));
c1cd(ldA<1/Gc) = c1cc;
c2cd(ldA<1/Gc) = c2cc;
%
szA = phiA(1)*c*(Ge(6)^2*lzA.^2-Ge(2)^2*lrA.^2) + ...
	  phiA(3)*c1c.*(Gc^2*lzA.^2-1).*exp(c2c.*(Gc^2*lzA.^2-1).^2)*Gc^2.*lzA.^2 + ...
	  phiA(4)*c1cd.*(Gc^2*ldA.^2-1).*exp(c2cd.*(Gc^2*ldA.^2-1).^2)*Gc^2.*lzA.^2*cos(alp)^2;
%
%** Pressure and transducer Force from corresponding equilibrium equations
%
PF(1,1) = (stM.*hM + stA.*hA)./ri;					% predicted luminal Pressure
PF(2,1) = szM*SM + szA*SA - pi*ri^2*PF(1,1);		% predicted transducer axial Force
%
end