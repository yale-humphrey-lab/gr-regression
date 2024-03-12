function PL = PreForGRTC(parC,dl,PAR,parT,emc1,rii1,N1)
% function PL = PreForGRTC(parC,dl,PAR,parT,emc1,rii1,emc2,rii2,N1)
%
%**	Function to determine Pressure (P) and Transducer axial force (L)
%   for all [outer diameter, axial stretch] included in [dl] using
%   a CMM-based bilayered 4-fiber thin-wall model of the arterial wall
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
%** parC
%
c1mc = parC(1);						% c1c muscle
c2mc = parC(2);						% c2c muscle
c1cc = parC(3);						% c1c collagen
c2cc = parC(4);						% c1c collagen
%
%** parT
%
c1mt = parT(1);						% c1t muscle
c2mt = parT(2);						% c2t muscle
c1ct = parT(3);						% c1t collagen
c2ct = parT(4);						% c2t collagen
Gm   = parT(5);						% circumferential deposition stretch (combined medial collagen and smc)
Gc   = parT(6);						% deposition stretch (collagen)
%
c    = parT(7);						% c elastin
%
%** PAR
%
% c   = PAR(1);						% c elastin
%
Get = PAR(1);						% circumferential deposition stretch elastin
Gez = PAR(2);						% axial deposition stretch elastin
Bt  = PAR(3);						% fraction of circumferential collagen within the adventitia
Bz  = PAR(4);						% fraction of axial collagen within the adventitia
alp = PAR(5);						% orientation of diagonal collagen wrt axial direction
%
betaM = [Bz 1-Bz];					% medial betas [bzM 2*bdM]
betaA = [Bt Bz 1-Bt-Bz];			% adventitial betas [btA bzA 2*bdA]
%
%** identify radii at original homeostatic state
%
if N1 <= size(dl,1)					% case setpar = 1, see <PreForGR.m>
	riio = rii1;
% else								% case setpar = 2, see <PreForGR.m>
% 	riio = rii2;
end
%
rio  = riio(1);						% inner radius at homeostatic configuration o
rMAo = riio(2);						% M-A radius at homeostatic configuration o
roo  = riio(3);						% outer radius at homeostatic configuration o
%
hMo = rMAo-rio;						% medial thickness at o
hAo = roo-rMAo;						% adventitial thickness at o
%
%---------------------  first homeostatic state [1]  ----------------------
%
%** dl
%
ro = dl(1:N1,1)/2;					% outer radius
lz = dl(1:N1,2);					% axial stretch from homeostatic configuration 1
%
%** emco and riio
%
phiM = [emc1(1:2) emc1(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
phiA = [emc1(4)   emc1(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
ri1  = rii1(1);						% inner radius at homeostatic configuration 1
rMA1 = rii1(2);						% M-A radius at homeostatic configuration 1
ro1  = rii1(3);						% outer radius at homeostatic configuration 1
%
hM1 = rMA1-ri1;						% medial thickness at 1
hA1 = ro1-rMA1;						% adventitial thickness at 1
%
ri  = sqrt(ro.^2+1./lz*(ri1^2-ro1^2));		% inner radii at P-d test
rMA = sqrt(ro.^2+1./lz*(rMA1^2-ro1^2));		% M-A radii at P-d test
%
hM = rMA-ri;						% medial thicknesses at P-d test
hA = ro-rMA;						% adventitial thicknesses at P-d test
%
SM = pi./lz*(rMA1^2-ri1^2);			% medial cross-sectional area at P-d test
SA = pi./lz*(ro1^2-rMA1^2);			% adventitial cross-sectional area at P-d test
%
ltM = (2*ri+hM)/(2*ri1+hM1);		% circumferential stretch media
ltA = (2*ro-hA)/(2*ro1-hA1);		% circumferential stretch adventitia
%
lrM = 1./ltM./lz;					% radial stretch media (approx. hM/hMh)
lrA = 1./ltA./lz;					% radial stretch adventitia (approx. hA/hAh)
%
ldM = sqrt(ltM.^2*sin(alp)^2+lz.^2*cos(alp)^2);		% diagonal stretch media
ldA = sqrt(ltA.^2*sin(alp)^2+lz.^2*cos(alp)^2);		% diagonal stretch adventitia
%
lzo1 = 1;							% incremental axial stretch between o and 1, media and adventitia
%
lrMo1 = hM1/hMo;					% incremental radial stretch between o and 1, media
lrAo1 = hA1/hAo;					% incremental radial stretch between o and 1, adventitia
%
ltMo1 = (2*ri1+hM1)/(2*rio+hMo);	% incremental circum. stretch between o and 1, media
ltAo1 = (2*ro1-hA1)/(2*roo-hAo);	% incremental circum. stretch between o and 1, adventitia
%
Ge = [1/Get/Gez*lrMo1 1/Get/Gez*lrAo1 Get*ltMo1 Get*ltAo1 Gez*lzo1 Gez*lzo1];	%* Fo1*Ge
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
PL = [ (stM.*hM+stA.*hA)./ri	szM.*SM+szA.*SA-pi*ri.*(stM.*hM+stA.*hA) ]; % P and L
%
% if N1 < size(dl,1)     % case setpar = 1, see <PreForGR.m>
% 	%
% 	%----------------  second homeostatic state [2 = h]  ------------------
% 	%
% 	%** dl
% 	%
% 	ro = dl(N1+1:end,1)/2;				% outer radius
% 	lz = dl(N1+1:end,2);				% axial stretch from homeostatic configuration 2 = h
% 	%
% 	%** emch and riih
% 	%
% 	phiM = [emc2(1:2) emc2(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
% 	phiA = [emc2(4)   emc2(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
% 	%
% 	ri2  = rii2(1);						% inner radius at homeostatic configuration 2 = h
% 	rMA2 = rii2(2);						% M-A radius at homeostatic configuration 2 = h
% 	ro2  = rii2(3);						% outer radius at homeostatic configuration 2 = h
% 	%
% 	hM2 = rMA2-ri2;						% medial thickness at 2 = h
% 	hA2 = ro2-rMA2;						% adventitial thickness at 2 = h
% 	%
% 	ri  = sqrt(ro.^2+1./lz*(ri2^2-ro2^2));		% inner radii at P-d test
% 	rMA = sqrt(ro.^2+1./lz*(rMA2^2-ro2^2));		% M-A radii at P-d test
% 	%
% 	hM = rMA-ri;						% medial thicknesses at P-d test
% 	hA = ro-rMA;						% adventitial thicknesses at P-d test
% 	%
% 	SM = pi./lz*(rMA2^2-ri2^2);			% medial cross-sectional area at P-d test
% 	SA = pi./lz*(ro2^2-rMA2^2);			% adventitial cross-sectional area at P-d test
% 	%
% 	ltM = (2*ri+hM)/(2*ri2+hM2);		% circumferential stretch media
% 	ltA = (2*ro-hA)/(2*ro2-hA2);		% circumferential stretch adventitia
% 	%
% 	lrM = 1./ltM./lz;					% radial stretch media (approx. hM/hMh)
% 	lrA = 1./ltA./lz;					% radial stretch adventitia (approx. hA/hAh)
% 	%
% 	ldM = sqrt(ltM.^2*sin(alp)^2+lz.^2*cos(alp)^2);		% diagonal stretch media
% 	ldA = sqrt(ltA.^2*sin(alp)^2+lz.^2*cos(alp)^2);		% diagonal stretch adventitia
% 	%
% 	lzo2 = 1;							% incremental axial stretch between o and 2 = h, media and adventitia
% 	%
% 	lrMo2 = hM2/hMo;					% incremental radial stretch between o and 2 = h, media
% 	lrAo2 = hA2/hAo;					% incremental radial stretch between o and 2 = h, adventitia
% 	%
% 	ltMo2 = (2*ri2+hM2)/(2*rio+hMo);	% incremental circum. stretch between o and 2 = h, media
% 	ltAo2 = (2*ro2-hA2)/(2*roo-hAo);	% incremental circum. stretch between o and 2 = h, adventitia
% 	%
% 	Ge = [1/Get/Gez*lrMo2 1/Get/Gez*lrAo2 Get*ltMo2 Get*ltAo2 Gez*lzo2 Gez*lzo2];	%* Fo2*Ge = Foh*Ge
%     %
%     %** circ. stress media
% 	%
% 	c1m = c1mt*ones(size(ro));
% 	c2m = c2mt*ones(size(ro));
% 	c1m(ltM<1/Gm) = c1mc;
% 	c2m(ltM<1/Gm) = c2mc;
% 	c1cd = c1ct*ones(size(ro));
% 	c2cd = c2ct*ones(size(ro));
% 	c1cd(ldM<1/Gc) = c1cc;
% 	c2cd(ldM<1/Gc) = c2cc;
%     %
%     %** circ. stress adventitia
% 	%
% 	stM = phiM(1)*c*(Ge(3)^2*ltM.^2-Ge(1)^2*lrM.^2) + ...
% 		  phiM(2)*c1m.*(Gm^2*ltM.^2-1).*exp(c2m.*(Gm^2*ltM.^2-1).^2)*Gm^2.*ltM.^2 + ...
% 		  phiM(4)*c1cd.*(Gc^2*ldM.^2-1).*exp(c2cd.*(Gc^2*ldM.^2-1).^2)*Gc^2.*ltM.^2*sin(alp)^2;
% 	%
% 	c1c = c1ct*ones(size(ro));
% 	c2c = c2ct*ones(size(ro));
% 	c1c(ltA<1/Gc) = c1cc;
% 	c2c(ltA<1/Gc) = c2cc;
% 	c1cd = c1ct*ones(size(ro));
% 	c2cd = c2ct*ones(size(ro));
% 	c1cd(ldA<1/Gc) = c1cc;
% 	c2cd(ldA<1/Gc) = c2cc;
% 	%
% 	stA = phiA(1)*c*(Ge(4)^2*ltA.^2-Ge(2)^2*lrA.^2) + ...
% 		  phiA(2)*c1c.*(Gc^2*ltA.^2-1).*exp(c2c.*(Gc^2*ltA.^2-1).^2)*Gc^2.*ltA.^2 + ...
% 		  phiA(4)*c1cd.*(Gc^2*ldA.^2-1).*exp(c2cd.*(Gc^2*ldA.^2-1).^2)*Gc^2.*ltA.^2*sin(alp)^2;
%     %
%     %** axial stress media
% 	%
% 	lzM = lz;
% 	%
% 	c1c = c1ct*ones(size(ro));
% 	c2c = c2ct*ones(size(ro));
% 	c1c(lzM<1/Gc) = c1cc;
% 	c2c(lzM<1/Gc) = c2cc;
% 	c1cd = c1ct*ones(size(ro));
% 	c2cd = c2ct*ones(size(ro));
% 	c1cd(ldM<1/Gc) = c1cc;
% 	c2cd(ldM<1/Gc) = c2cc;
% 	%
% 	szM = phiM(1)*c*(Ge(5)^2*lzM.^2-Ge(1)^2*lrM.^2) + ...
% 		  phiM(3)*c1c.*(Gc^2*lzM.^2-1).*exp(c2c.*(Gc^2*lzM.^2-1).^2)*Gc^2.*lzM.^2 + ...
% 		  phiM(4)*c1cd.*(Gc^2*ldM.^2-1).*exp(c2cd.*(Gc^2*ldM.^2-1).^2)*Gc^2.*lzM.^2*cos(alp)^2;
%     %
%     %** axial stress adventitia
% 	%
% 	lzA = lz;
% 	%
% 	c1c = c1ct*ones(size(ro));
% 	c2c = c2ct*ones(size(ro));
% 	c1c(lzA<1/Gc) = c1cc;
% 	c2c(lzA<1/Gc) = c2cc;
% 	c1cd = c1ct*ones(size(ro));
% 	c2cd = c2ct*ones(size(ro));
% 	c1cd(ldA<1/Gc) = c1cc;
% 	c2cd(ldA<1/Gc) = c2cc;
% 	%
% 	szA = phiA(1)*c*(Ge(6)^2*lzA.^2-Ge(2)^2*lrA.^2) + ...
% 		  phiA(3)*c1c.*(Gc^2*lzA.^2-1).*exp(c2c.*(Gc^2*lzA.^2-1).^2)*Gc^2.*lzA.^2 + ...
% 		  phiA(4)*c1cd.*(Gc^2*ldA.^2-1).*exp(c2cd.*(Gc^2*ldA.^2-1).^2)*Gc^2.*lzA.^2*cos(alp)^2;
% 	%
% 	PL = [ PL
% 		   (stM.*hM+stA.*hA)./ri	szM.*SM+szA.*SA-pi*ri.*(stM.*hM+stA.*hA) ]; % P and L
% 	%
% end
%
end