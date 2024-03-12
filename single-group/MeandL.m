function dL = MeandL(par,Pl,input)
%
%**	it generates diameter-Force data from a Pd test for given:
%
%	- mean mechanical properties (par)
%	- pairs of Pressure-stretch measurements (Pl)
%	- geometry at traction-free configuration (input)
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
%** par: material parameters for arterial MEAN response
%
c   = par(1);		% c elastin
c1t = par(2);		% c1 circumferential
c2t = par(3);		% c2 circumferential
c1z = par(4);		% c1 axial
c2z = par(5);		% c2 axial
c1d = par(6);		% c1 diagonal
c2d = par(7);		% c2 diagonal
alp = par(8);		% orientation of diagonal collagen wrt axial direction (in rad)
%
%** Pl: Pressure and axial stretch during Pd test
%
P  = Pl(:,1);		% distensional Pressure
lz = Pl(:,2);		% axial stretch from configuration tf
%
%** input: some auxiliary parameters
%
ritf = input(1);	% inner radius at traction-free configuration
rotf = input(2);	% outer radius at traction-free configuration
%
htf = rotf-ritf;	% thickness at traction-free configuration
%
%** current geometry is unknown at each deformation state
%   compute outer diameter from Laplace equilibrium equation
%
ro = ones(size(P));			% outer radius
exitflag = ones(size(P));	% convergence flag
funvalue = ones(size(P));	% residue at exit
%
r0 = rotf;			% initial guess for first equation resolution
%
inputLap = [c c1t c2t c1d c2d alp ritf rotf];	% input parameters for Laplace equation
%
for i = 1:length(P)
	%
	if (i>1), r0 = ro(i-1); end		% guess radius from previous computed value
	%
	%* solve non-linear equilibrium equation for given P and lz -> ro
	%
	[ro(i),fval,eF] = fsolve(@(r) P(i)-LaplaceMean(r,lz(i),inputLap),r0,optimset('Display','none'));
	exitflag(i) = eF;
	funvalue(i) = fval;
	%
end
%
%** current geometry is known at each deformation state from incompressibility
%
ri = sqrt(ro.^2+1./lz*(ritf^2-rotf^2));		% inner radius
h  = ro-ri;									% thickness
S  = pi./lz*(rotf^2-ritf^2);				% cross-sectional area
%
%** stretches
%
lt = (2*ri+h)/(2*ritf+htf);						% circumferential stretch
lr = 1./lt./lz;									% radial stretch (approx. h/htf)
ld = sqrt(lt.^2*sin(alp)^2+lz.^2*cos(alp)^2);	% diagonal stretch
%
% stored energy
W = c/2*(lr.^2 + lt.^2 + lz.^2 - 3) + ...
	c1t/(4*c2t)*(exp(c2t*(lt.^2-1).^2) - 1) + ...
	c1z/(4*c2z)*(exp(c2z*(lz.^2-1).^2) - 1) + ...
	2*c1d/(4*c2d)*(exp(c2d*(ld.^2-1).^2) - 1);
%
%** circumferential Cauchy stress
%
st =   c*(lt.^2-lr.^2) + ...									% elastin + Lagrange multiplier
	 c1t*(lt.^2-1).*exp(c2t*(lt.^2-1).^2).*lt.^2 + ...			% circumferential fibers
   2*c1d*(ld.^2-1).*exp(c2d*(ld.^2-1).^2).*lt.^2*sin(alp)^2;	% diagonal fibers
%** axial Cauchy stress
%
sz =   c*(lz.^2-lr.^2) + ...									% elastin + Lagrange multiplier
	 c1z*(lz.^2-1).*exp(c2z*(lz.^2-1).^2).*lz.^2 + ...			% axial fibers
   2*c1d*(ld.^2-1).*exp(c2d*(ld.^2-1).^2).*lz.^2*cos(alp)^2;	% diagonal fibers
%
%** outer diameter and transducer Force during the Pd test
%
dL(:,1) = 2*ro;						% predicted outer diameter
dL(:,2) = sz.*S - pi*ri.^2.*P;		% predicted transducer axial Force
%
end
