function ld = Meanld(par,LP,input)
%
%**	it generates stretch-diameter data during a Ll test for given:
%
%	- mean mechanical properties (par)
%	- pairs of Force-Pressure measurements (LP)
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
%** LP: transducer Force and Pressure during the Ll test
%
L = LP(:,1);		% transducer Force
P = LP(:,2);		% luminal pressure
%
%** input: some auxiliary parameters
%
ritf = input(1);	% inner radius at traction-free configuration
rotf = input(2);	% outer radius at traction-free configuration
%
%** current geometry is unknown at each deformation state
%   compute outer diameter from Laplace equilibrium equation
%
ro = ones(size(L));				% outer radius
lz = ones(size(L));				% axial stretch
exitflag = ones(size(L,1),2);	% convergence flag
funvalue = ones(size(L,1),2);	% residue at exit
%
r0 = rotf;			% initial radius guess for first system resolution
l0 = 1;				% initial stretch guess for first system resolution
%
inputEq = [c c1t c2t c1z c2z c1d c2d alp ritf rotf];	% input parameters for equilibrium equations
%
for i = 1:length(L)
	%
	if (i>1)
		r0 = ro(i-1);		% radius guess from previous computed value
		l0 = lz(i-1);		% stretch guess from previous computed value
	end
	%
	%* solve system of non-linear equilibrium equations for given L and P -> ro and lz
	%
	[rolz,fval,eF] = fsolve(@(rl) [P(i);L(i)]-LaplaceAxialMean(rl,inputEq),[r0,l0],optimset('Display','none'));
	exitflag(i,:) = eF;
	funvalue(i,:) = fval;
	%
	ro(i) = rolz(1);		% store solution radius
	lz(i) = rolz(2);		% store solution stretch
	%
end
%
%** axial stretch and outer diameter during the Ll test
%
ld(:,1) = lz;		% predicted axial stretch
ld(:,2) = 2*ro;		% predicted outer diameter
%
end
