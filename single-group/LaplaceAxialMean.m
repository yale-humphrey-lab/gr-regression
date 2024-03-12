function PL = LaplaceAxialMean(rolz,input)
%
%**	it computes the distensional Pressure (P) and transducer
%	Force (L) from both Laplace and axial equilibrium
%	equations for given:
%
%	- current outer radius (rolz(1))
%	- axial stretch from tf configuration (rolz(2))
%	- mean mechanical properties and tf geometry (input)
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
%** retrieve outer radius and axial stretch
%
ro = rolz(1);
lz = rolz(2);
%
%** retrieve parameters in 'input'
%
c   = input(1);		% c elastin
c1t = input(2);		% c1 circumferential
c2t = input(3);		% c2 circumferential
c1z = input(4);		% c1 axial
c2z = input(5);		% c2 axial
c1d = input(6);		% c1 diagonal
c2d = input(7);		% c2 diagonal
alp = input(8);		% orientation of diagonal collagen wrt axial direction (in rad)
%
ritf = input(9);	% inner radius at traction-free configuration
rotf = input(10);	% outer radius at traction-free configuration
%
htf = rotf-ritf;	% thickness at traction-free configuration
%
%** current geometry is known from incompressibility
%
ri = sqrt(ro^2+1/lz*(ritf^2-rotf^2));		% inner radius
h  = ro-ri;									% thickness
S  = pi./lz*(rotf^2-ritf^2);				% cross-sectional area
%
%** stretches
%
lt = (2*ri+h)/(2*ritf+htf);						% circumferential stretch
lr = 1/lt/lz;									% radial stretch (approx. h/htf)
ld = sqrt(lt^2*sin(alp)^2+lz^2*cos(alp)^2);		% diagonal stretch
%
%** circumferential Cauchy stress
%
st =   c*(lt^2-lr^2) + ...								% elastin + Lagrange multiplier
	 c1t*(lt^2-1)*exp(c2t*(lt^2-1)^2)*lt^2 + ...		% circumferential fibers
   2*c1d*(ld^2-1)*exp(c2d*(ld^2-1)^2)*lt^2*sin(alp)^2;	% diagonal fibers
%
%** axial Cauchy stress
%
sz =   c*(lz^2-lr.^2) + ...								% elastin + Lagrange multiplier
	 c1z*(lz^2-1)*exp(c2z*(lz^2-1)^2)*lz^2 + ...		% axial fibers
   2*c1d*(ld^2-1)*exp(c2d*(ld^2-1)^2)*lz^2*cos(alp)^2;	% diagonal fibers
%
%** Pressure and transducer Force from corresponding equilibrium equations
%
PL(1,1) = st*h/ri;						% predicted luminal Pressure
PL(2,1) = sz*S - pi*ri^2*PL(1,1);		% predicted transducer axial Force
%
end
