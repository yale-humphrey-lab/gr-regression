function P = LaplaceMean(ro,lz,input)
%
%**	it computes the distensional Pressure (P) from the Laplace
%	equilibrium equation for given:
%
%	- current outer radius (ro)
%	- axial stretch from tf configuration (lz)
%	- mean mechanical properties and tf geometry (input)
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
%** retrieve parameters in 'input'
%
c   = input(1);		% c elastin
c1t = input(2);		% c1 circumferential
c2t = input(3);		% c2 circumferential
c1d = input(4);		% c1 diagonal
c2d = input(5);		% c2 diagonal
alp = input(6);		% orientation of diagonal collagen wrt axial direction (in rad)
%
ritf = input(7);	% inner radius at traction-free configuration
rotf = input(8);	% outer radius at traction-free configuration
%
htf = rotf-ritf;	% thickness at traction-free configuration
%
%** current geometry is known from incompressibility
%
ri = sqrt(ro^2+1/lz*(ritf^2-rotf^2));		% inner radius
h  = ro-ri;									% thickness
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
%** Pressure from Laplace equilibrium equation
%
P = st*h/ri;
%
end