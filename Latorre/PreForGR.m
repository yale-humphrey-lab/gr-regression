function PL = PreForGR(PAR,dl,emco,riio,emch,riih,N0,setpar)
%
%**	it distinguishes between setpar = ...
%   
%   1: smc/coll. parameters do not evolve ([parT0 parC0]  = [parT4 parC4])
%
%	2: smc/coll. parameters do     evolve ([parT0 parC0] /= [parT4 parC4])
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
global parT0 parC0 parT4 parC4
global PLexp0 PLexp4
%
switch setpar	%* sets of parameters [parT parC] to determine
	%
	case 1		%* determine [parT0 parC0] = [parT4 parC4]
		%
		%**** update parT0 and parC0 through nonlinear regression with PAR fixed
		%
		%** call BiThinGRnlr with data from tests at BOTH days 0 and 28 (4 weeks)
		%
		dlexp = dl;						%* dl includes data from BOTH tests
		PLexp = [PLexp0;PLexp4];		%* data from BOTH tests
		%
		[parT0,parC0] = BiThinGRnlr(PAR,parT0,emco,riio,emch,riih,dlexp,PLexp,N0);
		%
		%** call PreForGRTC to compute values of P and L
		%
		PL = PreForGRTC(parC0,dlexp,PAR,parT0,emco,riio,emch,riih,N0);
		%
		%** same parameters at week 4
		%
		parT4 = parT0; parC4 = parC0;
		%
	case 2		%* determine [parT0 parC0] and [parT4 parC4]
		%
		%**** update parT0 and parC0 through nonlinear regression with PAR fixed
		%
		%** call BiThinGRnlr with data ONLY from tests at day 0
		%
		dlexp0 = dl(1:N0,:);			%* dl includes data from BOTH tests
		%
		[parT0,parC0] = BiThinGRnlr(PAR,parT0,emco,riio,emco,riio,dlexp0,PLexp0,N0);
		%
		%** call PreForTC to compute values of P0 and L0
		%
		PL0 = PreForGRTC(parC0,dlexp0,PAR,parT0,emco,riio,emco,riio,N0);
		%
		%**** update parT4 and parC4 through nonlinear regression with PAR fixed
		%
		%** call BiThinGRnlr with data ONLY from tests at day 28 (4 weeks)
		%
		dlexp4 = dl(N0+1:end,:);		%* dl includes data from BOTH tests
		N4 = size(dlexp4,1);			%* P-d and L-l measurements at day 28
		%
		[parT4,parC4] = BiThinGRnlr(PAR,parT4,emch,riih,emco,riio,dlexp4,PLexp4,N4);
		%
		%** call PreForTC4 to compute values of P4 and L4
		%
		PL4 = PreForGRTC(parC4,dlexp4,PAR,parT4,emch,riih,emco,riio,N4);
		%
		PL = [PL0; PL4];
	%
end
%
end