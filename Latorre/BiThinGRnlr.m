function [parT,parC] = BiThinGRnlr(PAR,parTi,emc1,rii1,emc2,rii2,dlexp,PLexp,N1)
%
%**	determines SMC and collagen parameters in a CMM-based bilayered
%   4-fiber thin-wall model of the arterial wall from biaxial data
%   included in dlexp (geometry) and PLexp (loads)
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

global Gmmax Gcmax Gmmin Gcmin
%
%** determine model parameters [c1m,c2m,c1c,c2c,Gm,Gc]
%
parT = parTi;						% initial guess
%
lb = [0.0 0.0 0.0 0.0 Gmmin Gcmin];
ub = [Inf Inf Inf Inf Gmmax Gcmax];
%
%---------------------  first homeostatic state [1]  ----------------------
%
%** dlexp
%
ro  = dlexp(1:N1,1)/2;				% outer radius
lz1 = dlexp(1:N1,2);				% axial stretch from homeostatic configuration 1
%
%** rii1
%
ri1  = rii1(1);						% inner radius at homeostatic configuration 1
rMA1 = rii1(2);						% M-A radius at homeostatic configuration 1
ro1  = rii1(3);						% outer radius at homeostatic configuration 1
%
hM1 = rMA1-ri1;						% medial thickness at 1
hA1 = ro1-rMA1;						% adventitial thickness at 1
%
ri  = sqrt(ro.^2+1./lz1*(ri1^2-ro1^2));		% inner radius
rMA = sqrt(ro.^2+1./lz1*(rMA1^2-ro1^2));	% M-A radii at P-d test
%
hM = rMA-ri;						% medial thicknesses at P-d test
hA = ro-rMA;						% adventitial thicknesses at P-d test
%
ltM1 = (2*ri+hM)/(2*ri1+hM1);		% circumferential stretch media
ltA1 = (2*ro-hA)/(2*ro1-hA1);		% circumferential stretch adventitia
%
ldM1 = sqrt(ltM1.^2*sin(PAR(6))^2+lz1.^2*cos(PAR(6))^2);	% diagonal stretch media, with PAR(6) = alpha
ldA1 = sqrt(ltA1.^2*sin(PAR(6))^2+lz1.^2*cos(PAR(6))^2);	% diagonal stretch adventitia, with PAR(6) = alpha
%
I1 = find(ltM1>1/parT(5) & ldM1>1/parT(6) & ltA1>1/parT(6) & ldA1>1/parT(6) & lz1>1/parT(6));	% parT(5:6) = [Gm,Gt]
I  = I1;
%
if N1 < size(dlexp,1)    % case setpar = 1, see <PreForGR.m>
	%
	%------------------  second homeostatic state [2]  --------------------
	%
	%** dlexp
	%
	ro  = dlexp(N1+1:end,1)/2;			% outer radius
	lz2 = dlexp(N1+1:end,2);			% axial stretch from homeostatic configuration 2 = h
	%
	%** rii2
	%
	ri2  = rii2(1);						% inner radius at homeostatic configuration 2 = h
	rMA2 = rii2(2);						% M-A radius at homeostatic configuration 2 = h
	ro2  = rii2(3);						% outer radius at homeostatic configuration 2 = h
	%
	hM2 = rMA2-ri2;						% medial thickness at 2 = h
	hA2 = ro2-rMA2;						% adventitial thickness at 2 = h
	%
	ri  = sqrt(ro.^2+1./lz2*(ri2^2-ro2^2));		% inner radius
	rMA = sqrt(ro.^2+1./lz2*(rMA2^2-ro2^2));	% M-A radii at P-d test
	%
	hM = rMA-ri;						% medial thicknesses at P-d test
	hA = ro-rMA;						% adventitial thicknesses at P-d test
	%
	ltM2 = (2*ri+hM)/(2*ri2+hM2);		% circumferential stretch media
	ltA2 = (2*ro-hA)/(2*ro2-hA2);		% circumferential stretch adventitia
	%
	ldM2 = sqrt(ltM2.^2*sin(PAR(6))^2+lz2.^2*cos(PAR(6))^2);	% diagonal stretch media, with PAR(6) = alpha
	ldA2 = sqrt(ltA2.^2*sin(PAR(6))^2+lz2.^2*cos(PAR(6))^2);	% diagonal stretch adventitia, with PAR(6) = alpha
	%
	I2 = find(ltM2>1/parT(5) & ldM2>1/parT(6) & ltA2>1/parT(6) & ldA2>1/parT(6) & lz2>1/parT(6));	% parT(5:6) = [Gm,Gt]
	I  = [I1;N1+I2];
	%
end
%
%------------------------------  REGRESSION  ------------------------------
%
options = optimoptions(@lsqcurvefit,'Display','None','FunctionTolerance',1e-6,'StepTolerance',1e-6,'OptimalityTolerance',1e-6);
%
%** determine model parameters parT = [c1mt,c2mt,c1ct,c2ct,Gm,Gc]
%
II = zeros(size(I));
%
iterations = 0;
%
while ([I - II ~= zeros(size(I)); length(I) > 2])
	%
	II = I;
	%
	iterations = iterations + 1;
	%
	NI1 = length(I1);					% tension states at [1]
	%
	parT = lsqcurvefit(@(parT,dl) PreForGRT(parT,dl,PAR,emc1,rii1,emc2,rii2,NI1),parT,dlexp(I,:),PLexp(I,:),lb,ub,options);
	%
	I1 = find(ltM1>1/parT(5) & ldM1>1/parT(6) & ltA1>1/parT(6) & ldA1>1/parT(6) & lz1>1/parT(6));		% update I1, with parT(5:6) = [Gm,Gt]
	I  = I1;
	%
	if N1 < size(dlexp,1)				% case setpar = 1, see <PreForGR.m>
		%
		I2 = find(ltM2>1/parT(5) & ldM2>1/parT(6) & ltA2>1/parT(6) & ldA2>1/parT(6) & lz2>1/parT(6));	% update I2, with parT(5:6) = [Gm,Gt]
		I  = [I1;N1+I2];
		%
	end
	%
	if (length(II) ~= length(I))
		II = zeros(size(I));
	end
	%
	if (iterations == 10), break, end
	%
end
%
%** determine model parameters parC = [c1mc,c2mc,c1cc,c2cc]
%
parC = parT(1:4);
%
lb = [0 0 0 0];
ub = parT(1:4);
%
parC = lsqcurvefit(@(parC,dl) PreForGRTC(parC,dl,PAR,parT,emc1,rii1,emc2,rii2,N1),parC,dlexp,PLexp,lb,ub,options);
%
end