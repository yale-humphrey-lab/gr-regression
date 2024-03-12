function PL = PreForGRT(parT,dl,PAR,emc1,rii1,NI1)
% function PL = PreForGRT(parT,dl,PAR,emc1,rii1,emc2,rii2,NI1)
%
%**	just call PreForGRTC to obtain values of P, L with fibers under tension
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
parC = parT(1:4);	% tension parameters for fibers under compression, if any
%
PL = PreForGRTC(parC,dl,PAR,parT,emc1,rii1,NI1);
% PL = PreForGRTC(parC,dl,PAR,parT,emc1,rii1,emc2,rii2,NI1);
%
end
