function [Pd,Ll] = BiaxialDataMean(parBDM,inputBDM,protocols)
%
%**	Function to generate biaxial data from arterial MEAN mechanical
%	properties
%
%	See, for example, Table S2 in:
%
%	[1] M Bersi et al. Excessive Adventitial Remodeling Leads
%		to Early Aortic Maladaptation in Angiotensin-Induced
%		Hypertension. Hypertension, 67(5), 890-896, 2016.
%
%  ------------  marcos.latorre@yale.edu (2017)  ------------

%
%** inputBDM: auxiliary parameters
%
ritf = inputBDM(1);				% inner radius at traction-free configuration
rotf = inputBDM(2);				% outer radius at traction-free configuration
lziv = inputBDM(3);				% in-vivo axial stretch (from tf)
%
inputMean(1) = ritf;			% inner radius at tf, used in <MeandL> and <Meanld> below
inputMean(2) = rotf;			% outer radius at tf, used in <MeandL> and <Meanld> below
%
mmHg_to_kPa = 0.13332;			% conversion factor: mmHg -> kPa
%
% figure						% get figure ready
% subplot(221), hold on, grid on, xlabel 'Outer Diameter [mm]', ylabel 'Pressure [mmHg]'
% subplot(222), hold on, grid on, xlabel 'Axial stretch', ylabel 'Outer diameter [mm]'
% subplot(223), hold on, grid on, xlabel 'Pressure [mmHg]', ylabel 'Transducer force [mN]'
% subplot(224), hold on, grid on, xlabel 'Axial stretch', ylabel 'Transducer force [mN]'
% legendPd = [];
% legendLl = [];
%
%** compute diameter-Force data during Pd test
%
Pd = [];						% P-d-L-l data from Pd tests
%
PPd = mmHg_to_kPa*(10:5:140)';	% Pressure at Pd tests (in kPa)  /  10:5:150  /  10:10:140
%
switch protocols				% protocols to consider for Pd tests
	case 1
		stretches = 1.00;
	case 2
		stretches = [1.00 1.05];
	case 3
		stretches = [0.95 1.00 1.05];
end
%
for l = stretches		% axial stretches relative to homeostatic configuration
	%
	lPd = l*lziv*ones(size(PPd));		% constant stretch at Pd test
	PlPd = [PPd lPd];					% Pressure-stretch pairs from Pd test, used in <MeandL> below
	%
	dLPd = MeandL(parBDM,PlPd,inputMean);		% diameter-Force data during Pd test
	%
	Pd = [Pd; PlPd(:,1) dLPd(:,1) dLPd(:,2) PlPd(:,2)];		% store P-d-L-l data from Pd test
	%
	% 	subplot(221)
	% 	plot(dLPd(:,1),PlPd(:,1)/mmHg_to_kPa)	% plot diameter-Pressure
	% 	%
	% 	subplot(223)
	% 	plot(PlPd(:,1)/mmHg_to_kPa,dLPd(:,2))	% plot Pressure-Force
	% 	%
	% 	legendPd = [legendPd; num2str(l,'%3.2f')];		% update legend
	%
end
%
Pd = [Pd; 0 2*rotf 0 1];		% add tf state to P-d tests, comment line otherwise
%
%** compute diameter-stretch data during Ll test
%
Ll = [];						% P-d-L-l data from Ll tests
%
LLl = (2:2:38)';				% transducer Force at Ll tests
%
switch protocols				% protocols to consider for Ll tests
	case 1
		pressures = 100;
	case 2
		pressures = [100 140];
	case 3
		pressures = [60 100 140];
end
%
for Pres = pressures			% different Pressures [mmHg]
	%
	PLl = mmHg_to_kPa*Pres*ones(size(LLl));		% constant pressure during Ll test [kPa]
	LPLl = [LLl PLl];							% Force-Pressure pairs from Ll test, used in <Meanld> below
	%
	ldLl = Meanld(parBDM,LPLl,inputMean);		% stretch-diameter data during Ll test
	%
	Ll = [Ll; LPLl(:,2) ldLl(:,2) LPLl(:,1) ldLl(:,1)];		% store P-d-L-l data from Ll test
	%
	% 	subplot(222)
	% 	plot(ldLl(:,1),ldLl(:,2))	% plot stretch-diameter
	% 	%
	% 	subplot(224)
	% 	plot(ldLl(:,1),LPLl(:,1))	% plot stretch-Force
	% 	%
	% 	legendLl = [legendLl; num2str(Pres,'%3.3d')];	% update legend
	%
end
%
Ll = [Ll; 0 2*rotf 0 1];		% add tf state to L-l tests, comment line otherwise
%
% subplot(221)
% hl = legend(legendPd);
% set(hl,'Location','northwest')
% subplot(222)
% hl = legend(legendLl);
% set(hl,'Location','northwest')
% subplot(223)
% hl = legend(legendPd);
% set(hl,'Location','northwest')
% subplot(224)
% hl = legend(legendLl);
% set(hl,'Location','northwest')
%
end