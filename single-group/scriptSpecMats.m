
% David Li | Yale University | Feb 2024
% Regenerates synthetic biaxial data from fitted 4FF parameters
% Requires .csv file containing fitted parameters for each specimen & mass
% fractions from Movat

% clear; clc; % close all;

save_indiv = 0;								% Save .mat for individual specimens
save_group = 1;								% Save .mat for group-averaged data
plot_data = 0;								% Plot sythetic data

% mmHg_to_kPa = 0.133322;                     % Convert mmHg to kPa
% protocols = 3;								% Number of protocols for both Pd and Fl tests: 1, 2, 3
% pres      = 3;								% Experimental pressures (1 dias, 2 MAP, 3 sys)
P = [80,100,120];							% Dia/MA/sys pressures (mmHg)

% Identify group and set layer-specific mass fractions
grp = 'C1041G'; sgrp = 'all';				% Group and subgroup
SM = 0.7271;								% Med wall percentage
emcM  = [0.4857,0.4097,0.1046];				% Med fractions
ecA   = [0.0037,0.9963];					% Adv fractions

% grp = 'mgR'; sgrp = 'all';					% 'all' 'sm' 'lg'
% SM = 0.6916;								% Med wall percentage
% emcM  = [0.3201,0.5022,0.1777];				% Med fractions
% ecA   = [0.0037,0.9963];					% Adv fractions

% Import data from .csv
allspec = csvread([grp,'_ATA_allspec.csv'],1,1);
% 1       2        3       4     5        6          7
% L (mm)  OD (um)  H (um)  lziv  c (kPa)  c11 (kPa)  c21
% 8          9    10           11     12        13
% c12 (kPa)  c22  c13,4 (kPa)  c23,4  ao (deg)  Eln por

% When averaging biaxial data, separate specimens into sub-groups
if strcmp(grp,'C1041G')
	if strcmp(num,'all'); N = 1:size(allspec,1);
	else; N = str2double(num);
	end
elseif strcmp(grp,'mgR')
	% Sort by elastin porosity (averaging eln por replicates)
	if strcmp(sgrp,'sm');      N = [1 2 4 7 8];			% <41% EP mgR
	elseif strcmp(sgrp,'lg');  N = [3 5 6];				% >41% EP mgR
	% Sort by mechanical behavior (Fl)
% 	if strcmp(sgrp,'sm');      N = [2 3 7 8 4];			% Small OD mgR
% 	elseif strcmp(sgrp,'lg');  N = [1 5 6];				% Large OD mgR
	elseif strcmp(sgrp,'all'); N = 1:size(allspec,1);	% All mgR
	end
end

% Uses scripts from Latorre EqG&R code
rii_all = zeros(length(N),3);					% Inner/M-A/outer radius
Pd_all  = zeros(length(N),81,4);				% Pressure-diameter
Ll_all  = zeros(length(N),57,4);				% Force-length
ct = 1;											% Specimen count
for sp = N
	dotf = allspec(sp,2)/1e3;					% Unloaded outer diam
	htf  = allspec(sp,3)/1e3;					% Unloaded thickness
	lziv = allspec(sp,4);						% In-vivo ax stretch
	meanpar = allspec(sp,5:12);					% Biaxial parameters (Excel order)
	meanpar = meanpar([1 4:5 2:3 6:8]);			% Sort parameters by code order
	
	rotf  = dotf/2;								% Outer radius at traction-free configuration (tf) [mm]
	ritf  = rotf-htf;							% Inner radius at tf [mm]
% 	lziv  = round(lziv,4);						% In-vivo axial stretch (from tf)
	SA    = 1-SM;								% Adventitial cross-sectional area fraction from histology
	rMAtf = sqrt(SA*ritf^2+SM*rotf^2);			% M-A interface radius at tf [mm]

	meanpar(end) = deg2rad(meanpar(end));		% Bulk parameters

	Pe = P(pres)*mmHg_to_kPa;						% Experimental pressure [kPa]
	dLPd = MeandL(meanpar,[Pe lziv],[ritf rotf]);	% In-vivo outer diameter and force

	ro  = dLPd(1)/2;								% Outer radius at o [mm]
	rMA = sqrt(ro^2+1/lziv*(rMAtf^2-rotf^2));		% M-A radius at o [mm]
	ri  = sqrt(ro^2+1/lziv*(ritf^2-rotf^2));		% Inner radius at o [mm]
	rii = [ri rMA ro];								% In-vivo radii at o
	rii_all(ct,:) = rii;
	
	% Pd, Ll: pres, diam, force, ax stretch
	[Pd,Ll] = BiaxialDataMean(meanpar,[ritf rotf lziv],protocols);		% Generate biaxial testing data
	Pd_all(ct,:,:) = Pd;
	Ll_all(ct,:,:) = Ll;
	
	Pd(:,4) = Pd(:,4)/lziv;					% Reset axial stretch from in-vivo state
	Ll(:,4) = Ll(:,4)/lziv;					% Reset axial stretch from in-vivo state
	
	PLexp = [Pd(:,[1 3]);Ll(:,[1 3])];		% Experimental Pressure-Force pairs from tests Pd and Ll
	dlexp = [Pd(:,[2 4]);Ll(:,[2 4])];		% Experimental diameter-stretch pairs from tests Pd and Ll
	
	if plot_data
		figure(1);
		set(gcf,'units','normalized','outerposition',[.3 .25 .4 .73]);
		subplot(221)
			hold on; grid on
			plot(Pd(:,2),PLexp(1:size(Pd,1),1)/mmHg_to_kPa,...
				'd','markersize',3,'linew',0.75)
			xlabel 'Outer Diameter [mm]'; ylabel 'Pressure [mmHg]'
			set(gca,'xlim',[0.2 3.4],'ylim',[-10 160])
			set(gca,'fontsize',11)
		subplot(222)
			hold on; grid on
			plot(PLexp(1:size(Pd,1),1)/mmHg_to_kPa,PLexp(1:size(Pd,1),2),...
				'd','markersize',3,'linew',0.75)
			xlabel 'Pressure [mmHg]'; ylabel 'Transducer axial force [mN]'
			set(gca,'xlim',[-10 160],'ylim',[-20 50])
			set(gca,'fontsize',11)
		subplot(223)
			hold on; grid on
			plot(lziv*Ll(:,4),PLexp(size(Pd,1)+1:end,2),...
				'd','markersize',3,'linew',0.75)
			xlabel 'Axial stretch'; ylabel 'Transducer axial force [mN]'
			set(gca,'xlim',[1.2 2.2],'ylim',[0 50])
			set(gca,'fontsize',11)
		subplot(224)
			hold on; grid on
			plot(lziv*Ll(:,4),PLexp(size(Pd,1)+1:end,1)/mmHg_to_kPa,...
				'd','markersize',3,'linew',0.75)
			xlabel 'Axial stretch'; ylabel 'Pressure [mmHg]'
			set(gca,'xlim',[1.2 2.0],'ylim',[0 160])
			set(gca,'fontsize',11)
	end
	
	if save_indiv
		save(['biax_',grp,'_',num2str(sp),'.mat'],...
			'dotf','htf','lziv','rii','SM','emcM','ecA',...
			'meanpar','P','Pd','Ll','PLexp','dlexp');
	end
	ct = ct+1;
end

% Compute averaged data
if isnan(str2double(num))
	dotf    = mean(allspec(N,2))/1e3;
	htf     = mean(allspec(N,3))/1e3;
	lziv    = mean(allspec(N,4));
	rii     = mean(rii_all,1);
	Pd = squeeze(mean(Pd_all,1));
	Ll = squeeze(mean(Ll_all,1));
	Pd(:,4) = Pd(:,4)/lziv;
	Ll(:,4) = Ll(:,4)/lziv;
	PLexp = [Pd(:,[1 3]);Ll(:,[1 3])];
	dlexp = [Pd(:,[2 4]);Ll(:,[2 4])];

	if plot_data
		figure(1);
		set(gcf,'units','normalized','outerposition',[.3 .25 .4 .73]);
		subplot(221)
			hold on; grid on
			plot(Pd(:,2),PLexp(1:size(Pd,1),1)/mmHg_to_kPa,...
				'ks','markersize',5,'linew',1)
			xlabel 'Outer Diameter [mm]'; ylabel 'Pressure [mmHg]'
			set(gca,'xlim',[0.2 3.4],'ylim',[0 150])
			set(gca,'fontsize',11)
		subplot(222)
			hold on; grid on
			plot(PLexp(1:size(Pd,1),1)/mmHg_to_kPa,PLexp(1:size(Pd,1),2),...
				'ks','markersize',5,'linew',1)
			xlabel 'Pressure [mmHg]'; ylabel 'Transducer axial force [mN]'
			set(gca,'xlim',[-10 160],'ylim',[0 40])
			set(gca,'fontsize',11)
		subplot(223)
			hold on; grid on
			plot(lziv*Ll(:,4),PLexp(size(Pd,1)+1:end,2),...
				'ks','markersize',5,'linew',1)
			xlabel 'Axial stretch'; ylabel 'Transducer axial force [mN]'
			set(gca,'xlim',[1.2 2.2],'ylim',[0 40])
			set(gca,'fontsize',11)
		subplot(224)
			hold on; grid on
			plot(lziv*Ll(:,4),PLexp(size(Pd,1)+1:end,1)/mmHg_to_kPa,...
				'ks','markersize',5,'linew',1)
			xlabel 'Axial stretch'; ylabel 'Pressure [mmHg]'
			set(gca,'xlim',[1.2 2.2],'ylim',[0 150])
			set(gca,'fontsize',11)
	end
end

if save_group
	meanpar = [];
	save(['biax_',grp,'_',sgrp,'.mat'],...
		'dotf','htf','lziv','rii','SM','emcM','ecA',...
		'meanpar','P','Pd','Ll','PLexp','dlexp');
end
