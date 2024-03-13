
clear; close all; clc;

group = 'MWT';
subgroup = '8+BAPN';
path = [group,'\'];
specs = dir([path,'*']); mats = {specs.name}'; mats = mats(3:end);

% allspec_data  = zeros(length(mats),12);
% allspec_param = zeros(length(mats),8);
allspec = zeros(length(mats),12);
all_unload_config = [];
all_load_config = [];

for ii = 1:length(mats)
	load([path,mats{ii}]);
	
	% Specimen ID
	spec_ID = mats{ii}(1:8);
	
	% Unloaded configuration quantities
	% 1       2        3
	% L(mm)   OD(um)   H(um)
	unload_config = properties.unloaded(end,:);
	all_unload_config = [all_unload_config; unload_config];
	% unload_config = unload_config(2:3);
	
	% Loaded configuration quantities
	% 1         2        3       4      5    6        7           8
	% P(mmHg)   od(um)   h(um)   lziv   lt   W(kPa)   sigz(kPa)   sigt(kPa)
	load_config = properties.loaded.run3;
	all_load_config = [all_load_config; load_config];
	load_od(ii) = load_config(2);
	load_config = load_config(4);
	% load_config = load_config(1:8);
	
	% Linearized stiffness
	% 1            2            3            4
	% ctttt(MPa)   czzzz(MPa)   cttzz(MPa)   ctztz(MPa)
	linear_stiff = stiffness_SoL.run3;
	load_circstiff(ii) = stiffness_SoL.run3(1);
	linear_stiff = linear_stiff(1:2);
	
	% Bulk 4-fiber family parameters
	% 1        2          3     4          5     6          7     8          9
	% c(kPa)   c1z(kPa)   c2z   c1t(kPa)   c2t   c1d(kPa)   c2d   alp(deg)   RMSE
	bulk_params = parameters(end,:);
	bulk_params = bulk_params(1:8);
	
	% allspec_data(ii,:)   = [unload_config load_config linear_stiff];
	% allspec_param(ii,:)  = bulk_params;
	allspec(ii,:) = [unload_config load_config bulk_params];
end

fprintf('Group: %s %s\n',group,subgroup);

fprintf('Unoaded:\t\t\t\t\tLoaded:\n');
fprintf('\tOD\t\tH\t\t\t\tod\t\th\t\tlziv\tsigt\tsigz\tW\n');
for ii = 1:size(all_unload_config,1)
	fprintf('\t%.1f\t%.1f\t\t\t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\n',...
		all_unload_config(ii,2),all_unload_config(ii,3),...
		all_load_config(ii,2),all_load_config(ii,3),...
		all_load_config(ii,4),all_load_config(ii,8),...
		all_load_config(ii,7),all_load_config(ii,6));
end
fprintf('----------------------\t----------------------------------------------------\n');
fprintf('\t%.1f\t%.1f\t\t\t%.1f\t%.1f\t%.2f\t%.1f\t%.1f\t%.1f\n',...
	mean(all_unload_config(:,2)),mean(all_unload_config(:,3)),...
	mean(all_load_config(:,2)),mean(all_load_config(:,3)),...
	mean(all_load_config(:,4)),mean(all_load_config(:,8)),...
	mean(all_load_config(:,7)),mean(all_load_config(:,6)));

% save([group,'-AllSpec.mat'],'allspec');
% scriptSpecMats
