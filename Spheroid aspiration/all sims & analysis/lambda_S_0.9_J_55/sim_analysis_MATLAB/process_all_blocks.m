
clearvars

block_scripts = {'force_06', 'force_08', 'force_1'};

for i = 1:length(block_scripts)
    script_name = block_scripts{i};
    fprintf('Running %s.m...\n', script_name);
    clear filePaths DM_string  
    run(script_name);       
end

dm_tags = {'0.6', '0.8', '1.0'};

k_stats{1} = NaN(1, numel(dm_tags));
k_stats{2} = NaN(1, numel(dm_tags));
k_stats{3} = NaN(1, numel(dm_tags));
tau_stats{1} = NaN(1, numel(dm_tags));   
tau_stats{2} = NaN(1, numel(dm_tags));   
tau_stats{3} = NaN(1, numel(dm_tags));   
xi_f_stats{1} = NaN(1, numel(dm_tags));   
xi_f_stats{2} = NaN(1, numel(dm_tags));   
xi_f_stats{3} = NaN(1, numel(dm_tags)); 
xi_t_stats{1} = NaN(1, numel(dm_tags));   
xi_t_stats{2} = NaN(1, numel(dm_tags));   
xi_t_stats{3} = NaN(1, numel(dm_tags));
L_p_stats{1} = NaN(1, numel(dm_tags));   
L_p_stats{2} = NaN(1, numel(dm_tags));   
L_p_stats{3} = NaN(1, numel(dm_tags)); 


for i = 1:numel(dm_tags)
    id = str2double(dm_tags{i}); 
    var_lk = sprintf('lower_k_%s', dm_tags{i});
    var_hk = sprintf('higher_k_%s', dm_tags{i});
    var_mk = sprintf('mean_k_%s', dm_tags{i});
    var_lt = sprintf('lower_tau_%s', dm_tags{i});
    var_ht = sprintf('higher_tau_%s', dm_tags{i});
    var_mt = sprintf('mean_tau_%s', dm_tags{i});
    var_lxi_f = sprintf('lower_xi_f_%s', dm_tags{i});
    var_hxi_f = sprintf('higher_xi_f_%s', dm_tags{i});
    var_mxi_f = sprintf('mean_xi_f_%s', dm_tags{i});
    var_lxi_t = sprintf('lower_xi_t_%s', dm_tags{i});
    var_hxi_t = sprintf('higher_xi_t_%s', dm_tags{i});
    var_mxi_t = sprintf('mean_xi_t_%s', dm_tags{i});
    var_lL_p = sprintf('lower_L_p_%s', dm_tags{i});
    var_hL_p = sprintf('higher_L_p_%s', dm_tags{i});
    var_mL_p = sprintf('mean_L_p_%s', dm_tags{i});

    if exist(var_lk,'var');  k_stats{1}(i) = eval(var_lk); end
    if exist(var_mk,'var');  k_stats{2}(i) = eval(var_mk); end
    if exist(var_hk,'var');  k_stats{3}(i) = eval(var_hk); end

    if exist(var_lt,'var');  tau_stats{1}(i) = eval(var_lt); end
    if exist(var_lt,'var');  tau_stats{2}(i) = eval(var_mt); end
    if exist(var_lt,'var');  tau_stats{3}(i) = eval(var_ht); end

    if exist(var_lxi_f,'var');  xi_f_stats{1}(i) = eval(var_lxi_f); end
    if exist(var_lxi_f,'var');  xi_f_stats{2}(i) = eval(var_mxi_f); end
    if exist(var_lxi_f,'var');  xi_f_stats{3}(i) = eval(var_hxi_f); end

    if exist(var_lxi_t,'var');  xi_t_stats{1}(i) = eval(var_lxi_t); end
    if exist(var_lxi_t,'var');  xi_t_stats{2}(i) = eval(var_mxi_t); end
    if exist(var_lxi_t,'var');  xi_t_stats{3}(i) = eval(var_hxi_t); end

     if exist(var_lL_p,'var');  L_p_stats{1}(i) = eval(var_lL_p); end
    if exist(var_lL_p,'var');  L_p_stats{2}(i) = eval(var_mL_p); end
    if exist(var_lL_p,'var');  L_p_stats{3}(i) = eval(var_hL_p); end

end

fprintf('All blocks processed and saved.\n');