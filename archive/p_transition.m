
sound_list = {'A','C','D','F','G'};

for list_i = 1:5

    p_ele_pos(1, list_i) = 100*(sum(strcmp(stimulusLog.sound_1_code(1:8),sound_list{list_i}))./8);
    p_ele_pos(2, list_i) = 100*(sum(strcmp(stimulusLog.sound_2_code(1:8),sound_list{list_i}))./8);
    p_ele_pos(3, list_i) = 100*(sum(strcmp(stimulusLog.sound_3_code(1:8),sound_list{list_i}))./8);
    p_ele_pos(4, list_i) = 100*(sum(strcmp(stimulusLog.sound_4_code(1:8),sound_list{list_i}))./8);
    p_ele_pos(5, list_i) = 100*(sum(strcmp(stimulusLog.sound_5_code(1:8),sound_list{list_i}))./8);


end

figure('Renderer', 'painters', 'Position', [100 100 500 500]); 
heatmap(p_ele_pos)