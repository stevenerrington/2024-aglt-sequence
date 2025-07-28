

clear glm_sig_count_table n_just_identity  n_just_position n_both n_neither
for area_i = 1:2


    neurons_in = [];
    switch area_i
        case 1
            neurons_in = neuron_class.auditory.all;
        case 2
            neurons_in = neuron_class.frontal.all;
    end

    % Logical masks
    clear just_identity just_position both neither
    just_identity = glm_encoding_flag(neurons_in,1) == 1 & glm_encoding_flag(neurons_in,2) == 0;
    just_position = glm_encoding_flag(neurons_in,1) == 0 & glm_encoding_flag(neurons_in,2) == 1;
    both = glm_encoding_flag(neurons_in,1) == 1 & glm_encoding_flag(neurons_in,2) == 1;
    neither = glm_encoding_flag(neurons_in,1) == 0 & glm_encoding_flag(neurons_in,2) == 0; % optional

    % Count neurons in each category
    n_just_identity(1, area_i) = sum(just_identity);
    n_just_position(1, area_i) = sum(just_position);
    n_both(1, area_i) = sum(both);
    n_neither(1, area_i) = sum(neither); % optional



end

glm_sig_count_table = [];
glm_sig_count_table = [n_just_identity, n_just_position, n_both];
glm_sig_count_table = glm_sig_count_table(:, [1 3 5 2 4 6]);


glm_sig_count_table  = array2table(glm_sig_count_table, 'VariableNames',{'n_aud_identity','n_aud_order',...
    'n_aud_both','n_frontal_identity','n_frontal_order','n_frontal_both'});

glm_counts = sum(glm_sig_count_table);

figuren('Renderer', 'painters', 'Position', [100 207 1200 400]);;
subplot(1,2,1)
donut([glm_counts.n_aud_identity, glm_counts.n_aud_order, glm_counts.n_aud_both]);
axis square; xlim([-2 2]); ylim([-2 2]); title('Auditory')
subplot(1,2,2)
donut([glm_counts.n_frontal_identity, glm_counts.n_frontal_order, glm_counts.n_frontal_both]);
axis square; xlim([-2 2]); ylim([-2 2]); title('Frontal')



%%
x = table([glm_counts.n_aud_identity; glm_counts.n_aud_order; glm_counts.n_aud_both],...
    [glm_counts.n_frontal_identity; glm_counts.n_frontal_order; glm_counts.n_frontal_both],...
    'VariableNames',{'auditory','frontal'},'RowNames',{'identity','order','both'});

% Replace with your actual data
data = [
    glm_counts.n_aud_identity, glm_counts.n_frontal_identity;  % Identity only
    glm_counts.n_aud_order, glm_counts.n_frontal_order;  % Position only
    glm_counts.n_aud_both, glm_counts.n_frontal_both   % Identity + Position
];

% Chi-squared test of independence
[chi2stat, p, stats] = chi2cont(data);

disp(['Chi-squared p-value: ', num2str(p)])
disp(stats)

data(:,1) = data(:,1) ./ sum(data(:,1))
data(:,2) = data(:,2) ./ sum(data(:,2))

%%
