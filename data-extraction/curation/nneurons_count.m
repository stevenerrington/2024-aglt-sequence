% Spike admin summary
spike_log = spike_log(spike_log.useNeuron == 1,:);
spike_log_troy = spike_log(strcmp(spike_log.monkey,'troy'),:);
spike_log_walt = spike_log(strcmp(spike_log.monkey,'walt'),:);

n_troy_sessions = sum(strcmp(ephysLog.monkey,'troy'));
n_walt_sessions = sum(strcmp(ephysLog.monkey,'walt'));

n_troy_auditory = sum(strcmp(spike_log_troy.area,'R') | strcmp(spike_log_troy.area,'dSTS'));
n_troy_frontal = sum(strcmp(spike_log_troy.area,'44') | strcmp(spike_log_troy.area,'45') | strcmp(spike_log_troy.area,'FOP'));

n_walt_auditory = sum(strcmp(spike_log_walt.area,'R') | strcmp(spike_log_walt.area,'dSTS'));
n_walt_frontal = sum(strcmp(spike_log_walt.area,'44') | strcmp(spike_log_walt.area,'45') | strcmp(spike_log_walt.area,'FOP'));


figuren;
subplot(1,2,1)
donut([n_troy_auditory, n_troy_frontal])
legend({'Aud','Frontal'})

subplot(1,2,2)
donut([n_walt_auditory, n_walt_frontal])
legend({'Aud','Frontal'})


for session_i = 1:size(ephysLog,1)
    session_name = ephysLog.session{session_i};
    count_session = [];
    count_session = [sum(strcmp(spike_log.area(strcmp(spike_log.session,session_name)),'R') | strcmp(spike_log.area(strcmp(spike_log.session,session_name)),'dSTS')),...
        sum(strcmp(spike_log.area(strcmp(spike_log.session,session_name)),'44') | strcmp(spike_log.area(strcmp(spike_log.session,session_name)),'45'))];

    n_pairs(session_i,1) = min(count_session);
end

sum(n_pairs(strcmp(ephysLog.monkey,'troy')))
sum(n_pairs(strcmp(ephysLog.monkey,'walt')))

area_labels = unique(spike_log.area)';
area_labels = area_labels([1,2,4,3,5,6,7,8]);

for area_i = 1:length(area_labels)
    counts(1,area_i) = sum(strcmp(spike_log_troy.area,area_labels{area_i}));
    counts(2,area_i) = sum(strcmp(spike_log_walt.area,area_labels{area_i}));
end



colors_in = [4 117 111; 26 83 105; 41 37 52; 242 96 12; 131 39 45; 115 34 43; 222 99 86; 160 88 161]./255


figuren('Renderer', 'painters', 'Position', [100 100 1200 500]); hold on;

for monkey_i = 1:2

    subplot(1,2,monkey_i)
    p_units_area = counts(monkey_i,:)/sum(counts(monkey_i,:));

    data = [p_units_area(1:3)',p_units_area(4:6)'];
    labels = [area_labels(1:3)',area_labels(4:6)'];

    for labels_i = 1:6
        labels{labels_i} = [labels{labels_i} ' | n = ' int2str(counts(monkey_i,labels_i))];
    end

    colors_seg = [{colors_in(1:3,:)}, {colors_in(4:6,:)}];

    % Lay out the column totals
    level1 = sum(data);

    cla reset
    r = treemap(level1);

    % Lay out each column within that column's rectangle from the overall
    % layout
    for j = 1:n
        colors = colors_seg{j};
        rNew = treemap(data(:,j),r(3,j),r(4,j));
        rNew(1,:) = rNew(1,:) + r(1,j);
        rNew(2,:) = rNew(2,:) + r(2,j);
        plotRectangles(rNew,labels(:,j),colors)
    end

    switch monkey_i
        case 1
            title('Monkey T')
        case 2
            title('Monkey W')
    end
end
outline(r)
axis([-0.01 1.01 -0.01 1.01])
