% Spike admin summary

spike_log = spike_log(spike_log.useNeuron == 1,:);
spike_log_troy = spike_log(strcmp(spike_log.monkey,'troy'),:);
spike_log_walt = spike_log(strcmp(spike_log.monkey,'walt'),:);


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
        labels{labels_i} = [labels{labels_i} ' | n = ' int2str(counts(monkey_i,labels_i))]
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
