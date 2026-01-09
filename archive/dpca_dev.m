id_labels = {'A','C','D','F','G'};
pos_labels = {'position_1','position_2','position_3','position_4','position_5'};

firingRatesAverage = zeros(size(sdf_soundAlign_data,2),5,5,length(-200:800));

for n = 1:size(sdf_soundAlign_data,2)

    neuron_data = sdf_soundAlign_data{n};
    neuron_sdf = cell2mat(neuron_data(:,1));

    valid_idx = strcmp(neuron_data(:,5), 'nonviol') & cell2mat(neuron_data(:,9)) == 1;

    for idx = 1:5
        id_sdf(idx,:) = smooth(nanmean(neuron_sdf(strcmp(neuron_data(:,3), id_labels{idx}) & valid_idx,:), 1), 50)';
        pos_sdf(idx,:) = smooth(nanmean(neuron_sdf(strcmp(neuron_data(:,4), pos_labels{idx}) & valid_idx,:), 1), 50)';
    end

    % Assign into firingRatesAverage
    % id_sdf → 2nd dimension (id)
    % pos_sdf → 3rd dimension (pos)
    % time → 4th dimension
    for i = 1:5
        for p = 1:5
            clear idxTrials
            idxTrials = valid_idx & ...
                strcmp(neuron_data(:,3), id_labels{i}) & ...
                strcmp(neuron_data(:,4), pos_labels{p});

            trialNum(n,i,p) = sum(idxTrials);

            if sum(idxTrials) > 0
                firingRatesAverage(n,i,p,:) = smooth(nanmean(neuron_sdf(idxTrials,:),1), 50);
            else
                firingRatesAverage(n,i,p,:) = zeros(1, size(neuron_sdf,2)); % handle empty condition
            end
        end
    end
end

%%

combinedParams = {{1, 2, [1 2]}};
margNames = {'identity', 'position', 'position-identity'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', -200:800,...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);



