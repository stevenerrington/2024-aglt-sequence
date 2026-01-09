function [Xb, Yb] = balance_classes(X, Y)
% Balance classes by subsampling to the minimum class count

classes = categories(Y);
nClasses = numel(classes);

counts = zeros(nClasses,1);
idx_by_class = cell(nClasses,1);

for c = 1:nClasses
    idx_by_class{c} = find(Y == classes{c});
    counts(c) = numel(idx_by_class{c});
end

minCount = min(counts);

sel_idx = [];

for c = 1:nClasses
    idx = idx_by_class{c};
    sel_idx = [sel_idx; idx(randperm(numel(idx), minCount))];
end

sel_idx = sel_idx(randperm(numel(sel_idx)));

Xb = X(sel_idx,:);
Yb = Y(sel_idx);

end
