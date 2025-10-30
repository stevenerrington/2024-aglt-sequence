function filters = design_filters(freqs_in)

% construct a cell array containing all of the bandpass filters
filters = {};

for f = 1:numel(freqs_in(:,1))

    l_bound = freqs_in(f,1);

    u_bound = freqs_in(f,2);

    filters{f} = designfilt('bandpassfir', ... response type
        'FilterOrder',1000, ... filter order
        'CutoffFrequency1', l_bound,... lower range
        'CutoffFrequency2', u_bound,... upper range
        'SampleRate', 1000); ... sampling rate

end % of looping over filters to construct
