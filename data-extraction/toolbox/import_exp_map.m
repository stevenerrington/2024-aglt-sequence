function [ephysLog, stimulusLog, spike_log] = import_exp_map()
ephysLog = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'agl_t'));

stimulusLog = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'agl_t_stimuli'));

spike_log = webread(sprintf('https://docs.google.com/spreadsheets/d/%s/gviz/tq?tqx=out:csv&sheet=%s',...
    '1_kpK6t0yXWO5wVneRrX4kspHJXAnouSg', 'agl_t_spikes'));

end




