function timestamp = get_rec_timestamps(signal, ops)

switch signal
    case 'spk'
        clear filepart_name NCSpath spk_ncs_in timestamp

        filepart_name = ['CSC1_' ops.session_n];
        NCSpath = [fullfile(ops.raw_dir,ops.filename,filepart_name) '.ncs'];

        spk_ncs_in = readncs([filepart_name '.ncs'],fullfile(ops.raw_dir,ops.filename));
        timestamp(1,1) = spk_ncs_in.TimeStamp(1);
        timestamp(1,2) = spk_ncs_in.TimeStamp(end);        

    case 'events'
        [evnt_raw] = ft_read_event_BA(fullfile(ops.raw_dir,ops.filename,['Events_' ops.session_n '.nev']));
        timestamp(1,1) = evnt_raw(1).timestamp;
        timestamp(1,2) = evnt_raw(end).timestamp;
        
end

timestamp = double(timestamp);

end

