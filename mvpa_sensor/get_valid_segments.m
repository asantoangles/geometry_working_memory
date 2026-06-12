function vals_out = get_valid_segments(data, p_thr, min_seg)

% data: (trial × time × stim)

vals_out = [];

[nTrial, nTime, nStim] = size(data);

for tr = 1:nTrial
    for st = 1:nStim

        ts = squeeze(data(tr,:,st));  % (time)

        % threshold mask
        mask = ts > p_thr;

        % find contiguous segments
        d = diff([0 mask 0]);
        start_idx = find(d == 1);
        end_idx   = find(d == -1) - 1;

        for i = 1:length(start_idx)

            seg_len = end_idx(i) - start_idx(i) + 1;

            if seg_len >= min_seg

                seg = ts(start_idx(i):end_idx(i));
                vals_out = [vals_out; seg(:)];

            end

        end

    end
end

end