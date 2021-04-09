function transient = isolate_transient(data, signal, n_events)
    headings = data(1,:);
    transient = cell(n_events+1,numel(headings));
    transient(1,:) = headings; 
    n = 30;
    data_loc = strcmpi(headings,signal);
    sig = data{2,data_loc};
    
    if strcmpi(signal, 'id')
% %         edges = findchangepts(sig, 'maxnumchanges', n_events, 'Statistic','linear');
        edges = findchangepts(smooth(sig,30)', 'maxnumchanges', n_events+1, 'Statistic','linear');
    else
        edges = findchangepts(sig, 'maxnumchanges', n_events+1); 
    end
    
    coarse_idx = round(mean([1, edges; edges, length(sig)]));
    for event = 1:n_events
        cur_sig = sig(coarse_idx(event):coarse_idx(event+1));
        samp = round(length(cur_sig)/20);   % consider five percent of the signal
        samp1 = cur_sig(1:samp); 
        samp2 = cur_sig(end+1-samp:end);
%         d1 = std(abs(diff(samp1)));
%         d2 = std(abs(diff(samp2)));
        d1 = max([std(abs(diff(samp1))),1e-1]);
        d2 = max([std(abs(diff(samp2))),1e-1]);
        sig1 = cur_sig(1:edges(event)-coarse_idx(event));
        sig2 = cur_sig(edges(event)-coarse_idx(event)+1:end);
        %%% Jan 28, 2020: update < to <= to account for no noise
        s1 = nlfilter(abs(diff(sig1))<=2*d1,[1,n],@all);
        s2 = nlfilter(abs(diff(sig2))<=2*d2,[1,n],@all);
        s1(1) = 1; 
        s2(end) = 1;
        tran_start = coarse_idx(event) - 1 + find(s1,1,'last');
        tran_end = edges(event) - 1 + find(s2,1);
        tran_idx = tran_start:tran_end; 
        for k = 1:numel(headings)
            transient{event+1,k} = data{2,k}(tran_idx);
        end
    end
end
