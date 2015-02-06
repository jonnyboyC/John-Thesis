function vars = load_data(type, direct_Gal)
% Loads correct data.results for the type of request for MOD_Pod

data = load(direct_Gal, 'results');
vars.OG_nm  = data.results.num_modesG;
vars.q      = data.results.q;
if any(strcmp(type, 'og'))
    vars.c = data.results.c;
    vars.l = data.results.l;
    vars.t = data.results.t1;
    vars.modal_amp = data.results.modal_amp_og;
elseif any(strcmp(type, 'vis1'))
    vars.c = data.results.ci;
    vars.l = data.results.li;
    vars.t = data.results.t2;
    vars.modal_amp = data.results.modal_amp_vis1;
elseif any(strcmp(type, 'vis2'))
    vars.c = data.results.ci_c;
    vars.l = data.results.li_c;
    vars.t = data.results.t3;
    vars.modal_amp = data.results.modal_amp_vis2;
end

end