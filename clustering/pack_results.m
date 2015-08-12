function results_scores = pack_results(results_scores, frob_km, frob_gmm, ...
    like_km, like_gmm, completed, km_steady, gmm_steady, model, sub_model, file)
% PACK_RESULTS store all variables into the results_scores struct, really
% just to clearn up the data
%
%   results_scores = PACK_RESULT(results_scores, completed, frob_km, ...
%       frob_gmm, like_km, like_gmm, fil, model, sub_model)


% Pack variables
results_scores.(file.name).(model).(sub_model).frob_km     = frob_km.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).frob_gmm    = frob_gmm.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).like_km     = like_km.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).like_gmm    = like_gmm.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).completed   = completed.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).steady   = km_steady.(model).(sub_model);


end