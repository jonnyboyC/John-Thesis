function results_scores = pack_results(results_scores, completed, frob_km, frob_gmm, ...
    prob_km, prob_gmm, model, sub_model, file)
% PACK_RESULTS store all variables into the results_scores struct, really
% just to clearn up the data
%
%   results_scores = PACK_RESULT(results_scores, completed, frob_km, ...
%       frob_gmm, prob_km, prob_gmm, fil, model, sub_model)


% Pack variables
results_scores.(file.name).(model).(sub_model).frob_km     = frob_km.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).frob_gmm    = frob_gmm.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).prob_km     = prob_km.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).prob_gmm    = prob_gmm.(model).(sub_model);
results_scores.(file.name).(model).(sub_model).completed   = completed.(model).(sub_model);

end