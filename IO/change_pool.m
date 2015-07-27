function pool = change_pool(request)
% CHANGE_POOL change the number of matlab workers currently running
%
% pool = CHANGE_POOL(request) change pool to request workers

pool = gcp('nocreate');
cluster = parcluster('local');

% Delete pool if reques is zero otherwise change to requested
if request == 0
    delete(gcp);
elseif pool.NumWorkers ~= request
    delete(gcp);
    parpool(cluster, request);
end
end