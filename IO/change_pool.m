function pool = change_pool(request)
% CHANGE_POOL change the number of matlab workers currently running
%
% pool = CHANGE_POOL(request) change pool to request workers

pool = gcp('nocreate');
cluster = parcluster;
cores = feature('numcores');

if ischar(request) && strcmp(request, 'auto')
    request = cores;
end

% Delete pool if request is zero otherwise change to requested
% If for some reason a negative number is request set to 1
if request <= 0
   request = 1;
end
    
% If more was requested than there are cores reduce request
if request > cores
    disp('requested to many cores reducing to system max');
    request = cores;
end

% Setup parallel pool to reflect requested cores
if request == 1
    delete(gcp);
else
    pool = gcp('nocreate');
    if isempty(pool);
        parpool(cluster, request);
    end

    if ~isempty(pool) && pool.NumWorkers ~= request
        delete(gcp);
        parpool(cluster, request);
    end
end

end