function setup_cores(request)
% SETUP_CORE setup the parent function to utilize a maximum number of cores
%
% setup_core(request) attempt to set max to request cores. Use 'auto' to
% set to system maximum

cores = feature('numcores');
cluster = parcluster;

% maxNumCompThread will be depreciated soon, but for now will use.
warning('off', 'MATLAB:maxNumCompThreads:Deprecated');

% set the max number of threads a single instance of MATLAB can use
if ischar(request) && strcmp(request, 'auto')
    request = cores;
    maxNumCompThreads('automatic');
else
    maxNumCompThreads(request);
end

% If for some reason a negative number is request set to 1
if request <= 0
   request = 1;
end
    
% If more was requested than there are cores reduce request
if request > cores
    disp('requested to many cores reducing to system max');
    request = cores;
    maxNumCompThreads(request);
end

% Setup parallel pool to reflect requested cores
if request == 1
    delete(gcp('nocreate'));
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

