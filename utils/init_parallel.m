function init_parallel(CoreNum)

% CoreNum: the number of CPU cores
feature('numCores')
delete(gcp('nocreate')) % open existing parallel
parpool('local', CoreNum); % open parallel
disp(['Initializing ', num2str(CoreNum), ' cores'])

end