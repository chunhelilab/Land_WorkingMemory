function [] = parsave(dir,r_end_all, S_end_all, init_all)
%save variables in dir
% so I can save in parfor loop
save(dir,'r_end_all', 'S_end_all', 'init_all')
end