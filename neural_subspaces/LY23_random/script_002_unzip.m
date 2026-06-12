%% unzip

cd /path_to_local/github/results/source_geometry_lm/LY23_random

for rand_i = 1:1000

    system(['unzip rand' num2str(rand_i) '.zip']);

end