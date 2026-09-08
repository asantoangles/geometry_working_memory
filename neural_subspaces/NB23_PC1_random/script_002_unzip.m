%% unzip

cd /path_to_local/results/neural_subspaces/NB23_PC1_random

for rand_i = 1:1000

    system(['unzip rand' num2str(rand_i) '.zip']);

end