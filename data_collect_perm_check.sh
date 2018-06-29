
#for size in {5,7,9,11,13}
p=1000
echo p = $p
for perm in {1..20..1}
do
    error_folder=data/toric/permutation/output_error/p${p}/perm${perm}
    filename_output=toric_S_size_5.mm_rate0.0${p}_input_output
    filename_output_bad=toric_S_size_5.mm_rate0.0${p}_input_output_perm${perm}_output_bad
    filename_gnudat=data/toric/permutation/gnudat/rate_versus_perm_p${p}.gnudat
    ./data_collect_perm_check.out ${error_folder}/${filename_output} ${error_folder}/${filename_output_bad} ${filename_gnudat} ${perm}
    
done
