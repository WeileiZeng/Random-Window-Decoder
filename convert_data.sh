error_folder=data/toric/bp_decoding

#./convert_data.out ${error_folder} 5
for size in {5,7,9,11,13}
do
    ./convert_data.out ${error_folder} ${size} &
done
