#test file for rand decode
#toric code, size 5, 9, 25, 35,..,7,11,13


#size=9
#run directly
error_file=input
#run after BP decoding
#error_file=output_nonconverge

error_folder=data/toric/bp_decoding
#error_folder=data/toric/temp

stabilizer_folder=data/toric/stabilizer

#for size in {5..6..2}
for size in {5,7,9,11,13}
do
    echo start rand decoding for size ${size} x ${size} >>${error_folder}/decode.log
    date  >>${error_folder}/decode.log

    #need 0.00 for 100..900 and 0.0 for 1000..5000
    for i in {100..900..100}
    do
	./rand_decode3.out ${stabilizer_folder}/toric_S_z_size_${size}.mm ${stabilizer_folder}/toric_S_x_size_${size}.mm ${error_folder}/toric_S_size_${size}.mm_rate0.00${i}_${error_file} $i &
    done
    #wait
    for i in {1000..2500..100}
    do
	./rand_decode3.out ${stabilizer_folder}/toric_S_z_size_${size}.mm ${stabilizer_folder}/toric_S_x_size_${size}.mm ${error_folder}/toric_S_size_${size}.mm_rate0.0${i}_${error_file} $i &
    done
    wait
    echo half time point >>${error_folder}/decode.log
    date >>${error_folder}/decode.log

    for i in {2600..5000..100}
    do
	./rand_decode3.out ${stabilizer_folder}/toric_S_z_size_${size}.mm ${stabilizer_folder}/toric_S_x_size_${size}.mm ${error_folder}/toric_S_size_${size}.mm_rate0.0${i}_${error_file} $i &
    done

    wait   
    echo finish rand decoding size ${size} x ${size} >>${error_folder}/decode.log
    date >>${error_folder}/decode.log

done


