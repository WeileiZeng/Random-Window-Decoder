#direct rand decoder with perm = 5, e_try=10,000,size = 5
#run it 20 times and save result in 20 folders
#test expotential decay with # of permutation

error_file=input
#error_file=output_nonconverge

error_folder=data/toric/permutation/input_error
#error_folder=data/toric/temp

stabilizer_folder=data/toric/stabilizer

#for size in {5..6..2}
#for size in {5,7,9,11,13}

size=5
for p in {1000,2000,3000}
do
    echo start rand decoding for size ${size} x ${size} when p = $p>>${error_folder}/decode.log
    date  >>${error_folder}/decode.log

    output_folder=data/toric/permutation/output_error/p${p}
    mkdir ${output_folder}
    
    #need 0.00 for 100..900 and 0.0 for 1000..5000
#    for i in {100..900..100}
    for j in {1..20..1}
    do
	output_folder=data/toric/permutation/output_error/p${p}/perm${j}
	mkdir ${output_folder}
	./rand_decode_perm.out ${stabilizer_folder}/toric_S_z_size_${size}.mm ${stabilizer_folder}/toric_S_x_size_${size}.mm ${error_folder}/toric_S_size_${size}.mm_rate0.0${p}_${error_file} ${output_folder}/toric_S_size_${size}.mm_rate0.0${p}_${error_file} $p &
    done
   
    wait   
    echo finish rand decoding size ${size} x ${size} when p = $p >>${error_folder}/decode.log
    date >>${error_folder}/decode.log

done


