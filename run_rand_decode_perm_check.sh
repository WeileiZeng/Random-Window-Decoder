#test expotential decay with # of permutation
#size =5, p={1000,2000,3000}, perm ={5..100..5}
#direct rand decode on toric code

#switch input file
error_file=input
#error_file=output_nonconverge

error_folder=data/toric/permutation/input_error
stabilizer_folder=data/toric/stabilizer


#for size in {5,7,9,11,13}
size=5
for p in {1000,2000,3000}
#for p in {3000..3100..500}
do
    echo start check for size ${size} x ${size} when p = $p>>${error_folder}/decode.log
    date  >>${error_folder}/decode.log

    output_folder=data/toric/permutation/output_error/p${p}
 #   mkdir ${output_folder}
    
    #name rule: need 0.00 for 100..900 and 0.0 for 1000..5000
    for j in {2..20..1}
    do
	#to make this work, one need rename the first file in perm1/ from *output to *output_perm1
	output_folder1=data/toric/permutation/output_error/p${p}/perm$(($j-1))
	output_folder2=data/toric/permutation/output_error/p${p}/perm${j}
	output_file1=${error_file}_output_perm$(($j-1))
	output_file2=${error_file}_output
	output_file=${error_file}_output_perm${j}
	./rand_decode_perm_check.out \
	    ${stabilizer_folder}/toric_S_z_size_${size}.mm \
	    ${error_folder}/toric_S_size_${size}.mm_rate0.0${p}_${error_file} \
	    ${output_folder1}/toric_S_size_${size}.mm_rate0.0${p}_${output_file1} \
	    ${output_folder2}/toric_S_size_${size}.mm_rate0.0${p}_${output_file2} \
	    ${output_folder2}/toric_S_size_${size}.mm_rate0.0${p}_${output_file} \
	    ${output_folder}/toric_S_size_${size}.mm_rate0.0${p}_${error_file} $p 

	echo p = ${p}, j = $j  >>${error_folder}/decode.log
    done
   
    
    echo finish rand decoding size ${size} x ${size} when p = $p >>${error_folder}/decode.log
    date >>${error_folder}/decode.log

done


