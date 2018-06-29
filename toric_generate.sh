#this one can generate toric code with x and z check matix seperately or togetther. the switch is inside the code_generate.c
#genrating some toric codes with different dimension, saved into data_folder
data_folder=data/toric/stabilizer

for size in {2,3,5,7,9,11,13,25,35}	    

do
    #need to recompile code_generate.out when switch between x and z    
#    ./code_generator.out ${size} ${data_folder}/toric_S_x_size_${size}.mm &
    ./code_generator.out ${size} ${data_folder}/toric_S_z_size_${size}.mm &
done
#size=25
#./code_generator.out ${size} ${data_folder}/toric_S_z_size_${size}.mm
echo done
date
