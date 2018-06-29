#This file print the out put of bp_decoding2.c, trying to figureout the local convegence in the case where it doesn't reach overall convergence

#try bp decoding for toric code with different size, 5,9,25,35
#p=0.01% to 0.20%
#100,000 division, 2000 means 2%
#./bp_decoding2.out data/toric/toric_S_size_5.mm data/toric/bp_converge3/toric_S_size_5.mm 500

./bp_decoding3.out data/toric/stabilizer/toric_S_x_size_5.mm data/toric/temp/toric_S_size_5.mm 2000

#wait
#for i in {26..50..1}
#do
    
#    ./bp_decoding.out data/toric/toric_S_size_35.mm data/toric/bp_converge/toric_S_size_35.mm $i &
#done
#wait

echo finish bp_decoding 10000 cycle for size 5 x 5 when p is \in 100,000 division
#>> data/toric/bp_converge3/bp_decoding.log &
#date
#>> data/toric/bp_converge2/bp_decoding.log &

