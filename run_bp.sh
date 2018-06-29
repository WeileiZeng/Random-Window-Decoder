#try bp decoding for toric code with different size, 5,9,25,35
#p=0.01% to 0.20%
for i in {10..200..10}
do    
    ./bp_decoding.out data/toric/toric_S_size_9.mm data/toric/bp_converge2/toric_S_size_9.mm $i &
done
wait
#for i in {26..50..1}
#do
    
#    ./bp_decoding.out data/toric/toric_S_size_35.mm data/toric/bp_converge/toric_S_size_35.mm $i &
#done
#wait

echo finish bp_decoding 1000 cycle for size 9 x 9 when p is \in 100,000 division >> data/toric/bp_converge2/bp_decoding.log &
date >> data/toric/bp_converge2/bp_decoding.log &

