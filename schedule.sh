#to run a job after some time
date
sleep 1200
date
echo start convert data
./convert_data.sh >>a.txt
date
echo start rand decode
./run_rand3.sh >>a.txt
date
echo finish schedule >>a.txt
