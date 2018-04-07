for i in `seq 0 40`;
do
	omega=$(python -c "print(1.93+$i/2000.0)")

	./magnetic -w $omega -c 100 -r 100 --SOR -n 0
	sleep 1
done
