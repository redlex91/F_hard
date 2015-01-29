#!/bin/bash

for(( count = 20; count<=240; count = count +1 ))
do
	echo "-" $count
	echo $count > temp
	./main.out < temp
done

rm temp
