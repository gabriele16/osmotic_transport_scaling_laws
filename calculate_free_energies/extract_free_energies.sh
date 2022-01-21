#! /bin/bash

for dd in graphene_I_CV graphene_K_CV hbn_I_CV hbn_K_CV
do	
	cd $dd 
	./run_ubint.sh
	cd ..
done
