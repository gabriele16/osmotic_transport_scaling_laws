#! /bin/bash

echo  "Extract free energy of I- on hBN"

export BNIDIR=$PWD
export GRAIDIR=$PWD/../graphene_I_CV
export PATH=$PATH:"$PWD/../"

met_file="meta_data_file.dat"
colvar_file='colvar_data.dat'

echo "# metadata_file" > ${met_file}


for q in {5..30..1}
do

        cd ${BNIDIR}/us_$q

        loc_win_min=`grep -A 1 "COLVAR_RESTART" aimd_water_membrane.inp | tail -n 1 | awk '{print 0.0529177*$NF}'`
        colv_data_path=`readlink -f *COLVAR*`
        spring=`grep "K \[" aimd_water_membrane.inp | awk '{print 2*$NF}'`

        colv_data_path=`readlink -f *COLVAR.metadynLog`
        echo ${colv_data_path} ${loc_win_min} ${spring} >> ${BNIDIR}/${met_file}

done


for q in {20..30..1}
do
	cd ${GRAIDIR}/us_${q} 

	coord_surf=`head -n 3 mod_gra_wat_i.xyz | awk '{print $4}' | tail -n 1`
	loc_win_min=`tail -n 1 mod_gra_wat_i.xyz | awk '{print ($4 - "'"${coord_surf}"'")*0.1}'`
	colv_data_path=`readlink -f *COLVAR*`
	spring=`grep "K \[" aimd_water_membrane.inp | awk '{print 2*$NF}'`

	cat *COLVAR.metadynLog | awk '{i+=1; print i, 0.0529177*$2}'  > ${colvar_file}

        colv_data_path=`readlink -f *COLVAR.metadynLog`

	echo ${colv_data_path} ${loc_win_min} ${spring} >> ${BNIDIR}/${met_file}

done


cd ${BNIDIR}

umbrellaint.py ${met_file}  26 10000 70000 -T 300 -o FreeEnergyBNI.dat -c 3 --split no

