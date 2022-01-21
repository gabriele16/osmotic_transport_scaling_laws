#! /bin/bash


echo "Extract free energy of K+ on hBN"

export BNKDIR=$PWD
export GRAKDIR=$PWD/../graphene_K_CV
export PATH=$PATH:"$PWD/../"

met_file="meta_data_file.dat"

echo "# metadata_file" > ${met_file}



for q in {5..30..1}
do

        cd ${BNKDIR}/us_$q

        loc_win_min=`grep -A 1 "COLVAR_RESTART" aimd_water_membrane.inp | tail -n 1 | awk '{print 0.0529177*$NF}'`
        spring=`grep "K \[" aimd_water_membrane.inp | awk '{print 2*$NF}'`
        colv_data_path=`readlink -f *COLVAR.metadynLog`
        
	echo ${colv_data_path} ${loc_win_min} ${spring} >> ${BNKDIR}/${met_file}

done

for q in {20..30..1}
do
	cd ${GRAKDIR}/us_${q} 

	loc_win_min=`grep -A1 "COLVAR_RESTART" GRAS-1.restart  | head -n 2 | tail -n 1 | awk '{print $0*0.0529177 }' `
	spring=`grep "K \[" aimd_water_membrane.inp | awk '{print 2*$NF}'`
	colv_data_path=`readlink -f *COLVAR.metadynLog`

	echo ${colv_data_path} ${loc_win_min} ${spring} >> ${BNKDIR}/${met_file}

done

cd ${BNKDIR}

umbrellaint.py ${met_file}  26 10000 40000 -T 300 -o FreeEnergyBNK.dat -c 3 --split no

