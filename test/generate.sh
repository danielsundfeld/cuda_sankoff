for i in ../seqs/RFAM*/*fasta; do
    NAME="$(basename $i)"
    FILE="${NAME}_gpu.sh"
    cp base.sh $FILE
    OUT=$(echo ${i:1} | sed 's/\//\\\//g')
    sed -i "s/__SEQ__/$OUT/" $FILE
    sed -i "s/__CMD__/.\/bin\/cuda_sankoff/" $FILE
    sed -i "s/__OPT__//" $FILE
    sed -i "s/__OUT__/gpu/" $FILE

    FILE="${NAME}_1_cpu.sh"
    cp base_cpu.sh $FILE
    OUT=$(echo ${i:1} | sed 's/\//\\\//g')
    sed -i "s/__SEQ__/$OUT/" $FILE
    sed -i "s/__CMD__/.\/bin\/sankoff/" $FILE
    sed -i "s/__OPT__/-t 1/" $FILE
    sed -i "s/__OUT__/cpu1/" $FILE
    
    FILE="${NAME}_16_cpu.sh"
    cp base_cpu.sh $FILE
    OUT=$(echo ${i:1} | sed 's/\//\\\//g')
    sed -i "s/__SEQ__/$OUT/" $FILE
    sed -i "s/__CMD__/.\/bin\/sankoff/" $FILE
    sed -i "s/__OPT__/-t 16/" $FILE
    sed -i "s/__OUT__/cpu16/" $FILE
done
