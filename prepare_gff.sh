
chr="000000F"
#--------------------------------------------------
# awk -v chr=$chr '$5==chr{print $5,$6,$7,$11}' ~/wlflab/Ming/genome_assembly/EonycterisSpelaea/falcon.1000.2.400.nscc/2-asm-falcon.14000.md300.xc400.nc2.b20.i80/eonSpe2.mask2/eonSpe2.fa.out > rep.data 
#--------------------------------------------------
for gff in ~/batlab/Ming/genome_assembly/EonycterisSpelaea/methylation/ipd_out/*$chr.gff* 
do
    if [[ ! -e `basename $gff`.m6c ]]; then
        sed '1,4d' $gff  > `basename $gff`.m6c
    fi
done
for gff in ~/batlab/Ming/genome_assembly/EonycterisSpelaea/methylation/AgIn_out/*$chr.gff*
do
    if [[ ! -e `basename $gff`.m5c ]]; then
        sed '1,5d' $gff  > `basename $gff`.m5c
    fi
done
if [[ ! -e kidney.000000F_class.wig ]]; then
    ln -s /Users/gmswenm/batlab/Ming/genome_assembly/EonycterisSpelaea/methylation/AgIn_out/kidney.000000F_class.wig
    ln -s /Users/gmswenm/batlab/Ming/genome_assembly/EonycterisSpelaea/methylation/AgIn_out/lung.000000F_class.wig
fi

