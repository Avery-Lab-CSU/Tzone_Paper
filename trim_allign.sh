# trim
trim_galore -j 4 --paired \
../${names[${SLURM_ARRAY_TASK_ID}]}1.fq.gz \
../${names[${SLURM_ARRAY_TASK_ID}]}2.fq.gz

#Align count
STAR \
--genomeDir ../str_cf3.100 \
--outSAMtype None \
--quantMode GeneCounts \
--runThreadN 5 \
--readFilesCommand zcat \
--outFileNamePrefix ${sid[${SLURM_ARRAY_TASK_ID}]}_ \
--readFilesIn \ 
../trim/${sid[${SLURM_ARRAY_TASK_ID}]}_1_val_1.fq.gz \
../trim/${sid[${SLURM_ARRAY_TASK_ID}]}_2_val_2.fq.gz