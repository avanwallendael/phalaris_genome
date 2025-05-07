# Circos Generation for _Phalaris minor_ IWGC Genome

## Setup
Use singularity or Docker to pull the required containers
  
### iwgc_circos_tracks.sh
```
singularity pull bedtools_v2.27.1dfsg-4-deb_cv1.sif docker://biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1
singularity pull samtools_1.21.sif docker://staphb/samtools:1.21
singularity pull quartet_1.1.7.sif docker://changshengw6/quartet:1.1.7

```

### Circos
Ran locally (2021 Macbook M1 Pro, 10 core, 32GB memory, Sequoia 15.4.1 (24E263)) with Docker
```

docker pull staphb/circos
```

## Run
### iwgc_circos_tracks.sh
```
./pminor_circos_tracks.sh \
Avena_fatua_Chromosomes.fasta AveFa_v01.0.gff Avena_fatua_Chromosomes.fasta.mod.EDTA.TEanno.gff3 Avena_fatua_Chromosomes.fasta.mod.EDTA.intact.gff3 \
samtools_1.21.sif bedtools_v2.27.1dfsg-4-deb_cv1.sif quartet.v1.1.7.sif \
3000000
```

  
### Circos

```
docker run --rm -v "$PWD:/data" staphb/circos \                       
  circos -conf /data/pminor_circos/pminor_circos.conf -outputdir /data/pminor_circos/tmp -noparanoid
```
