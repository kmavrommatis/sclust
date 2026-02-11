# Sclust
Code for SCLUST and modifications

- user provides directory with partition data (--dir)
- bamprocess is multithreading (--threads)


### Protocol for use:

[Nature Protocols](https://www.nature.com/articles/nprot.2018.033)


Step 1. Extract the read ratio and SNP information of the chromosome (<chr>) from the .bam-files by following this command:

```
Sclust bamprocess -t <sample>_T.bam -n <sample>_N.bam -o <sample> -part 1 -build hg38 -r <chr>
```
After completion of all chromosomes, generate temporary data files using this command:
```
Sclust bamprocess -build hg19 -i <sample>_chr1_bamprocess_data.txt,...,<sample>_chrY_bamprocess_data.txt
```

### Update use

Sclust bamprocess:

```
Sclust bamprocess -dir /lab/ni/mavrommk/tests/sclust/Sclust_1.1/ -part 1 -build hg38 -t <sample>_T.bam -n <sample>_N.bam --threads 18 -r all -o <sample>
```

The second call to bamprocess is no longer required.


