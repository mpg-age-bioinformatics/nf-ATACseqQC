# nf-flexbar

Create the test directory:
```
mkdir -p ~/nf_atacseq_test_paired/raw_data
```

Download the demo data:
```
cd ~/nf_atacseq_test_paired/
curl -J -O https://datashare.mpcdf.mpg.de/s/ftg8K5610cIHGpR/download
unzip bowtie2_output.zip 
```

Change the parameters in params_paired.json accordingly, e.g. change "project_folder" : "/nexus/posix0/MAGE-flaski/service/hpc/home/wangy/nf_atacseq_test_paired/" to "project_folder" : Users/YOURNAME/nf_atacseq_test_paired/"


Run the workflow:

```
PROFILE=studio
nextflow run nf-ATACseqQC -params-file ~/nf-ATACseqQC/params.slurm.json -entry images -profile ${PROFILE}  && \
nextflow run nf-ATACseqQC -params-file ~/nf-ATACseqQC/params.slurm.json -profile ${PROFILE}
```

```

## Contributing

Make a commit, check the last tag, add a new one, push it and make a release:
```
git add -A . && git commit -m "<message>" && git push
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push origin --tags
gh release create <tag> 
```