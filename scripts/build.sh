#!/bin/bash

mkdir -p build
cd build
cmake \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_CLANG_TIDY=ON \
  -DENABLE_HEADER_GUARDS_CHECK=ON \
  -DENABLE_CLANG_FORMAT_CHECK=ON \
  -DKSEQPP_READ_BUILD_TESTS=ON \
  -DKSEQPP_READ_BUILD_BENCHMARKS=ON \
  ..
cmake --build . -j8
cd ..


# download benchmark_objects
mkdir -p benchmark_objects
cd benchmark_objects

# Download Fasta
FASTA="FASTA.fna"
## Source: https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml
wget -nc "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz" -O "${FASTA}.gz"
if test -f "${FASTA}"; then
  echo "${FASTA} already unzipped, skipping "
else
  echo "Unzipping ${FASTA}"
  gzip -dkf "${FASTA}.gz" > "${FASTA}"
fi

# Download Fastq
FASTQ="FASTQ.fnq"
## Source: https://www.ebi.ac.uk/ena/browser/view/SRX11174563?show=reads
wget -nc "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/030/SRR14856430/SRR14856430_1.fastq.gz" -O "${FASTQ}.gz"
if test -f "${FASTQ}"; then
  echo "${FASTQ} already unzipped, skipping "
else
  echo "Unzipping ${FASTQ}"
  gzip -dkf "${FASTQ}.gz" > "${FASTQ}"
fi

cd ..

