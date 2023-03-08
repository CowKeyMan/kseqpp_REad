# kseqpp_REad
A rewrite of kseqpp to only read sequences, and only read until a certain amount of characters, making it more scalable. Original kseqpp_REad can be found here: https://github.com/cartoonist/kseqpp

## Introduction

This implementation of kseqpp works by setting a maximum for the number of maximum characters that can be read from the filestream at one time. Thus, each Seq will contain as many characters as specified or less (in the case of the last read sequences).

Since it would not be possible to tell where one sequence starts and when does the next one end, the vector of positive numbers 'chars_before_newline' contain the indexes of the last character for each sequence within the 'seq' string.

So for example:
* given the seq 'ABCDEFG' and chars_before_newline [3, 7], we would have 2 strings which are 'ABC' and 'DEFG'
* given the seq 'ABCDEFG' and chars_before_newline [1, 1, 4], we would have 3 strings which are 'A', 'BCD' and 'EFG', where the final string (EFG) will continue when we continue reading from the file.

I have also stripped everything else unrelated to reading the sequences. Error checking has been left as is.

## Installation

This library is a single header, as such you can easily just download the header and use it. However it can also be integrated with CMake, and in order to do this, you can include the following in your project:

```Cmake
find_package(ZLIB)

include(FetchContent)
FetchContent_Declare(
  reklibpp
  QUIET
  GIT_REPOSITORY       "https://github.com/CowKeyMan/kseqpp_REad"
  GIT_TAG              "v1.4.0"
  GIT_SHALLOW          TRUE
)
FetchContent_MakeAvailable(reklibpp)

add_executable(main main.cpp)
target_link_libraries(main kseqpp_read)
```

## Usage
Too use this in your project, you can do the following:

```c++
Seq record(bufsize);
SeqStreamIn iss(filename.c_str());
while (iss >> record) {
  // do stuff with the record described in the introduction
  record.clear(); // remember to clear the record at the end if you wish to reuse it
}
```

For some example usage, checkout `src/test.cpp` which contains some basic unit tests to make sure the program works on well formed fasta and fastq files.

## Benchmarks

To run the benchmarks, it is recommended to download this repository and run `./scripts/build.sh`. This will build the project and download the benchmark objects. To run the benchmakrs, you may then run `./build/bin/benchmark`

Unfortunately, as a result to make it more scalable, more checks need to be made and thus performance takes a small hit when compared to the original kseqpp. This was one benchmark which was made using the AMD Rome 7H12 CPU and a 3.8TB NVME drive.

```
reklibpp fasta iteration 1: 3211ms
reklibpp fasta iteration 2: 3144ms
reklibpp fasta iteration 3: 3140ms
reklibpp fasta iteration 4: 3140ms
reklibpp fasta iteration 5: 3140ms
Average: 3155ms

reklibpp fastq iteration 1: 3877ms
reklibpp fastq iteration 2: 3745ms
reklibpp fastq iteration 3: 3729ms
reklibpp fastq iteration 4: 3735ms
reklibpp fastq iteration 5: 3726ms
Average: 3762.4ms

klibpp fasta iteration 1: 2380ms
klibpp fasta iteration 2: 2380ms
klibpp fasta iteration 3: 2380ms
klibpp fasta iteration 4: 2379ms
klibpp fasta iteration 5: 2380ms
Average: 2379.8ms

klibpp fastq iteration 1: 3304ms
klibpp fastq iteration 2: 3310ms
klibpp fastq iteration 3: 3321ms
klibpp fastq iteration 4: 3340ms
klibpp fastq iteration 5: 3440ms
Average: 3343ms

reklibpp fasta.gz iteration 1: 14239ms
reklibpp fasta.gz iteration 2: 14236ms
reklibpp fasta.gz iteration 3: 14195ms
reklibpp fasta.gz iteration 4: 14218ms
reklibpp fasta.gz iteration 5: 14208ms
Average: 14219.2ms

reklibpp fastq.gz iteration 1: 15393ms
reklibpp fastq.gz iteration 2: 15349ms
reklibpp fastq.gz iteration 3: 15343ms
reklibpp fastq.gz iteration 4: 15357ms
reklibpp fastq.gz iteration 5: 15364ms
Average: 15361.2ms

klibpp fasta.gz iteration 1: 13557ms
klibpp fasta.gz iteration 2: 13572ms
klibpp fasta.gz iteration 3: 13598ms
klibpp fasta.gz iteration 4: 13608ms
klibpp fasta.gz iteration 5: 13618ms
Average: 13590.6ms

klibpp fastq.gz iteration 1: 15085ms
klibpp fastq.gz iteration 2: 15101ms
klibpp fastq.gz iteration 3: 15081ms
klibpp fastq.gz iteration 4: 15078ms
klibpp fastq.gz iteration 5: 15080ms
Average: 15085ms
```

The advantage of this implementation is that the strings read will be sequential in memory, which means that it can easily be transferred to the gpu, or parallelised in some form or another. Processing in general becomes easier.

## Test

To run the tests, use the same script to build as the benchmarks (`./scripts/build.sh`), and then run `./build/bin/test`. The project uses googletest.
