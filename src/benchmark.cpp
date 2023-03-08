#include <chrono>
std::chrono::system_clock::time_point TIME_IT_START_TIME, TIME_IT_END_TIME;
unsigned long long TIME_IT_TOTAL = 0;
#define TIME_IT(code_block)                                              \
  TIME_IT_START_TIME = std::chrono::high_resolution_clock::now();        \
  code_block;                                                            \
  TIME_IT_END_TIME = std::chrono::high_resolution_clock::now();          \
  TIME_IT_TOTAL = std::chrono::duration_cast<std::chrono::milliseconds>( \
                    TIME_IT_END_TIME - TIME_IT_START_TIME                \
  )                                                                      \
                    .count();

#include <iostream>

#include "kseq++/seqio.hpp"
#include "kseqpp_read.hpp"
using std::cout;
using std::endl;

auto benchmark_reklibpp_fasta() -> void;
auto benchmark_reklibpp_fastq() -> void;
auto benchmark_klibpp_fasta() -> void;
auto benchmark_klibpp_fastq() -> void;
auto benchmark_reklibpp_fasta_gz() -> void;
auto benchmark_reklibpp_fastq_gz() -> void;
auto benchmark_klibpp_fasta_gz() -> void;
auto benchmark_klibpp_fastq_gz() -> void;

const int iterations = 5;

int main() {
  benchmark_reklibpp_fasta();
  cout << endl;
  benchmark_reklibpp_fastq();
  cout << endl;
  benchmark_klibpp_fasta();
  cout << endl;
  benchmark_klibpp_fastq();
  cout << endl;
  benchmark_reklibpp_fasta_gz();
  cout << endl;
  benchmark_reklibpp_fastq_gz();
  cout << endl;
  benchmark_klibpp_fasta_gz();
  cout << endl;
  benchmark_klibpp_fastq_gz();
}

auto benchmark_reklibpp_fasta() -> void {
  double total = 0;
  for (int i = 0; i < iterations + 1; ++i) {
    auto record = reklibpp::Seq(16 * 1024 * 1024);
    auto ssi = reklibpp::SeqStreamIn("benchmark_objects/FASTA.fna");
    TIME_IT({
      while (ssi >> record) { record.clear(); }
    })
    if (i > 0) {
      cout << "reklibpp fasta iteration " << i << ": " << TIME_IT_TOTAL << "ms"
           << endl;
      total += double(TIME_IT_TOTAL) / iterations;
    }
  }
  cout << "Average: " << total << "ms" << endl;
}

auto benchmark_reklibpp_fastq() -> void {
  double total = 0;
  for (int i = 0; i < iterations + 1; ++i) {
    auto record = reklibpp::Seq(16 * 1024 * 1024);
    auto ssi = reklibpp::SeqStreamIn("benchmark_objects/FASTQ.fnq");
    TIME_IT({
      while (ssi >> record) { record.clear(); }
    })
    if (i > 0) {
      cout << "reklibpp fastq iteration " << i << ": " << TIME_IT_TOTAL << "ms"
           << endl;
      total += double(TIME_IT_TOTAL) / iterations;
    }
  }
  cout << "Average: " << total << "ms" << endl;
}

auto benchmark_klibpp_fasta() -> void {
  double total = 0;
  for (int i = 0; i < iterations + 1; ++i) {
    auto record = klibpp::KSeq();
    auto ssi = klibpp::SeqStreamIn("benchmark_objects/FASTA.fna");
    TIME_IT({
      while (ssi >> record) {}
    })
    if (i > 0) {
      cout << "klibpp fasta iteration " << i << ": " << TIME_IT_TOTAL << "ms"
           << endl;
      total += double(TIME_IT_TOTAL) / iterations;
    }
  }
  cout << "Average: " << total << "ms" << endl;
}

auto benchmark_klibpp_fastq() -> void {
  double total = 0;
  for (int i = 0; i < iterations + 1; ++i) {
    auto record = klibpp::KSeq();
    auto ssi = klibpp::SeqStreamIn("benchmark_objects/FASTQ.fnq");
    TIME_IT({
      while (ssi >> record) {}
    })
    if (i > 0) {
      cout << "klibpp fastq iteration " << i << ": " << TIME_IT_TOTAL << "ms"
           << endl;
      total += double(TIME_IT_TOTAL) / iterations;
    }
  }
  cout << "Average: " << total << "ms" << endl;
}

auto benchmark_reklibpp_fasta_gz() -> void {
  double total = 0;
  for (int i = 0; i < iterations + 1; ++i) {
    auto record = reklibpp::Seq(16 * 1024 * 1024);
    auto ssi = reklibpp::SeqStreamIn("benchmark_objects/FASTA.fna.gz");
    TIME_IT({
      while (ssi >> record) { record.clear(); }
    })
    if (i > 0) {
      cout << "reklibpp fasta.gz iteration " << i << ": " << TIME_IT_TOTAL << "ms"
           << endl;
      total += double(TIME_IT_TOTAL) / iterations;
    }
  }
  cout << "Average: " << total << "ms" << endl;
}

auto benchmark_reklibpp_fastq_gz() -> void {
  double total = 0;
  for (int i = 0; i < iterations + 1; ++i) {
    auto record = reklibpp::Seq(16 * 1024 * 1024);
    auto ssi = reklibpp::SeqStreamIn("benchmark_objects/FASTQ.fnq.gz");
    TIME_IT({
      while (ssi >> record) { record.clear(); }
    })
    if (i > 0) {
      cout << "reklibpp fastq.gz iteration " << i << ": " << TIME_IT_TOTAL << "ms"
           << endl;
      total += double(TIME_IT_TOTAL) / iterations;
    }
  }
  cout << "Average: " << total << "ms" << endl;
}

auto benchmark_klibpp_fasta_gz() -> void {
  double total = 0;
  for (int i = 0; i < iterations + 1; ++i) {
    auto record = klibpp::KSeq();
    auto ssi = klibpp::SeqStreamIn("benchmark_objects/FASTA.fna.gz");
    TIME_IT({
      while (ssi >> record) {}
    })
    if (i > 0) {
      cout << "klibpp fasta.gz iteration " << i << ": " << TIME_IT_TOTAL << "ms"
           << endl;
      total += double(TIME_IT_TOTAL) / iterations;
    }
  }
  cout << "Average: " << total << "ms" << endl;
}

auto benchmark_klibpp_fastq_gz() -> void {
  double total = 0;
  for (int i = 0; i < iterations + 1; ++i) {
    auto record = klibpp::KSeq();
    auto ssi = klibpp::SeqStreamIn("benchmark_objects/FASTQ.fnq.gz");
    TIME_IT({
      while (ssi >> record) {}
    })
    if (i > 0) {
      cout << "klibpp fastq.gz iteration " << i << ": " << TIME_IT_TOTAL << "ms"
           << endl;
      total += double(TIME_IT_TOTAL) / iterations;
    }
  }
  cout << "Average: " << total << "ms" << endl;
}
