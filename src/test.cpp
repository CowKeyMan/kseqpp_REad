#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "kseqpp_read.hpp"

using namespace reklibpp;

using std::max;
using std::string;
using std::vector;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

const vector<string> expected = { "1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC2",
                                  "3TCTAGCTACTACTACTGATGGATGGAATGTGATG4",
                                  "5TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6" };

class KseqppFullReader {
  public:
    vector<string> seqs = { "" };
    KseqppFullReader(const string filename, const size_t bufsize) {
      Seq record(bufsize);
      SeqStreamIn iss(filename.c_str());
      while (iss >> record) {
        size_t seq_start = 0;
        for (auto str_break: record.string_breaks) {
          seqs.back().append(
            record.seq.substr(seq_start, str_break - seq_start + 1)
          );
          seq_start = str_break + 1;
          if (seqs.back().size() > 0) { seqs.push_back(""); }
        }
        if (record.seq.size() - seq_start > 0) {
          seqs.back().append(
            record.seq.substr(seq_start, record.seq.size() - seq_start)
          );
        }
        record.clear();
      }
      if (seqs.back().size() == 0) { seqs.resize(seqs.size() - 1); }
    }

    void assert_correct() {
      for (size_t i = 0; i < max(seqs.size(), expected.size()); ++i) {
        ASSERT_EQ(seqs[i], expected[i]);
      }
    }
};

vector<Seq> get_seqs(const string filename, const size_t bufsize) {
  vector<Seq> ret;
  Seq record(bufsize);
  SeqStreamIn iss(filename.c_str());
  while (iss >> record) {
    ret.push_back(record);
    record.clear();
  }
  return ret;
}

const string fasta_file = "test_objects/queries.fna";
const string fastq_file = "test_objects/queries.fnq";

const vector<Seq> small_expected = {
  { "1ACTGCAATGGGCAAT", {} },    { "ATGTCTCTGTGTGGAT", {} },
  { "TAC23TCTAGCTACTA", { 3 } }, { "CTACTGATGGATGGAA", {} },
  { "TGTGATG45TGAGTGA", { 7 } }, { "GATGAGGTGATAGTGA", {} },
  { "CGTAGTGAGGA6", { 11 } },
};

const vector<Seq> half_expected
  = { { "1ACTGCAATGGGCAATAT", {} }, { "GTCTCTGTGTGGATTAC2", { 17 } },
      { "3TCTAGCTACTACTACTG", {} }, { "ATGGATGGAATGTGATG4", { 17 } },
      { "5TGAGTGAGATGAGGTGA", {} }, { "TAGTGACGTAGTGAGGA6", { 17 } } };

const vector<Seq> equal_expected
  = { { "1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC2", { 35 } },
      { "3TCTAGCTACTACTACTGATGGATGGAATGTGATG4", { 35 } },
      { "5TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6", { 35 } } };

const vector<Seq> big_expected
  = { { "1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCT", { 35 } },
      { "ACTACTACTGATGGATGGAATGTGATG45TGAGTGAGATGAGGT", { 27 } },
      { "GATAGTGACGTAGTGAGGA6", { 19 } } };

const vector<Seq> common_multiple_expected
  = { { "1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTG", { 35 } },
      { "ATGGATGGAATGTGATG45TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6",
        { 17, 53 } } };

const vector<Seq> full_expected
  = { { "1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTGATGGATGGAATGTGAT"
        "G45TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6",
        { 36 - 1, 36 * 2 - 1, 36 * 3 - 1 } } };

TEST(TestFASTA, SmallBufferFASTA) {
  KseqppFullReader(fasta_file, 16).assert_correct();
  ASSERT_EQ(get_seqs(fasta_file, 16), small_expected);
}

TEST(TestFASTA, HalfBufferFASTA) {
  KseqppFullReader(fasta_file, 18).assert_correct();
  ASSERT_EQ(get_seqs(fasta_file, 18), half_expected);
}

TEST(TestFASTA, EqualBufferFASTA) {
  KseqppFullReader(fasta_file, 36).assert_correct();
  ASSERT_EQ(get_seqs(fasta_file, 36), equal_expected);
}

TEST(TestFASTA, BigBufferFASTA) {
  KseqppFullReader(fasta_file, 44).assert_correct();
  ASSERT_EQ(get_seqs(fasta_file, 44), big_expected);
}

TEST(TestFASTA, CommonMultipleBufferFASTA) {
  // 36 * 1.5 = 54
  KseqppFullReader(fasta_file, 36 * 1.5).assert_correct();
  ASSERT_EQ(get_seqs(fasta_file, 36 * 1.5), common_multiple_expected);
}

TEST(TestFASTA, BiggestBufferFASTA) {
  KseqppFullReader(fasta_file, 9999).assert_correct();
  ASSERT_EQ(get_seqs(fasta_file, 9999), full_expected);
}

TEST(TestFASTQ, HalfBufferFASTQ) {
  KseqppFullReader(fastq_file, 18).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 44), half_expected);
}

TEST(TestFASTQ, SmallBufferFASTQ) {
  KseqppFullReader(fastq_file, 30).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 44), small_expected);
}

TEST(TestFASTQ, EqualBufferFASTQ) {
  KseqppFullReader(fastq_file, 36).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 44), equal_expected);
}

TEST(TestFASTQ, BigBufferFASTQ) {
  KseqppFullReader(fastq_file, 44).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 44), big_expected);
}

TEST(TestFASTQ, CommonMultipleBufferFASTQ) {
  // 36 * 1.5 = 54
  KseqppFullReader(fastq_file, 36 * 1.5).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 36 * 1.5), common_multiple_expected);
}

TEST(TestFASTQ, BiggestBufferFASTQ) {
  KseqppFullReader(fastq_file, 9999).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 9999), full_expected);
}
