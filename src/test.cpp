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

const vector<string> expected = {
  "1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC2",
  "3TCTAGCTACTACTACTGATGGATGGAATGTGATG4",
  "5TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6"};

class KseqppFullReader {
public:
  vector<string> seqs = {""};
  KseqppFullReader(
    const string filename,
    const size_t bufsize,
    const size_t file_bufsize = DEFAULT_BUFSIZE
  ) {
    Seq record(bufsize);
    SeqStreamIn iss(filename.c_str(), file_bufsize);
    while (iss >> record) {
      size_t seq_start = 0;
      for (auto str_break : record.chars_before_new_read) {
        const auto str_size = seqs.back().size();
        seqs.back().resize(str_size + str_break - seq_start);
        std::copy(
          record.seq.begin() + seq_start,
          record.seq.begin() + str_break,
          seqs.back().begin() + str_size
        );
        seq_start = str_break;
        if (seqs.back().size() > 0) { seqs.push_back(""); }
      }
      if (record.seq.size() - seq_start > 0) {
        const auto str_size = seqs.back().size();
        seqs.back().resize(str_size + record.seq.size() - seq_start);
        std::copy(
          record.seq.begin() + seq_start,
          record.seq.begin() + record.seq.size(),
          seqs.back().begin() + str_size
        );
        /* seqs.back().append( */
        /*   record.seq.substr(seq_start, record.seq.size() - seq_start) */
        /* ); */
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

vector<Seq> get_seqs(
  const string filename,
  const size_t bufsize,
  const size_t file_bufsize = DEFAULT_BUFSIZE
) {
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
  {"1ACTGCAATGGGCAAT", {}},
  {"ATGTCTCTGTGTGGAT", {}},
  {"TAC23TCTAGCTACTA", {4}},
  {"CTACTGATGGATGGAA", {}},
  {"TGTGATG45TGAGTGA", {8}},
  {"GATGAGGTGATAGTGA", {}},
  {"CGTAGTGAGGA6", {12}},
};

const vector<Seq> half_expected = {
  {"1ACTGCAATGGGCAATAT", {}},
  {"GTCTCTGTGTGGATTAC2", {18}},
  {"3TCTAGCTACTACTACTG", {}},
  {"ATGGATGGAATGTGATG4", {18}},
  {"5TGAGTGAGATGAGGTGA", {}},
  {"TAGTGACGTAGTGAGGA6", {18}}};

const vector<Seq> equal_expected = {
  {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC2", {36}},
  {"3TCTAGCTACTACTACTGATGGATGGAATGTGATG4", {36}},
  {"5TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6", {36}}};

const vector<Seq> big_expected = {
  {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCT", {36}},
  {"ACTACTACTGATGGATGGAATGTGATG45TGAGTGAGATGAGGT", {28}},
  {"GATAGTGACGTAGTGAGGA6", {20}}};

const vector<Seq> common_multiple_expected = {
  {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTG", {36}},
  {"ATGGATGGAATGTGATG45TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6", {18, 54}}};

const vector<Seq> full_expected = {
  {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTGATGGATGGAATGTGAT"
   "G45TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6",
   {36, 36 * 2, 36 * 3}}};

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

TEST(TestFASTQ, SmallBufferFASTQ) {
  KseqppFullReader(fastq_file, 16).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 16), small_expected);
}

TEST(TestFASTQ, HalfBufferFASTQ) {
  KseqppFullReader(fastq_file, 18).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 18), half_expected);
}

TEST(TestFASTQ, EqualBufferFASTQ) {
  KseqppFullReader(fastq_file, 36).assert_correct();
  ASSERT_EQ(get_seqs(fastq_file, 36), equal_expected);
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

const vector<Seq> full_expected_empty_line = {
  {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTGATGGATGGAATGTGAT"
   "G45TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6",
   {36, 36, 36 * 2, 36 * 3}}};

TEST(TestEmptyLine, Test) {
  KseqppFullReader("test_objects/fasta_empty_line.fna", 9999).assert_correct();
  ASSERT_EQ(
    get_seqs("test_objects/fasta_empty_line.fna", 9999),
    full_expected_empty_line
  );
}

TEST(TestFASTA, CommonMultipleBufferFASTASmallFileBuffer) {
  // 36 * 1.5 = 54
  KseqppFullReader(fasta_file, 36 * 1.5, 36).assert_correct();
  ASSERT_EQ(get_seqs(fasta_file, 36 * 1.5, 36), common_multiple_expected);
}

TEST(TestFASTA, BiggestBufferFASTASmallFileBuffer) {
  KseqppFullReader(fasta_file, 9999, 36).assert_correct();
  ASSERT_EQ(get_seqs(fasta_file, 9999, 36), full_expected);
}
