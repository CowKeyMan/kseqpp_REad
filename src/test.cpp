#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "kseqpp_read.hpp"

namespace reklibpp {

using std::max;
using std::string;
using std::vector;

auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

class KseqppFullReader {
private:
  vector<string> expected = {
    "1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC2",
    "3TCTAGCTACTACTACTGATGGATGGAATGTGATG4",
    "5TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6"};
  vector<string> seqs = {""};

public:
  KseqppFullReader(
    const string &filename,
    const size_t max_chars,
    const size_t max_reads = DEFAULT_BUFSIZE / 100,
    const size_t file_bufsize = DEFAULT_BUFSIZE
  ) {
    Seq record(max_chars, max_reads);
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
        if (!seqs.back().empty()) { seqs.emplace_back(""); }
      }
      if (record.seq.size() - seq_start > 0) {
        const auto str_size = seqs.back().size();
        seqs.back().resize(str_size + record.seq.size() - seq_start);
        std::copy(
          record.seq.begin() + seq_start,
          record.seq.begin() + record.seq.size(),
          seqs.back().begin() + str_size
        );
      }
      record.clear();
    }
    if (seqs.back().empty()) { seqs.resize(seqs.size() - 1); }
  }

  void assert_correct() {
    for (size_t i = 0; i < max(seqs.size(), expected.size()); ++i) {
      ASSERT_EQ(seqs[i], expected[i]) << "differs at index " << i;
    }
  }
};

auto get_seqs(
  const string &filename,
  const size_t max_chars,
  const size_t max_reads = DEFAULT_BUFSIZE / 100,
  const size_t file_bufsize = DEFAULT_BUFSIZE
) -> vector<Seq> {
  vector<Seq> ret;
  Seq record(max_chars, max_reads);
  SeqStreamIn iss(filename.c_str(), file_bufsize);
  while (iss >> record) {
    ret.push_back(record);
    record.clear();
  }
  return ret;
}

class Test: public ::testing::Test {
public:
  string fasta_file = "test_objects/queries.fna";
  string fastq_file = "test_objects/queries.fnq";
  vector<Seq> small_expected = {
    {"1ACTGCAATGGGCAAT", {}},
    {"ATGTCTCTGTGTGGAT", {}},
    {"TAC23TCTAGCTACTA", {4}},
    {"CTACTGATGGATGGAA", {}},
    {"TGTGATG45TGAGTGA", {8}},
    {"GATGAGGTGATAGTGA", {}},
    {"CGTAGTGAGGA6", {12}},
  };
  vector<Seq> half_expected = {
    {"1ACTGCAATGGGCAATAT", {}},
    {"GTCTCTGTGTGGATTAC2", {18}},
    {"3TCTAGCTACTACTACTG", {}},
    {"ATGGATGGAATGTGATG4", {18}},
    {"5TGAGTGAGATGAGGTGA", {}},
    {"TAGTGACGTAGTGAGGA6", {18}}};
  vector<Seq> equal_expected = {
    {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC2", {36}},
    {"3TCTAGCTACTACTACTGATGGATGGAATGTGATG4", {36}},
    {"5TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6", {36}}};
  vector<Seq> big_expected = {
    {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCT", {36}},
    {"ACTACTACTGATGGATGGAATGTGATG45TGAGTGAGATGAGGT", {28}},
    {"GATAGTGACGTAGTGAGGA6", {20}}};
  vector<Seq> common_multiple_expected = {
    {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTG", {36}},
    {"ATGGATGGAATGTGATG45TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6", {18, 54}}};
  vector<Seq> full_expected = {
    {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTGATGGATGGAATGTGAT"
     "G45TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6",
     {36, 36ULL * 2, 36ULL * 3}}};
  vector<Seq> full_expected_empty_line = {
    {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTGATGGATGGAATGTGAT"
     "G45TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6",
     {36, 36, 36ULL * 2, 36ULL * 3}}};
  vector<Seq> double_expected = {
    {"1ACTGCAATGGGCAATATGTCTCTGTGTGGATTAC23TCTAGCTACTACTACTGATGGATGGAATGTGATG4",
     {36, 72}},
    {"5TGAGTGAGATGAGGTGATAGTGACGTAGTGAGGA6", {36}}};
};

auto assert_seqs_equal(const vector<Seq> &seqs1, const vector<Seq> &seqs2)
  -> void {
  for (size_t i = 0; i < seqs1.size(); ++i) {
    EXPECT_EQ(seqs1[i].chars_before_new_read, seqs2[i].chars_before_new_read)
      << "differs at index " << i;
    EXPECT_EQ(seqs1[i].seq, seqs2[i].seq) << "differs at index " << i;
  }
}

TEST_F(Test, TestFASTASmallBufferFASTA) {
  KseqppFullReader(fasta_file, 16).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 16), small_expected);
}

TEST_F(Test, TestFASTAHalfBufferFASTA) {
  KseqppFullReader(fasta_file, 18).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 18), half_expected);
}

TEST_F(Test, TestFASTAEqualBufferFASTA) {
  KseqppFullReader(fasta_file, 36).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 36), equal_expected);
}

TEST_F(Test, TestFASTABigBufferFASTA) {
  KseqppFullReader(fasta_file, 44).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 44), big_expected);
}

TEST_F(Test, TestFASTACommonMultipleBufferFASTA) {
  // 36 * 1.5 = 54
  KseqppFullReader(fasta_file, 36 * 1.5).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 36 * 1.5), common_multiple_expected);
}

TEST_F(Test, TestFASTABiggestBufferFASTA) {
  KseqppFullReader(fasta_file, 9999).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 9999), full_expected);
}

TEST_F(Test, TestFASTQSmallBufferFASTQ) {
  KseqppFullReader(fastq_file, 16).assert_correct();
  assert_seqs_equal(get_seqs(fastq_file, 16), small_expected);
}

TEST_F(Test, TestFASTQHalfBufferFASTQ) {
  KseqppFullReader(fastq_file, 18).assert_correct();
  assert_seqs_equal(get_seqs(fastq_file, 18), half_expected);
}

TEST_F(Test, TestFASTQEqualBufferFASTQ) {
  KseqppFullReader(fastq_file, 36).assert_correct();
  assert_seqs_equal(get_seqs(fastq_file, 36), equal_expected);
}

TEST_F(Test, TestFASTQBigBufferFASTQ) {
  KseqppFullReader(fastq_file, 44).assert_correct();
  assert_seqs_equal(get_seqs(fastq_file, 44), big_expected);
}

TEST_F(Test, TestFASTQCommonMultipleBufferFASTQ) {
  // 36 * 1.5 = 54
  KseqppFullReader(fastq_file, 36 * 1.5).assert_correct();
  assert_seqs_equal(get_seqs(fastq_file, 36 * 1.5), common_multiple_expected);
}

TEST_F(Test, TestFASTQBiggestBufferFASTQ) {
  KseqppFullReader(fastq_file, 9999).assert_correct();
  assert_seqs_equal(get_seqs(fastq_file, 9999), full_expected);
}

TEST_F(Test, TestEmptyLineTest) {
  KseqppFullReader("test_objects/fasta_empty_line.fna", 9999).assert_correct();
  assert_seqs_equal(
    get_seqs("test_objects/fasta_empty_line.fna", 9999),
    full_expected_empty_line
  );
}

TEST_F(Test, TestFASTACommonMultipleBufferFASTASmallFileBuffer) {
  // 36 * 1.5 = 54
  KseqppFullReader(fasta_file, 36 * 1.5, 999, 36).assert_correct();
  assert_seqs_equal(
    get_seqs(fasta_file, 36 * 1.5, 999, 36), common_multiple_expected
  );
}

TEST_F(Test, TestFASTABiggestBufferFASTASmallFileBuffer) {
  KseqppFullReader(fasta_file, 9999, 999, 36).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 9999, 999, 36), full_expected);
}

TEST_F(Test, TestLimitedMaxReads) {
  KseqppFullReader(fasta_file, 9999, 1).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 9999, 1), equal_expected);
}

TEST_F(Test, TestLimitedMaxReadsDouble) {
  KseqppFullReader(fasta_file, 9999, 2).assert_correct();
  assert_seqs_equal(get_seqs(fasta_file, 9999, 2), double_expected);
}

}  // namespace reklibpp
