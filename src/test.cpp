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

TEST(TestFASTA, SmallBufferFASTA) {
  KseqppFullReader("test_objects/queries.fna", 16).assert_correct();
}

TEST(TestFASTA, HalfBufferFASTA) {
  KseqppFullReader("test_objects/queries.fna", 18).assert_correct();
}

TEST(TestFASTA, EqualBufferFASTA) {
  KseqppFullReader("test_objects/queries.fna", 36).assert_correct();
}

TEST(TestFASTA, BigBufferFASTA) {
  KseqppFullReader("test_objects/queries.fna", 36 * 1.2).assert_correct();
}

TEST(TestFASTA, CommonMultipleBufferFASTA) {
  KseqppFullReader("test_objects/queries.fna", 36 * 1.5).assert_correct();
}

TEST(TestFASTA, BiggestBufferFASTA) {
  KseqppFullReader("test_objects/queries.fna", 9999).assert_correct();
}

TEST(TestFASTQ, HalfBufferFASTQ) {
  KseqppFullReader("test_objects/queries.fnq", 18).assert_correct();
}

TEST(TestFASTQ, SmallBufferFASTQ) {
  KseqppFullReader("test_objects/queries.fnq", 30).assert_correct();
}

TEST(TestFASTQ, EqualBufferFASTQ) {
  KseqppFullReader("test_objects/queries.fnq", 36).assert_correct();
}

TEST(TestFASTQ, BigBufferFASTQ) {
  KseqppFullReader("test_objects/queries.fnq", 36 * 1.2).assert_correct();
}

TEST(TestFASTQ, CommonMultipleBufferFASTQ) {
  KseqppFullReader("test_objects/queries.fnq", 36 * 1.5).assert_correct();
}

TEST(TestFASTQ, BiggestBufferFASTQ) {
  KseqppFullReader("test_objects/queries.fnq", 9999).assert_correct();
}
