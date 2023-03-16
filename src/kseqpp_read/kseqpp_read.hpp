#ifndef KSEQPP_READ_HPP
#define KSEQPP_READ_HPP

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <ios>
#include <stdexcept>
#include <string>
#include <vector>
#include <zlib.h>

namespace reklibpp {

using std::vector;
using size_t = uint64_t;

const size_t DEFAULT_BUFSIZE = 16ULL * 1024;

class Seq {  // kseq_t
public:
  size_t max_chars;
  size_t max_reads;
  vector<size_t> chars_before_new_read;
  vector<char> seq;

  explicit Seq(
    size_t max_chars_ = DEFAULT_BUFSIZE,
    size_t max_reads_ = DEFAULT_BUFSIZE / 100
  ):
      max_chars(max_chars_), max_reads(max_reads_) {
    seq.reserve(max_chars);
    chars_before_new_read.reserve(max_reads);
  }
  Seq(
    vector<char> seq_,
    vector<size_t> chars_before_new_read_,
    size_t max_chars_ = DEFAULT_BUFSIZE,
    size_t max_reads_ = DEFAULT_BUFSIZE / 100
  ):
      max_chars(max_chars_),
      seq(std::move(seq_)),
      chars_before_new_read(std::move(chars_before_new_read_)),
      max_reads(max_reads_) {
    seq.reserve(max_chars);
  }
  Seq(
    const std::string &seq_,
    vector<size_t> chars_before_new_read_,
    size_t max_chars_ = DEFAULT_BUFSIZE,
    size_t max_reads_ = DEFAULT_BUFSIZE / 100
  ):
      max_chars(max_chars_),
      chars_before_new_read(std::move(chars_before_new_read_)),
      max_reads(max_reads_) {
    seq.reserve(max_chars);
    seq.resize(seq_.size());
    std::copy(seq_.begin(), seq_.end(), seq.begin());
  }

  inline void clear() {
    chars_before_new_read.clear();
    seq.clear();
  }

  auto operator==(const Seq &other) const -> bool {
    return seq == other.seq
      && chars_before_new_read == other.chars_before_new_read;
  }
};

template <typename TFile, typename TFunc>
class KStream {  // kstream_t
public:
  /* Typedefs */
  using close_type = int (*)(TFile);

private:
  /* Separators */
  enum class SEP {
    SPACE = 0,  // isspace(): \t, \n, \v, \f, \r
    LINE = 2    // line separator: "\n" (Unix) or "\r\n" (Windows)
  };
  /* Data members */
  Seq *rec;
  char *buf;
  size_t bufsize;
  size_t buf_begin;
  size_t buf_end;
  size_t current_seq_size;
  bool eof;
  bool finished_reading_seq;
  char next_char;
  size_t qual_size;
  TFile file_handle;
  TFunc load_buf;
  close_type close_func;

public:
  // NOLINTNEXTLINE (cppcoreguidelines-pro-type-member-init,hicpp-member-init)
  KStream(
    TFile file_handle_,
    TFunc load_bufile_handle_,
    size_t _bufsize = DEFAULT_BUFSIZE,
    close_type close_func_ = nullptr
  )  // ks_init
      :
      buf(new char[_bufsize]),
      bufsize(_bufsize),
      file_handle(std::move(file_handle_)),
      load_buf(std::move(load_bufile_handle_)),
      close_func(close_func_) {
    this->buf_begin = 0;
    this->buf_end = 0;
    this->eof = false;
    this->finished_reading_seq = true;
    this->current_seq_size = true;
    this->next_char = 0;
  }

  // NOLINTNEXTLINE (cppcoreguidelines-pro-type-member-init,hicpp-member-init)
  KStream(TFile file_handle_, TFunc load_buf_, close_type close_func_):
      KStream(
        std::move(file_handle_),
        std::move(load_buf_),
        DEFAULT_BUFSIZE,
        close_func_
      ) {}

  KStream(KStream &) = delete;
  KStream(KStream &&other) = delete;
  auto operator=(KStream &) = delete;
  auto operator=(KStream &&) = delete;

  ~KStream() noexcept {
    delete[] this->buf;
    if (this->close_func != nullptr) { this->close_func(this->file_handle); }
  }

  inline auto operator>>(Seq &rec_) -> bool {
    rec = &rec_;
    size_t initial_rec_size
      = rec->seq.size() + rec->chars_before_new_read.size();
    while (!this->eof) {
      go_to_next_char();
      if (this->eof) { break; }
      // skip header
      if (this->finished_reading_seq) {
        this->current_seq_size = 0;
        read_header();
      }

      if (this->eof) { break; }
      read_read();
      if (this->eof) { break; }
      if (!this->finished_reading_seq) { return true; }
      read_quality_string();
      if (
        this->rec->seq.size() == rec->max_chars
        || this->rec->chars_before_new_read.size() == rec->max_reads
        || this->eof
      ) {
        break;
      }
    }
    return rec->seq.size() + rec->chars_before_new_read.size()
      > initial_rec_size;
  }

  inline auto read_header() {
    go_to_next_char();
    skip_to_next_line();
  }

  inline auto read_read() -> void {
    this->finished_reading_seq = false;
    char c = 0;
    while (true) {
      go_to_next_char();
      c = getc();
      if (this->eof || c == '+' || c == '>' || c == '@') {
        --buf_begin;
        rec->chars_before_new_read.push_back(rec->seq.size());
        this->finished_reading_seq = true;
        return;
      }
      rec->seq.push_back(c);
      ++this->current_seq_size;
      if (rec->seq.size() == rec->max_chars) {
        c = getc();
        while (c == '\r' || c == '\n') { c = getc(); }
        if (this->eof || c == '+' || c == '>' || c == '@') {
          rec->chars_before_new_read.push_back(rec->seq.size());
          this->finished_reading_seq = true;
        }
        --buf_begin;
        return;
      }
    }
  }

  inline auto read_quality_string() {
    char c = getc();
    if (c == '+') {
      skip_to_next_line();
      read_n_chars(this->current_seq_size);
    } else {
      --buf_begin;
    }
  }

  /* Low-level methods */
  inline auto getc() noexcept -> char {
    if (this->buf_begin >= this->buf_end) {
      this->fetch_buffer();
      if (this->buf_end <= 0) { this->eof = true; }
      if (this->eof) { return 0; }
    }
    return this->buf[this->buf_begin++];
  }

  inline auto fetch_buffer() noexcept -> void {
    this->buf_begin = 0;
    this->buf_end = this->load_buf(this->file_handle, this->buf, this->bufsize);
  }

  inline auto skip_to_next_line() -> void {
    while (!this->eof && getc() != '\n') {}
  }

  inline auto read_n_chars(size_t n) -> void {
    char c = 0;
    for (size_t i = 0; (c = getc()) && i < n; ++i) {
      while (c == '\r' || c == '\n') { c = getc(); }
    }
  }

  inline auto go_to_next_char() -> void {
    char c = 0;
    while ((c = getc()) && (c == '\r' || c == '\n')) {}
    --buf_begin;
  }
};

class SeqStreamIn:
    public KStream<gzFile, int (*)(gzFile_s *, void *, unsigned int)> {
public:
  using base_type = KStream<gzFile, int (*)(gzFile_s *, void *, unsigned int)>;

  explicit SeqStreamIn(
    const char *filename, const size_t bufsize = DEFAULT_BUFSIZE
  ):
      base_type(gzopen(filename, "r"), gzread, bufsize, gzclose) {}
  explicit SeqStreamIn(int fd, const size_t bufsize = DEFAULT_BUFSIZE):
      base_type(gzdopen(fd, "r"), gzread, bufsize, gzclose) {}
};

}  // namespace reklibpp
#endif
