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
  size_t buf_begin = 0;
  size_t buf_end = 0;
  size_t current_seq_size = 0;
  bool eof = false;
  bool finished_reading_seq = true;
  char next_char = 0;
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
      close_func(close_func_) {}

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
    while (
      !(this->eof || this->rec->seq.size() == rec->max_chars
        || this->rec->chars_before_new_read.size() == rec->max_reads)
    ) {
      // skip header
      if (this->finished_reading_seq) {
        this->current_seq_size = 0;
        this->finished_reading_seq = false;
        skip_to_next_line();
      }

      // populate read
      read_read();
      peek_next_char();
      if (!this->finished_reading_seq) { return true; }
      if (this->eof) { break; }
    }
    return rec->seq.size() + rec->chars_before_new_read.size()
      > initial_rec_size;
  }

  inline auto read_read() -> void {
    char c = 0;
    // for each line
    while (rec->seq.size() != rec->max_chars) {
      // check if this lines starts with a shit character
      c = peek_next_char();
      if (this->eof || c == '+' || c == '>' || c == '@') {
        rec->chars_before_new_read.push_back(rec->seq.size());
        this->finished_reading_seq = true;
        return;
      }

      fill_read();
      // 3 stopping conditions: buf_end, eof, max_chars, newline
      c = peek_next_char();
      if (c == '\r' || c == '\n') {
        skip_to_next_line();
        c = peek_next_char();
      }

      if (this->eof || c == '+' || c == '@' || c == '>') {
        if (c == '+') {
          skip_to_next_line();
          read_quality_string();
        }
        rec->chars_before_new_read.push_back(rec->seq.size());
        this->finished_reading_seq = true;
        break;
      }
    }
  }

  inline auto fill_read() {
    size_t seq_start = rec->seq.size();
    size_t seq_size = seq_start;
    size_t start = buf_begin;
    char c = 0;
    // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-pointer-arithmetic)
    while (this->buf_begin < buf_end && buf[this->buf_begin] != '\r'
           && buf[this->buf_begin] != '\n' && seq_size < rec->max_chars) {
      ++this->buf_begin;
      ++seq_size;
    }
    rec->seq.resize(rec->seq.size() + this->buf_begin - start);
    // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-pointer-arithmetic)
    std::copy(buf + start, buf + this->buf_begin, rec->seq.data() + seq_start);
    this->current_seq_size += this->buf_begin - start;
  }

  inline auto read_quality_string() { read_n_chars(this->current_seq_size); }

  /* Low-level methods */
  inline auto getc() noexcept -> char {
    if (this->buf_begin >= this->buf_end) {
      this->fetch_buffer();
      if (this->buf_end <= 0) { this->eof = true; }
      if (this->eof) { return 0; }
    }
    return this->buf[this->buf_begin++];
  }

  inline auto peek_next_char() noexcept -> char {
    char c = getc();
    --buf_begin;
    return c;
  }

  inline auto fetch_buffer() noexcept -> void {
    this->buf_begin = 0;
    this->buf_end = this->load_buf(this->file_handle, this->buf, this->bufsize);
  }

  inline auto skip_to_next_line() -> void {
    // stop once you find a newline, next getc() will be the next char
    while (!this->eof && getc() != '\n') {}
  }

  inline auto read_n_chars(size_t n) -> void {
    char c = 0;
    for (size_t i = 0; (c = getc()) && i < n; ++i) {
      while (c == '\r' || c == '\n') { c = getc(); }
    }
    skip_to_next_line();
  }

  inline auto get_next_char() -> char {
    char c = 0;
    while ((c = getc()) && (c == '\r' || c == '\n')) {}
    return c;
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
