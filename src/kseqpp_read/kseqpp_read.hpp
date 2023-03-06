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
  vector<size_t> chars_before_new_read;
  vector<char> seq;

  Seq(size_t max_chars_): max_chars(max_chars_) { seq.reserve(max_chars); }
  Seq(
    const vector<char> seq_,
    const vector<size_t> chars_before_new_read_,
    size_t max_chars_ = DEFAULT_BUFSIZE
  ):
      max_chars(max_chars_),
      seq(std::move(seq_)),
      chars_before_new_read(chars_before_new_read_) {
    seq.reserve(max_chars);
  }
  Seq(
    const std::string &seq_,
    const vector<size_t> chars_before_new_read_,
    size_t max_chars_ = DEFAULT_BUFSIZE
  ):
      max_chars(max_chars_), chars_before_new_read(chars_before_new_read_) {
    seq.reserve(max_chars);
    seq.resize(seq_.size());
    std::copy(seq_.begin(), seq_.end(), seq.begin());
  }

  inline void clear() {
    chars_before_new_read.clear();
    seq.clear();
  }

  bool operator==(const Seq &other) const {
    return seq == other.seq
      && chars_before_new_read == other.chars_before_new_read;
  }
};

template <typename TFile, typename TFunc>
class KStream {  // kstream_t
public:
  /* Typedefs */
  using close_type = int (*)(TFile);

protected:
  /* Separators */
  enum class SEP {
    SPACE = 0,  // isspace(): \t, \n, \v, \f, \r
    LINE = 2    // line separator: "\n" (Unix) or "\r\n" (Windows)
  };
  /* Data members */
  char *buf;               /**< @brief character buffer */
  size_t bufsize;          /**< @brief buffer size */
  size_t buf_begin;        /**< @brief begin buffer index */
  size_t buf_end;          /**< @brief end buffer index or error flag if -1 */
  size_t current_seq_size; /**< @brief size of the sequence being read */
  bool is_eof;             /**< @brief eof flag */
  bool is_tqs;             /**< @brief truncated quality string flag */
  bool first_header_char_read; /**< @brief next record ready flag */
  bool finished_reading_seq;   /**< @brief Sequence is done reading */
  bool has_read_something;     /**< @brief Current seq contains something */
  char next_char;    /** @brief A character to add at start of next sequence */
  size_t qual_size;  /**< @brief quality string size */
  TFile file_handle; /**< @brief file handler */
  TFunc load_buf;    /**< @brief read function */
  close_type close_func; /**< @brief close function */
public:
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
    this->is_eof = false;
    this->is_tqs = false;
    this->first_header_char_read = false;
    this->finished_reading_seq = true;
    this->current_seq_size = 0;
    this->next_char = 0;
  }

  KStream(TFile file_handle_, TFunc load_buf_, close_type close_func_):
      KStream(
        std::move(file_handle_),
        std::move(load_buf_),
        DEFAULT_BUFSIZE,
        close_func_
      ) {}

  KStream(KStream const &) = delete;
  KStream &operator=(KStream const &) = delete;

  KStream(KStream &&other) noexcept {
    this->buf = other.buf;
    other.buf = nullptr;
    this->bufsize = other.bufsize;
    this->buf_begin = other.buf_begin;
    this->buf_end = other.buf_end;
    this->is_eof = other.is_eof;
    this->is_tqs = other.is_tqs;
    this->first_header_char_read = other.first_header_char_read;
    this->file_handle = std::move(other.file_handle);
    this->load_buf = std::move(other.load_buf);
    this->close_func = other.close_func;
    this->finished_reading_seq = other.finished_reading_seq;
    this->current_seq_size = other.current_seq_size;
    this->next_char = other.next_char;
  }

  KStream &operator=(KStream &&other) noexcept {
    if (this == &other) { return *this; }
    delete[] this->buf;
    this->buf = other.buf;
    other.buf = nullptr;
    this->bufsize = other.bufsize;
    this->buf_begin = other.buf_begin;
    this->buf_end = other.buf_end;
    this->is_eof = other.is_eof;
    this->is_tqs = other.is_tqs;
    this->first_header_char_read = other.first_header_char_read;
    this->file_handle = std::move(other.file_handle);
    this->load_buf = std::move(other.load_buf);
    this->close_func = other.close_func;
    this->finished_reading_seq = other.finished_reading_seq;
    this->current_seq_size = other.current_seq_size;
    this->next_char = other.next_char;
    return *this;
  }

  ~KStream() noexcept {
    delete[] this->buf;
    if (this->close_func != nullptr) { this->close_func(this->file_handle); }
  }
  /* Methods */
  inline bool err() const {  // ks_err
    return this->buf_end == size_t(-1);
  }

  inline bool eof() const {  // ks_eof
    return this->is_eof && this->buf_begin >= this->buf_end;
  }

  inline bool tqs() const { return this->is_tqs; }

  inline bool fail() const {
    return this->err() || this->tqs() || (this->eof() && !has_read_something);
  }

  inline KStream &operator>>(Seq &rec) {
    char c = -1;
    this->has_read_something = false;
    while (true) {
      if (this->finished_reading_seq) { this->read_header(); }
      if (this->fail() || this->eof()) { return *this; }

      if (read_sequence(rec, c) || this->eof()) { return *this; }

      // Read + line if FASTQ
      if (c != '+') {
        continue;  // FASTA if no +
      }

      // Read Quality string line
      this->read_quality_string();
      if (this->is_tqs || this->err()) { return *this; }
      this->first_header_char_read = false;

      if (this->eof()) { return *this; }
    }
  }

  inline bool read_sequence(Seq &rec, char &last_char) {
    char c;
    size_t previous_length = rec.seq.size();
    while ((c = this->next_char ? this->next_char : this->getc())) {
      if (c == '\n' || c == '\r') {
        next_char = 0;
        continue;
      } else if (is_ready_char(c) || c == '+') {
        break;
      }
      if (rec.seq.size() == rec.max_chars) {
        this->next_char = c;
        this->current_seq_size += rec.seq.size() - previous_length;
        return true;
      } else {
        rec.seq.push_back(c);
        this->next_char = 0;
      }
      this->getuntil(
        SEP::LINE, &rec.seq, &this->next_char, nullptr, rec.max_chars
      );
      this->has_read_something = rec.seq.size() > 0;
    }
    this->current_seq_size += rec.seq.size() - previous_length;
    this->first_header_char_read = is_ready_char(c);
    this->finished_reading_seq = is_ready_char(c) || c == '+' || this->eof();
    last_char = c;
    if (this->finished_reading_seq) {
      rec.chars_before_new_read.push_back(rec.seq.size());
    }
    return false;
  }

  inline void read_header() {
    char last_char;
    this->finished_reading_seq = false;
    if (!this->first_header_char_read) {
      this->read_until_next_ready_char();
      if (this->fail()) { return; }
      this->first_header_char_read = true;
    }
    if (!this->read_name(last_char)) { return; }
    if (last_char != '\n') { this->read_comment(); }
    this->current_seq_size = 0;
  }

  inline void read_until_next_ready_char() {
    char c;
    while ((c = this->getc()) && !is_ready_char(c)) {}
  }

  inline bool read_name(char &last_char) {
    return this->getuntil(SEP::SPACE, nullptr, &last_char);
  }

  inline bool read_comment() {
    return this->getuntil(SEP::LINE, nullptr, nullptr);
  }

  inline bool is_ready_char(char c) { return c == '>' || c == '@'; }

  inline void skip_to_next_line() {
    char c;
    while ((c = this->getc()) && c != '\n') {}
  }

  inline void read_quality_string() {
    if (this->eof()) {  // no quality string
      this->is_tqs = true;
      return;
    }
    size_t qual_size = 0;
    while (qual_size < this->current_seq_size
           && this->getuntil(SEP::LINE, nullptr, nullptr, &qual_size)) {}
    if (this->err()) { return; }
    if (this->current_seq_size != qual_size) { this->is_tqs = true; }
  }

  operator bool() const { return !this->fail(); }

  /* Low-level methods */
  inline char getc() noexcept {
    if (this->err() || this->eof()) { return 0; }
    if (this->buf_begin >= this->buf_end) {
      this->fetch_buffer();
      this->check_eof();
      if (this->is_eof) { return 0; }
    }
    return this->buf[this->buf_begin++];
  }

  inline void fetch_buffer() noexcept {
    this->buf_begin = 0;
    this->buf_end = this->load_buf(this->file_handle, this->buf, this->bufsize);
  }

  inline void check_eof() {
    if (this->buf_end <= 0) { this->is_eof = true; }
  }

  inline bool getuntil(
    const SEP delimiter,
    vector<char> *str,
    char *last_char,
    size_t *str_size = nullptr,
    const size_t max_str_size = ULLONG_MAX
  ) noexcept {
    bool gotany = false;
    if (last_char != nullptr) { *last_char = 0; }
    size_t i = ULLONG_MAX;
    char c = -1;
    do {
      if (!this->getc()) { break; }
      --this->buf_begin;
      i = this->get_idx_at_delim_or_bufend(delimiter, str, max_str_size);
      gotany = true;
      this->append_to_string(str, str_size, i - this->buf_begin);
      this->buf_begin = i + 1;
      if (i > 0) {
        if (delimiter == SEP::LINE && this->buf[i - 1] == '\r') {
          remove_last_element(str, str_size);
        }
      }
    } while (i >= this->buf_end);
    if (this->err() || (this->eof() && !gotany)) { return false; }
    if (!this->eof() && last_char != nullptr) { *last_char = this->buf[i]; }
    return true;
  }

  inline size_t get_idx_at_delim_or_bufend(
    SEP delimiter, vector<char> *str, const size_t max_str_size
  ) {
    size_t i;
    if (delimiter == SEP::LINE) {
      for (i = this->buf_begin; i < this->buf_end && this->buf[i] != '\n'
           && (str == nullptr
               || str->size() + (i - this->buf_begin) < max_str_size);
           ++i) {}
    } else {  // delimiter == SEP::SPACE
      for (i = this->buf_begin;
           i < this->buf_end && !std::isspace(this->buf[i]);
           ++i) {}
    }
    return i;
  }

  inline void
  append_to_string(vector<char> *str, size_t *str_size, const size_t amount) {
    if (str != nullptr) {
      const size_t str_size = str->size();
      str->resize(str_size + amount);
      std::copy(
        this->buf + this->buf_begin,
        this->buf + this->buf_begin + amount,
        str->begin() + str_size
      );
    }
    if (str_size != nullptr) { *str_size += amount; }
  }

  inline void remove_last_element(vector<char> *str, size_t *str_size) {
    if (str != nullptr) { str->pop_back(); }
    if (str_size != nullptr) { *str_size -= 1; }
  }
};

class SeqStreamIn:
    public KStream<gzFile, int (*)(gzFile_s *, void *, unsigned int)> {
public:
  using base_type = KStream<gzFile, int (*)(gzFile_s *, void *, unsigned int)>;

  SeqStreamIn(const char *filename, const size_t bufsize = DEFAULT_BUFSIZE):
      base_type(gzopen(filename, "r"), gzread, bufsize, gzclose) {}
  SeqStreamIn(int fd, const size_t bufsize = DEFAULT_BUFSIZE):
      base_type(gzdopen(fd, "r"), gzread, bufsize, gzclose) {}
};

}  // namespace reklibpp
#endif
