#ifndef KSEQPP_READ_HPP
#define KSEQPP_READ_HPP

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <ios>
#include <string>
#include <vector>
#include <zlib.h>

using std::string;
using std::vector;

namespace reklibpp {

using size_type = uint64_t;

class Seq {  // kseq_t
  public:
    Seq(size_type max_seq_size) {
      this->max_seq_size = max_seq_size;
      this->seq.reserve(max_seq_size);
    }
    size_type max_seq_size;
    vector<size_type> string_breaks;
    string seq;
    size_type qual_size = 0;

    inline void clear() {
      string_breaks.resize(0);
      seq.resize(0);
      qual_size = 0;
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
      LINE = 2  // line separator: "\n" (Unix) or "\r\n" (Windows)
    };
    /* Consts */
    constexpr static std::make_unsigned_t<size_type> DEFAULT_BUFSIZE
      = 16384;  // 16kB buffer
    /* Data members */
    char *buf; /**< @brief character buffer */
    size_type bufsize; /**< @brief buffer size */
    size_type begin; /**< @brief begin buffer index */
    size_type buf_end; /**< @brief end buffer index or error flag if -1 */
    size_type current_seq_size; /**< @brief size of the sequence being read */
    bool is_eof; /**< @brief eof flag */
    bool is_tqs; /**< @brief truncated quality string flag */
    bool is_ready; /**< @brief next record ready flag */
    bool finished_reading_seq; /**< @brief Sequence is done reading */
    bool has_read_something; /**< @brief Current seq contains something */
    TFile f; /**< @brief file handler */
    TFunc func; /**< @brief read function */
    close_type close; /**< @brief close function */
  public:
    KStream(
      TFile f_,
      TFunc func_,
      std::make_unsigned_t<size_type> bs_ = DEFAULT_BUFSIZE,
      close_type cfunc_ = nullptr
    )  // ks_init
        :
        buf(new char[bs_]),
        bufsize(bs_),
        f(std::move(f_)),
        func(std::move(func_)),
        close(cfunc_) {
      this->begin = 0;
      this->buf_end = 0;
      this->is_eof = false;
      this->is_tqs = false;
      this->is_ready = false;
      this->finished_reading_seq = true;
      this->current_seq_size = 0;
    }

    KStream(TFile f_, TFunc func_, close_type cfunc_):
        KStream(std::move(f_), std::move(func_), DEFAULT_BUFSIZE, cfunc_) {}

    KStream(KStream const &) = delete;
    KStream &operator=(KStream const &) = delete;

    KStream(KStream &&other) noexcept {
      this->buf = other.buf;
      other.buf = nullptr;
      this->bufsize = other.bufsize;
      this->begin = other.begin;
      this->buf_end = other.buf_end;
      this->is_eof = other.is_eof;
      this->is_tqs = other.is_tqs;
      this->is_ready = other.is_ready;
      this->f = std::move(other.f);
      this->func = std::move(other.func);
      this->close = other.close;
      this->finished_reading_seq = other.finished_reading_seq;
      this->current_seq_size = other.current_seq_size;
    }

    KStream &operator=(KStream &&other) noexcept {
      if (this == &other) return *this;
      delete[] this->buf;
      this->buf = other.buf;
      other.buf = nullptr;
      this->bufsize = other.bufsize;
      this->begin = other.begin;
      this->buf_end = other.buf_end;
      this->is_eof = other.is_eof;
      this->is_tqs = other.is_tqs;
      this->is_ready = other.is_ready;
      this->f = std::move(other.f);
      this->func = std::move(other.func);
      this->close = other.close;
      this->finished_reading_seq = other.finished_reading_seq;
      this->current_seq_size = other.current_seq_size;
      return *this;
    }

    ~KStream() noexcept {
      delete[] this->buf;
      if (this->close != nullptr) this->close(this->f);
    }
    /* Methods */
    inline bool err() const {  // ks_err
      return this->buf_end == -1;
    }

    inline bool eof() const {  // ks_eof
      return this->is_eof && this->begin >= this->buf_end;
    }

    inline bool tqs() const { return this->is_tqs; }

    inline bool fail() const {
      return this->err() || this->tqs()
          || (this->eof() && !has_read_something);
    }

    inline KStream &operator>>(Seq &rec) {
      char c;
      this->has_read_something = false;
      while (true) {
        if (this->finished_reading_seq) {
          this->current_seq_size = 0;
          this->finished_reading_seq = false;
          if (!this->is_ready) {  // then jump to the next header line
            while ((c = this->getc()) && c != '>' && c != '@') {}
            if (this->fail()) return *this;
            this->is_ready = true;
          }  // else: the first header char has been read in the previous call
          if (!this->getuntil(SEP::SPACE, nullptr, &c))
            return *this;  // read name
          if (c != '\n') {
            this->getuntil(SEP::LINE, nullptr, nullptr);  // read comment
          }
        }

        // read Sequence
        size_type previous_length = rec.seq.size();
        while (rec.seq.size() < rec.max_seq_size - 1 && (c = this->getc())
               && c != '>' && c != '@' && c != '+') {
          if (c == '\n') continue;  // skip empty lines
          rec.seq += c;
          // read the rest of the line
          this->getuntil(
            SEP::LINE, &rec.seq, &c, nullptr, rec.max_seq_size - 1
          );
          this->has_read_something = true;
        }
        this->finished_reading_seq
          = c == '>' || c == '@' || c == '+' || std::isspace(c) || this->eof();
        this->current_seq_size += rec.seq.size() - previous_length;
        if (this->finished_reading_seq) {
          rec.string_breaks.push_back(rec.seq.size() - 1);
        } else {
          rec.seq += c;
          ++this->current_seq_size;
          return *this;
        }
        // the first header character has been read
        if (c == '>' || c == '@') this->is_ready = true;

        // Read + line if FASTQ
        if (c != '+') return *this;  // FASTA
        while ((c = this->getc()) && c != '\n') {}  // skip the rest of '+' line

        // Read Quality string line
        if (this->eof()) {  // error: no quality string
          this->is_tqs = true;
          return *this;
        }
        while (this->getuntil(SEP::LINE, nullptr, nullptr, &rec.qual_size)
               && rec.qual_size < this->current_seq_size) {}
        if (this->err()) return *this;
        this->is_ready = false;  // we have not come to the next header line
        if (this->current_seq_size != rec.qual_size) {
          this->is_tqs = true;  // error: qual string is of a different length
          return *this;
        }
        rec.qual_size = 0;
        if (this->eof()) { return *this; }
      }
    }

    operator bool() const { return !this->fail(); }

    /* Low-level methods */
    inline char getc() noexcept {  // ks_getc
      if (this->err() || this->eof()) return 0;  // error
      if (this->begin >= this->buf_end) {  // fetch
        this->begin = 0;
        this->buf_end = this->func(this->f, this->buf, this->bufsize);
        if (this->buf_end <= 0) {  // err if end == -1 and eof if 0
          this->is_eof = true;
          return 0;
        }
      }
      // ready
      return this->buf[this->begin++];
    }

    inline bool getuntil(
      const SEP delimiter,
      std::string *str,
      char *last_char,
      size_type *str_size = nullptr,
      const size_type max_str_size = -1  // max value
    ) noexcept {
      bool gotany = false;
      if (last_char != nullptr) *last_char = 0;
      size_type i = -1;
      char c = -1;
      do {
        if (!(c = this->getc())) break;
        --this->begin;
        // find index i where there is the delimeter or until the buffer's end
        if (delimiter == SEP::LINE) {
          for (i = this->begin;
               i < this->buf_end && this->buf[i] != '\n'
               && (str == nullptr || str->size() + (i - this->begin) < max_str_size);
               ++i) {}
        } else {  // delimiter == SEP::SPACE
          for (i = this->begin;
               i < this->buf_end && !std::isspace(this->buf[i]);
               ++i) {}
        }
        gotany = true;
        if (str != nullptr) {
          str->append(this->buf + this->begin, i - this->begin);
        }
        if (str_size != nullptr) {
          *str_size += i - this->begin;
        }
        if (i > 0) { c = this->buf[i - 1]; }
        this->begin = i + 1;
      } while (i >= this->buf_end);
      if (this->err() || (this->eof() && !gotany)) return false;
      assert(i != -1);
      if (!this->eof() && last_char != nullptr) *last_char = this->buf[i];
      if (delimiter == SEP::LINE && c == '\r') {
        if (str != nullptr) { str->pop_back(); }
        if (str_size != nullptr) { *str_size -= 1; }
      }
      return true;
    }
};

class SeqStreamIn:
    public KStream<gzFile, int (*)(gzFile_s *, void *, unsigned int)> {
  public:
    using base_type
      = KStream<gzFile, int (*)(gzFile_s *, void *, unsigned int)>;

    SeqStreamIn(const char *filename):
        base_type(gzopen(filename, "r"), gzread, gzclose) {}
    SeqStreamIn(int fd): base_type(gzdopen(fd, "r"), gzread, gzclose) {}
};

}  // namespace reklibpp
#endif
