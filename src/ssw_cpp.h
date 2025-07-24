// ssw_cpp.h
// Created by Wan-Ping Lee
// Last revision by Mengyao Zhao on 2023-Apr-21

#ifndef COMPLETE_STRIPED_SMITH_WATERMAN_CPP_H_
#define COMPLETE_STRIPED_SMITH_WATERMAN_CPP_H_

#include <stddef.h>
#include <stdint.h>
#include <string>
#include <vector>

namespace StripedSmithWaterman {

struct Alignment {
  uint16_t    sw_score;            // The best alignment score
  uint16_t    sw_score_next_best;  // The next best alignment score
  int32_t     ref_begin;           // Reference begin position of the best alignment
  int32_t     ref_end;             // Reference end position of the best alignment
  int32_t     query_begin;         // Query begin position of the best alignment
  int32_t     query_end;           // Query end position of the best alignment
  int32_t     ref_end_next_best;   // Reference end position of the next best alignment
  int32_t     mismatches;          // Number of mismatches of the alignment
  std::string cigar_string;        // Cigar string of the best alignment
                                   // Cigar stored in the BAM format
                                   //   high 28 bits: length
  std::vector<uint32_t> cigar;     //   low 4 bits: M/I/D/S/X (0/1/2/4/8);

  Alignment()
      : sw_score(0)
      , sw_score_next_best(0)
      , ref_begin(0)
      , ref_end(0)
      , query_begin(0)
      , query_end(0)
      , ref_end_next_best(0)
      , mismatches(0) {}
};

struct Filter {
  // NOTE: No matter the filter, those five fields of Alignment will be given anyway.
  //       sw_score; sw_score_next_best; ref_end; query_end; ref_end_next_best.
  // NOTE: Only need score of alignments, please set 'report_begin_position'
  //       and 'report_cigar' false.

  bool report_begin_position;  // Give ref_begin and query_begin.
                               //   If it is not set, ref_begin and query_begin are -1.
  bool report_cigar;           // Give cigar_string and cigar.
                               //   report_begin_position is automatically TRUE.

  // When *report_cigar* is true and alignment passes these two filters,
  //   cigar_string and cigar will be given.
  uint16_t score_filter;     // score >= score_filter
  uint16_t distance_filter;  // ((ref_end - ref_begin) < distance_filter) &&

  // ((query_end - read_begin) < distance_filter)

  Filter()
      : report_begin_position(true)
      , report_cigar(true)
      , score_filter(0)
      , distance_filter(32767) {}
};

class Aligner {
public:
  // =========
  // @function Construct an Aligner on default values.
  //             The function will build the {A,C,G,T,N} aligner.
  //             If you target for other character aligners, then please
  //             use the other constructor and pass the corresponding matrix in.
  // =========
  Aligner();

  // =========
  // @function Construct an Aligner by assigning scores.
  //             The function will build the {A,C,G,T,N} aligner.
  //             If you target for other character aligners, then please
  //             use the other constructor and pass the corresponding matrix in.
  // =========
  Aligner(
      uint8_t match_score,
      uint8_t mismatch_penalty,
      uint8_t gap_opening_penalty,
      uint8_t gap_extending_penalty);

  // =========
  // @function Construct an Aligner by the specific matrices.
  // =========
  Aligner(
      const int8_t* score_matrix,
      int           score_matrix_size,
      const int8_t* translation_matrix,
      int           translation_matrix_size);

  // =========
  // @function Build the reference sequence and thus make
  //             Align(const char* query, size_t query_len, const Filter&
  //             filter, Alignment& alignment, int32_t maskLen) function;
  //             otherwise the reference should be given when aligning.
  //           [NOTICE] If there exists a sequence, that one will be deleted
  //                    and replaced.
  // @param    ref     The reference bases;
  //                   [NOTICE] Does not have to be null terminated.
  // @param    ref_len The number of bases in ref.
  // @return   The length of the built bases.
  // =========
  size_t SetReferenceSequence(const char* ref, size_t ref_len);

  // =========
  // @function Like SetReferenceSequence(const char* ref, size_t ref_len), but
  // requires null-terminated ref.
  // =========
  size_t SetReferenceSequence(const char* ref);

  void ClearReferenceSequence();

  // =========
  // @function Set penalties for opening and extending gaps
  //           [NOTICE] The defaults are 3 and 1 respectively.
  // =========
  void SetGapPenalty(uint8_t opening, uint8_t extending);

  // =========
  // @function Align the query against the reference set by SetReferenceSequence.
  // @param    query     The query sequence.
  //                     [NOTICE] Does not have to be null terminated.
  // @param    query_len The number of bases in query.
  // @param    filter    The filter for the alignment.
  // @param    alignment The container contains the result.
  // @param    maskLen   The distance between the optimal and suboptimal alignment ending position will >= maskLen. We suggest to
  //                     use query_len/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function
  //                     will NOT return the suboptimal alignment information.
  //                     If the value is not provided, it will default to 15.
  // @return   If the alignment path is accurate (or has missing part). 0: accurate; 1: banded_sw is totally failed; 2: banded_sw returned path has missing part
  // =========
  uint16_t Align(
      const char*   query,
      size_t        query_len,
      const Filter& filter,
      Alignment&    alignment,
      int32_t       maskLen = 0) const;

  // =========
  // @function Like Align(const char* query, size_t query_len, ...), but
  // requires null-terminated query.
  // =========
  uint16_t Align(
      const char* query, const Filter& filter, Alignment& alignment, int32_t maskLen = 0) const;

  // =========
  // @function Align the query against the reference.
  //           [NOTICE] The reference won't replace the reference set by SetReferenceSequence.
  // @param    query     The query sequence.
  //                     [NOTICE] Does not have to be null terminated.
  // @param    query_len The number of bases in query.
  // @param    ref       The reference sequence.
  //                     [NOTICE] Does not have to be null terminated.
  // @param    ref_len   The number of bases in ref.
  // @param    filter    The filter for the alignment.
  // @param    alignment The container contains the result.
  // @param    maskLen   The distance between the optimal and suboptimal alignment ending position will >= maskLen. We suggest to
  //                     use query_len/2, if you don't have special concerns. Note: maskLen has to be >= 15, otherwise this function
  //                     will NOT return the suboptimal alignment information.
  //                     If the value is not provided, it will default to 15.
  // @return   If the alignment path is accurate (or has missing part). 0: accurate; 1: banded_sw is totally failed; 2: banded_sw returned path has missing part
  // =========
  uint16_t Align(
      const char*   query,
      size_t        query_len,
      const char*   ref,
      size_t        ref_len,
      const Filter& filter,
      Alignment&    alignment,
      int32_t       maskLen = 0) const;

  // =========
  // @function Like Align(const char* query, size_t query_len, const char* ref,
  // size_t ref_len, ...), but requires null-terminated query and ref.
  // =========
  uint16_t Align(
      const char*   query,
      const char*   ref,
      const Filter& filter,
      Alignment&    alignment,
      int32_t       maskLen = 0) const;

  // @function Clear up all containers and thus the aligner is disabled.
  //             To rebuild the aligner please use Build functions.
  void Clear();

  // =========
  // @function Rebuild the aligner's ability on default values.
  //           [NOTICE] If the aligner is not cleaned, rebuilding will fail.
  // @return   True: succeed; false: fail.
  // =========
  bool ReBuild();

  // =========
  // @function Rebuild the aligner's ability by the specific matrixs.
  //           [NOTICE] If the aligner is not cleaned, rebuilding will fail.
  // @return   True: succeed; false: fail.
  // =========
  bool ReBuild(
      uint8_t match_score,
      uint8_t mismatch_penalty,
      uint8_t gap_opening_penalty,
      uint8_t gap_extending_penalty);

  // =========
  // @function Construct an Aligner by the specific matrixs.
  //           [NOTICE] If the aligner is not cleaned, rebuilding will fail.
  // @return   True: succeed; false: fail.
  // =========
  bool ReBuild(
      const int8_t* score_matrix,
      int           score_matrix_size,
      const int8_t* translation_matrix,
      int           translation_matrix_size);

private:
  uint16_t AlignImpl(
      const char* const          query,
      const size_t               query_len,
      const std::vector<int8_t>& translated_ref,
      const Filter&              filter,
      Alignment&                 alignment,
      int32_t                    maskLen) const;

  std::vector<int8_t> TranslateBase(const char* str, size_t strLen) const;

  void SetAllDefault();

  void BuildDefaultMatrix();

  void ClearMatrices();

private:
  struct Parameters {
    Parameters(
        uint8_t match_score           = 2,
        uint8_t mismatch_penalty      = 2,
        uint8_t gap_opening_penalty   = 3,
        uint8_t gap_extending_penalty = 1)
        : match_score_(match_score)
        , mismatch_penalty_(mismatch_penalty)
        , gap_opening_penalty_(gap_opening_penalty)
        , gap_extending_penalty_(gap_extending_penalty) {}

    uint8_t match_score_;
    uint8_t mismatch_penalty_;
    uint8_t gap_opening_penalty_;
    uint8_t gap_extending_penalty_;
  } p_;

  int                 score_matrix_size_;
  std::vector<int8_t> score_matrix_;
  std::vector<int8_t> translation_matrix_;
  std::vector<int8_t> translated_reference_;
};  // class Aligner

}  // namespace StripedSmithWaterman

#endif  // COMPLETE_STRIPED_SMITH_WATERMAN_CPP_H_
