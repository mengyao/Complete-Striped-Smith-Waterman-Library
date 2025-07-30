// ssw_cpp.cpp
// Created by Wan-Ping Lee
// Last revision by Mengyao Zhao on 2023-Apr-21

#include "ssw_cpp.h"
#include "ssw.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <sstream>

namespace StripedSmithWaterman {
namespace {

const int default_score_matrix_size = 5;

const int8_t kBaseTranslation[128] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    // A     C           G                                      T
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    // a     c           g                                      t
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

void BuildSwScoreMatrix(
    const uint8_t match_score, const uint8_t mismatch_penalty, std::vector<int8_t>& matrix) {
  // The score matrix looks like
  //                 // A,  C,  G,  T,  N
  //  score_matrix_ = { 2, -2, -2, -2, -2, // A
  //                   -2,  2, -2, -2, -2, // C
  //                   -2, -2,  2, -2, -2, // G
  //                   -2, -2, -2,  2, -2, // T
  //                   -2, -2, -2, -2, -2};// N

  int id = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      matrix[id] = ((i == j) ? match_score : -mismatch_penalty);
      ++id;
    }
    matrix[id] = -mismatch_penalty;  // For N
    ++id;
  }

  for (int i = 0; i < 5; ++i) {
    matrix[id] = -mismatch_penalty;  // For N
    ++id;
  }
}

Alignment ConvertAlignment(const s_align& s_al, const int query_len) {
  Alignment result;
  result.sw_score           = s_al.score1;
  result.sw_score_next_best = s_al.score2;
  result.ref_begin          = s_al.ref_begin1;
  result.ref_end            = s_al.ref_end1;
  result.query_begin        = s_al.read_begin1;
  result.query_end          = s_al.read_end1;
  result.ref_end_next_best  = s_al.ref_end2;

  if (s_al.cigarLen > 0) {
    std::ostringstream cigar_string;
    if (result.query_begin > 0) {
      uint32_t cigar = to_cigar_int(result.query_begin, 'S');
      result.cigar.push_back(cigar);
      cigar_string << result.query_begin << 'S';
    }

    for (int i = 0; i < s_al.cigarLen; ++i) {
      const uint32_t cigar_op = s_al.cigar[i];
      result.cigar.push_back(cigar_op);
      cigar_string << cigar_int_to_len(cigar_op) << cigar_int_to_op(cigar_op);
    }

    int end = query_len - result.query_end - 1;
    if (end > 0) {
      uint32_t cigar = to_cigar_int(end, 'S');
      result.cigar.push_back(cigar);
      cigar_string << end << 'S';
    }

    result.cigar_string = cigar_string.str();
  }  // end if

  return result;
}

// @Function:
//     Calculate the length of the previous cigar operator
//     and store it in new_cigar and new_cigar_string.
//     Clean up in_M (false), in_X (false), length_M (0), and length_X(0).
void CleanPreviousMOperator(
    bool&                  in_M,
    bool&                  in_X,
    uint32_t&              length_M,
    uint32_t&              length_X,
    std::vector<uint32_t>& new_cigar,
    std::ostringstream&    new_cigar_string) {
  if (in_M) {
    uint32_t match = to_cigar_int(length_M, '=');
    new_cigar.push_back(match);
    new_cigar_string << length_M << '=';
  } else if (in_X) {
    uint32_t match = to_cigar_int(length_X, 'X');
    new_cigar.push_back(match);
    new_cigar_string << length_X << 'X';
  }

  // Clean up
  in_M     = false;
  in_X     = false;
  length_M = 0;
  length_X = 0;
}

// @Function:
//     1. Calculate the number of mismatches.
//     2. Modify the cigar string:
//         differentiate matches (M) and mismatches(X).
// @Return:
//     The number of mismatches.
int CalculateNumberMismatch(
    Alignment& al, const int8_t* ref, const int8_t* query, const int query_len) {
  ref += al.ref_begin;
  query += al.query_begin;
  int mismatch_length = 0;

  std::vector<uint32_t> new_cigar;
  std::ostringstream    new_cigar_string;

  if (al.query_begin > 0) {
    uint32_t cigar = to_cigar_int(al.query_begin, 'S');
    new_cigar.push_back(cigar);
    new_cigar_string << al.query_begin << 'S';
  }

  bool     in_M     = false;  // the previous is match
  bool     in_X     = false;  // the previous is mismatch
  uint32_t length_M = 0;
  uint32_t length_X = 0;

  for (std::vector<uint32_t>::const_iterator cigar_op_it     = al.cigar.begin(),
                                             cigar_op_it_end = al.cigar.end();
       cigar_op_it != cigar_op_it_end; ++cigar_op_it) {
    const uint32_t cigar_op = *cigar_op_it;
    char           op       = cigar_int_to_op(cigar_op);
    uint32_t       length   = cigar_int_to_len(cigar_op);
    if (op == 'M') {
      for (uint32_t j = 0; j < length; ++j) {
        if (*ref != *query) {
          ++mismatch_length;
          if (in_M) {  // the previous is match; however the current one is mismatched
            uint32_t match = to_cigar_int(length_M, '=');
            new_cigar.push_back(match);
            new_cigar_string << length_M << '=';
          }
          length_M = 0;
          ++length_X;
          in_M = false;
          in_X = true;
        } else {       // *ref == *query
          if (in_X) {  // the previous is mismatch; however the current one is matched
            uint32_t match = to_cigar_int(length_X, 'X');
            new_cigar.push_back(match);
            new_cigar_string << length_X << 'X';
          }
          ++length_M;
          length_X = 0;
          in_M     = true;
          in_X     = false;
        }  // end of if (*ref != *query)
        ++ref;
        ++query;
      }
    } else if (op == 'I') {
      query += length;
      mismatch_length += length;
      CleanPreviousMOperator(in_M, in_X, length_M, length_X, new_cigar, new_cigar_string);
      new_cigar.push_back(cigar_op);
      new_cigar_string << length << 'I';
    } else if (op == 'D') {
      ref += length;
      mismatch_length += length;
      CleanPreviousMOperator(in_M, in_X, length_M, length_X, new_cigar, new_cigar_string);
      new_cigar.push_back(cigar_op);
      new_cigar_string << length << 'D';
    }
  }

  CleanPreviousMOperator(in_M, in_X, length_M, length_X, new_cigar, new_cigar_string);

  int end = query_len - al.query_end - 1;
  if (end > 0) {
    uint32_t cigar = to_cigar_int(end, 'S');
    new_cigar.push_back(cigar);
    new_cigar_string << end << 'S';
  }

  al.cigar_string = new_cigar_string.str();
  al.cigar.swap(new_cigar);

  return mismatch_length;
}

void SetFlag(const Filter& filter, uint8_t& flag) {
  if (filter.report_begin_position) {
    flag |= 0x08;
  }
  if (filter.report_cigar) {
    flag |= 0x0f;
  }
}

}  // namespace

Aligner::Aligner()
    : score_matrix_size_(default_score_matrix_size) {
  BuildDefaultMatrix();
}

Aligner::Aligner(
    const uint8_t match_score,
    const uint8_t mismatch_penalty,
    const uint8_t gap_opening_penalty,
    const uint8_t gap_extending_penalty)
    : p_(match_score, mismatch_penalty, gap_opening_penalty, gap_extending_penalty)
    , score_matrix_size_(default_score_matrix_size) {
  BuildDefaultMatrix();
}

Aligner::Aligner(
    const int8_t* score_matrix,
    const int     score_matrix_size,
    const int8_t* translation_matrix,
    const int     translation_matrix_size)
    : score_matrix_size_(score_matrix_size)
    , score_matrix_(score_matrix, score_matrix + score_matrix_size_ * score_matrix_size_)
    , translation_matrix_(translation_matrix, translation_matrix + translation_matrix_size) {}

size_t Aligner::SetReferenceSequence(const char* const ref, const size_t ref_len) {
  ClearReferenceSequence();
  if (!translation_matrix_.empty()) {
    translated_reference_ = TranslateBase(ref, ref_len);
  }

  return translated_reference_.size();
}

size_t Aligner::SetReferenceSequence(const char* const ref) {
  return SetReferenceSequence(ref, std::strlen(ref));
}

void Aligner::ClearReferenceSequence() { translated_reference_.clear(); }

void Aligner::SetGapPenalty(const uint8_t opening, const uint8_t extending) {
  p_.gap_opening_penalty_   = opening;
  p_.gap_extending_penalty_ = extending;
}

std::vector<int8_t> Aligner::TranslateBase(const char* const str, const size_t str_len) const {
  assert(!translation_matrix_.empty());
  std::vector<int8_t> result;
  result.reserve(str_len);
  for (size_t i = 0; i < str_len; ++i) {
    result.push_back(translation_matrix_[static_cast<unsigned char>(str[i])]);
  }
  return result;
}

uint16_t Aligner::Align(
    const char* const query,
    const size_t      query_len,
    const Filter&     filter,
    Alignment&        alignment,
    int32_t           maskLen) const {
  if (translated_reference_.empty() || (query_len == 0)) {
    return false;
  }

  return AlignImpl(query, query_len, translated_reference_, filter, alignment, maskLen);
}

uint16_t Aligner::Align(
    const char* const query,
    const Filter&     filter,
    Alignment&        alignment,
    const int32_t     maskLen) const {
  return Align(query, std::strlen(query), filter, alignment, maskLen);
}

uint16_t Aligner::Align(
    const char* const query,
    const size_t      query_len,
    const char* const ref,
    const size_t      ref_len,
    const Filter&     filter,
    Alignment&        alignment,
    int32_t           maskLen) const {
  if (translation_matrix_.empty() || (ref_len == 0) || (query_len == 0)) {
    return false;
  }

  const std::vector<int8_t> translated_ref(TranslateBase(ref, ref_len));
  assert(translated_ref.size() == ref_len);

  return AlignImpl(query, query_len, translated_ref, filter, alignment, maskLen);
}

uint16_t Aligner::Align(
    const char* const query,
    const char* const ref,
    const Filter&     filter,
    Alignment&        alignment,
    const int32_t     maskLen) const {
  return Align(query, std::strlen(query), ref, std::strlen(ref), filter, alignment, maskLen);
}

uint16_t Aligner::AlignImpl(
    const char* const          query,
    const size_t               query_len,
    const std::vector<int8_t>& translated_ref,
    const Filter&              filter,
    Alignment&                 alignment,
    int32_t                    maskLen) const {
  if (translation_matrix_.empty()) {
    return false;
  }

  maskLen = std::max(maskLen, 15);

  const std::vector<int8_t> translated_query(TranslateBase(query, query_len));
  assert(translated_query.size() == query_len);

  const int8_t score_size = 2;
  s_profile*   profile    = ssw_init(
      translated_query.data(), translated_query.size(), score_matrix_.data(), score_matrix_size_,
      score_size);

  uint8_t flag = 0;
  SetFlag(filter, flag);
  s_align* s_al = ssw_align(
      profile, translated_ref.data(), translated_ref.size(),
      static_cast<int>(p_.gap_opening_penalty_), static_cast<int>(p_.gap_extending_penalty_), flag,
      filter.score_filter, filter.distance_filter, maskLen);

  alignment            = ConvertAlignment(*s_al, translated_query.size());
  alignment.mismatches = CalculateNumberMismatch(
      alignment, translated_ref.data(), translated_query.data(), translated_query.size());
  uint16_t align_flag = s_al->flag;

  // Free memory
  align_destroy(s_al);
  init_destroy(profile);

  return align_flag;
}

void Aligner::Clear() {
  ClearMatrices();
  ClearReferenceSequence();
}

void Aligner::SetAllDefault() {
  p_                 = Parameters();
  score_matrix_size_ = default_score_matrix_size;
  translated_reference_.clear();
}

bool Aligner::ReBuild() {
  if (!translation_matrix_.empty()) {
    return false;
  }

  SetAllDefault();
  BuildDefaultMatrix();

  return true;
}

bool Aligner::ReBuild(
    const uint8_t match_score,
    const uint8_t mismatch_penalty,
    const uint8_t gap_opening_penalty,
    const uint8_t gap_extending_penalty) {
  if (!translation_matrix_.empty()) {
    return false;
  }

  SetAllDefault();

  p_ = Parameters(match_score, mismatch_penalty, gap_opening_penalty, gap_extending_penalty);

  BuildDefaultMatrix();

  return true;
}

bool Aligner::ReBuild(
    const int8_t* score_matrix,
    const int     score_matrix_size,
    const int8_t* translation_matrix,
    const int     translation_matrix_size) {
  score_matrix_size_ = score_matrix_size;
  score_matrix_.assign(score_matrix, score_matrix + score_matrix_size_ * score_matrix_size_);
  translation_matrix_.assign(translation_matrix, translation_matrix + translation_matrix_size);

  return true;
}

void Aligner::BuildDefaultMatrix() {
  score_matrix_.resize(score_matrix_size_ * score_matrix_size_);
  BuildSwScoreMatrix(p_.match_score_, p_.mismatch_penalty_, score_matrix_);
  translation_matrix_.assign(kBaseTranslation, kBaseTranslation + sizeof(kBaseTranslation));
}

void Aligner::ClearMatrices() {
  score_matrix_.clear();
  translation_matrix_.clear();
}

}  // namespace StripedSmithWaterman
