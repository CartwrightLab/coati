#include <doctest/doctest.h>

#include <coati/align_pair.hpp>

namespace coati {

// m -> m      (1-g)*(1-g)*P(b[j] | a[i])/P(b[j])
// m -> d      (1-g)*g
// m -> i      g
// m -> END    1-g
// d -> m      1-e
// d -> d      e
// d -> i      0
// d -> END    1
// i -> m      (1-e)*(1-g)
// i -> d      (1-e)*g
// i -> i      e
// i -> END    (1-e)

// a is "ancestor"
// b is "descendant"
// match matrix holds log(P(b[j]|a[i]))-log(P(b[j]))
// gap_open = log(g)
// gap_extend = log(e)

/**
 * \brief Implementation of Forward algorithm.
 *
 * @param[in,out] work align_pair_work_t dynamic programming matrices.
 * @param[in] a coati::seq_view_t encoded reference/ancestor sequence.
 * @param[in] b coati::seq_view_t encoded descendant sequence.
 * @param[in] aln coati::alignment_t alignment parameters.
 */
template <class S, class W>
void forward_impl(W &work, const seq_view_t &a, const seq_view_t &b,
                  const alignment_t &aln) {
    // calculate log(1-g) log(1-e) log(g) log(e)
    float_t no_gap = S::from_linear_1mf(aln.gap.open);
    float_t gap_stop = S::from_linear_1mf(aln.gap.extend);
    float_t gap_open = S::from_linearf(aln.gap.open);
    float_t gap_extend = S::from_linearf(aln.gap.extend);

    size_t look_back = aln.gap.len;
    size_t start = look_back - 1;

    const float_t lowest = std::numeric_limits<float_t>::lowest();

    // create matrices
    size_t len_a = a.length() + look_back;  // length of ancestor
    size_t len_b = b.length() + look_back;  // length of descendant
    work.resize(len_a, len_b, lowest);      // resize work matrices

    // initialize the margins of the matrices
    work.mch(start, start) = 0.0;

    for(size_t i = start + look_back; i < len_a; i += look_back) {
        work.del(i, start) =
            no_gap + gap_open + gap_extend * static_cast<float_t>(i - 1);
    }
    for(size_t j = start + look_back; j < len_b; j += look_back) {
        work.ins(start, j) =
            gap_open + gap_extend * static_cast<float_t>(j - 1);
    }
    work.init_margins();

    // fill the body of the matrices
    for(size_t i = look_back; i < len_a; ++i) {
        for(size_t j = look_back; j < len_b; ++j) {
            //  from match, ins, or del to match
            work.subst = aln.subst_matrix(a[i - look_back], b[j - look_back]);
            work.mch2mch = work.mch(i - 1, j - 1) + 2 * no_gap + work.subst;
            work.del2mch = work.del(i - 1, j - 1) + gap_stop + work.subst;
            work.ins2mch =
                work.ins(i - 1, j - 1) + gap_stop + no_gap + work.subst;

            // from match, del, or ins to del
            work.mch2del = work.mch(i - look_back, j) + no_gap + gap_open +
                           gap_extend * static_cast<float_t>(look_back - 1);
            work.ins2del = work.ins(i - look_back, j) + gap_stop + gap_open +
                           gap_extend * static_cast<float_t>(look_back - 1);
            work.del2del = work.del(i - look_back, j) +
                           gap_extend * static_cast<float_t>(look_back);

            // from match  or ins to ins
            work.mch2ins = work.mch(i, j - look_back) + gap_open +
                           gap_extend * static_cast<float_t>(look_back - 1);
            work.ins2ins = work.ins(i, j - look_back) +
                           gap_extend * static_cast<float_t>(look_back - 1);

            // save score
            work.mch(i, j) = S::plus(work.mch2mch, work.del2mch, work.ins2mch);
            work.del(i, j) = S::plus(work.mch2del, work.del2del, work.ins2del);
            work.ins(i, j) = S::plus(work.mch2ins, work.ins2ins);

            // save intermediate values
            work.save_values(i, j);
        }
    }
    {
        // adjust the terminal state
        work.mch(len_a - 1, len_b - 1) += no_gap;
        work.ins(len_a - 1, len_b - 1) += gap_stop;
    }
}

/**
 * \brief Forward algorithm.
 *
 * @param[in] a coati::seq_view_t encoded reference/ancestor sequence.
 * @param[in] b coati::seq_view_t encoded descendant sequence.
 * @param[in] aln coati::alignment_t alignment parameters.
 */
void forward(align_pair_work_t &work, const seq_view_t &a, const seq_view_t &b,
             const alignment_t &aln) {
    coati::forward_impl<coati::semiring::log>(work, a, b, aln);
}

/**
 * \brief Memory efficient Forward algorithm.
 *
 * @param[in] a coati::seq_view_t encoded reference/ancestor sequence.
 * @param[in] b coati::seq_view_t encoded descendant sequence.
 * @param[in] aln coati::alignment_t alignment parameters.
 */
void forward_mem(align_pair_work_mem_t &work, const seq_view_t &a,
                 const seq_view_t &b, const alignment_t &aln) {
    coati::forward_impl<coati::semiring::log>(work, a, b, aln);
}

/**
 * \brief Dynamic programming of pairwise alignment.
 *
 * Gotoh-like algorithm for finding the best alignment of two sequences.
 *  Fill matrices using dynamic programming O(n*m).
 *
 * @param[in] a coati::seq_view_t encoded reference/ancestor sequence.
 * @param[in] b coati::seq_view_t encoded descendant sequence.
 * @param[in] aln coati::alignment_t alignment parameters.
 */
void viterbi(align_pair_work_t &work, const seq_view_t &a, const seq_view_t &b,
             const alignment_t &aln) {
    coati::forward_impl<coati::semiring::tropical>(work, a, b, aln);
}

/**
 * \brief Memory efficient dynamic programming of pairwise alignment.
 *
 * Gotoh-like algorithm for finding the best alignment of two sequences.
 *  Fill matrices using dynamic programming O(n*m).
 *
 * @param[in] a coati::seq_view_t encoded reference/ancestor sequence.
 * @param[in] b coati::seq_view_t encoded descendant sequence.
 * @param[in] aln coati::alignment_t alignment parameters.
 */
void viterbi_mem(align_pair_work_mem_t &work, const seq_view_t &a,
                 const seq_view_t &b, const alignment_t &aln) {
    coati::forward_impl<coati::semiring::tropical>(work, a, b, aln);
}

enum struct AlnState { MATCH, DELETION, INSERTION };

/** \brief Index of max value between match, deletion, insertion coati::float_t
 * values in log space.*/
AlnState max_mdi(float_t mch, float_t del, float_t ins) {
    auto i = AlnState::MATCH;
    float_t val = mch;
    if(del > val) {
        val = del;
        i = AlnState::DELETION;
    }
    if(ins > val) {
        return AlnState::INSERTION;
    }
    return i;
}

/** \brief Index of max value between match and insertion coati::float_t values
 * in log space.*/
AlnState max_mi(float_t mch, float_t ins) {
    return mch > ins ? AlnState::MATCH : AlnState::INSERTION;
}

/**
 * \brief Implementation of traceback algorithm.
 *
 * Trace back a set of filled matrices from coati::align_pair to retrieve
 *  the pairwise aligned sequences.
 *
 * @param[in] work coati::align_pair_work_mem_t filled dynamic programming
 * matrices.
 * @param[in] a std::string reference/ancestor sequence.
 * @param[in] b std::string descendant sequence.
 * @param[in,out] aln coati::alignment_t aligned sequences data.
 * @param[in] look_back std::size_t gap unit size.
 *
 */
void traceback(const align_pair_work_mem_t &work, const std::string &a,
               const std::string &b, alignment_t &aln, size_t look_back) {
    size_t i = static_cast<int>(work.mch.rows() - 1);
    size_t j = static_cast<int>(work.mch.cols() - 1);

    coati::semiring::log s_ring;
    float_t no_gap = s_ring.from_linear_1mf(aln.gap.open);
    float_t gap_stop = s_ring.from_linear_1mf(aln.gap.extend);
    float_t gap_open = s_ring.from_linearf(aln.gap.open);
    float_t gap_extend = s_ring.from_linearf(aln.gap.extend);

    aln.data.seqs.clear();
    aln.data.seqs.resize(2);
    aln.seq(0).reserve(i + j);
    aln.seq(1).reserve(i + j);

    aln.data.weight = maximum(work.mch(i, j), work.del(i, j), work.ins(i, j));
    auto m = max_mdi(work.mch(i, j), work.del(i, j), work.ins(i, j));

    while((j > (look_back - 1)) || (i > (look_back - 1))) {
        switch(m) {
        case AlnState::MATCH:  // match
            aln.seq(0).push_back(a[i - look_back]);
            aln.seq(1).push_back(b[j - look_back]);
            i--;
            j--;
            m = max_mdi(work.mch(i, j) + 2 * no_gap, work.del(i, j) + gap_stop,
                        work.ins(i, j) + gap_stop + no_gap);
            break;
        case AlnState::DELETION:  // deletion
            for(size_t k = i; k > (i - look_back); k--) {
                aln.seq(0).push_back(a[k - look_back]);
                aln.seq(1).push_back('-');
            }
            i -= look_back;
            m = max_mdi(work.mch(i, j) + no_gap + gap_open,
                        work.del(i, j) + gap_extend,
                        work.ins(i, j) + gap_stop + gap_open);
            break;
        case AlnState::INSERTION:  // insertion
            for(size_t k = j; k > (j - look_back); k--) {
                aln.seq(0).push_back('-');
                aln.seq(1).push_back(b[k - look_back]);
            }
            j -= look_back;
            m = max_mi(work.mch(i, j) + gap_open, work.ins(i, j) + gap_extend);
            break;
        }
    }

    std::reverse(aln.seq(0).begin(), aln.seq(0).end());
    std::reverse(aln.seq(1).begin(), aln.seq(1).end());
}

/** \brief Sample from match, insertion, and deletion.
 *
 * @param[in] log_mch float log match value.
 * @param[in] log_del float log deletion value.
 * @param[in] log_ins float log insertion value.
 * @param[in] p float random value.
 *
 * \return std::pair of match, deletion, or insertion (AlnState) and value
 * (float).
 */
std::pair<AlnState, float> sample_mdi(float log_mch, float log_del,
                                      float log_ins, float p) {
    float mch = ::expf(log_mch);
    float del = ::expf(log_del);
    float ins = ::expf(log_ins);
    float scale = mch + del + ins;
    p *= scale;
    AlnState ret{AlnState::MATCH};
    float weight{0.f};
    if(p < mch) {
        ret = AlnState::MATCH;
        weight = log_mch;
    } else if(p < del + mch) {
        ret = AlnState::DELETION;
        weight = log_del;
    } else {
        ret = AlnState::INSERTION;
        weight = log_ins;
    }
    weight = weight - ::logf(scale);
    return {ret, weight};
}

/** \brief Sample from match and insertion.
 *
 * @param[in] log_mch float log match value.
 * @param[in] log_del float log deletion value.
 * @param[in] log_ins float log insertion value.
 * @param[in] p float random value.
 *
 * \return std::pair of match or insertion (AlnState) and value (float).
 */

std::pair<AlnState, float> sample_mi(float log_mch, float log_ins, float p) {
    float mch = ::expf(log_mch);
    float ins = ::expf(log_ins);
    float scale = mch + ins;
    p *= scale;
    AlnState ret{AlnState::MATCH};
    float weight{0.f};
    if(p < mch) {
        ret = AlnState::MATCH;
        weight = log_mch;
    } else {
        ret = AlnState::INSERTION;
        weight = log_ins;
    }
    weight = weight - ::logf(scale);
    return {ret, weight};
}

/**
 * \brief Traceback with sampling.
 *
 * Get alignment from dynamic programming matrices with sampling (i.e. not
 * necessarily best alignment).
 *
 * @param[in,out] work coati::align_pair_work_t dynamic programming matrices.
 * @param[in] a std::string ancestor (ref) sequence.
 * @param[in] b std::string descendant sequence.
 * @param[in,out] aln coati::alignment_t alignment object.
 * @param[in] look_back std::size_t size of unit gap.
 * @param[in] rand coati::random_t random number generator object.
 *
 */
void sampleback(const align_pair_work_t &work, const std::string &a,
                const std::string &b, alignment_t &aln, size_t look_back,
                random_t &rand) {
    size_t i = work.mch_mch.rows() - 1;
    size_t j = work.mch_mch.cols() - 1;

    aln.data.seqs.resize(2);
    aln.seq(0).clear();
    aln.seq(1).clear();
    aln.seq(0).reserve(i + j);
    aln.seq(1).reserve(i + j);
    aln.data.weight = 0.0f;

    float w = maximum(work.mch(i, j), work.del(i, j), work.ins(i, j));
    auto pick = sample_mdi(work.mch(i, j) - w, work.del(i, j) - w,
                           work.ins(i, j) - w, rand.f24());
    aln.data.weight += pick.second;

    while((j > (look_back - 1)) || (i > (look_back - 1))) {
        switch(pick.first) {
        case AlnState::MATCH:
            aln.seq(0).push_back(a[i - look_back]);
            aln.seq(1).push_back(b[j - look_back]);
            w = work.mch(i, j);
            pick = sample_mdi(work.mch_mch(i, j) - w, work.del_mch(i, j) - w,
                              work.ins_mch(i, j) - w, rand.f24());
            aln.data.weight += pick.second;
            i--;
            j--;
            break;
        case AlnState::DELETION:
            for(size_t k = i; k > (i - look_back); k--) {
                aln.seq(0).push_back(a[k - look_back]);
                aln.seq(1).push_back('-');
            }
            w = work.del(i, j);
            pick = sample_mdi(work.mch_del(i, j) - w, work.del_del(i, j) - w,
                              work.ins_del(i, j) - w, rand.f24());
            aln.data.weight += pick.second;
            i -= look_back;
            break;
        case AlnState::INSERTION:
            for(size_t k = j; k > (j - look_back); k--) {
                aln.seq(0).push_back('-');
                aln.seq(1).push_back(b[k - look_back]);
            }
            w = work.ins(i, j);
            pick = sample_mi(work.mch_ins(i, j) - w, work.ins_ins(i, j) - w,
                             rand.f24());
            aln.data.weight += pick.second;
            j -= look_back;
            break;
        }
    }

    std::reverse(aln.seq(0).begin(), aln.seq(0).end());
    std::reverse(aln.seq(1).begin(), aln.seq(1).end());
}
}  // namespace coati
