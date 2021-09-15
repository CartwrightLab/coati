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
void align_pair(align_pair_work_t &work, const seq_view_t &a,
                const seq_view_t &b, const Matrixf &match, float_t gap_open,
                float_t gap_extend, std::vector<float_t> pi) {
    // calculate log(1-g) log(1-e) log(g) log(e) log(pi)
    float_t no_gap = std::log1pf(-gap_open);
    float_t gap_stop = std::log1pf(-gap_extend);
    gap_open = ::logf(gap_open);
    gap_extend = ::logf(gap_extend);
    for(size_t i = 0; i < pi.size(); i++) {
        pi[i] = ::logf(pi[i]);
    }

    const float_t lowest = std::numeric_limits<float_t>::lowest();

    // create matrices
    size_t len_a = a.length() + 1;              // length of ancestor
    size_t len_b = b.length() + 1;              // length of descendant
    work.mch.resize(len_a, len_b, lowest);      // match matrix
    work.del.resize(len_a, len_b, lowest);      // deletion matrix
    work.ins.resize(len_a, len_b, lowest);      // insertion matrix
    work.mch_mch.resize(len_a, len_b, lowest);  // match to match matrix
    work.mch_del.resize(len_a, len_b, lowest);  // match to del matrix
    work.mch_ins.resize(len_a, len_b, lowest);  // match to ins matrix
    work.del_mch.resize(len_a, len_b, lowest);  // del to match matrix
    work.del_del.resize(len_a, len_b, lowest);  // del to del matrix
    work.ins_mch.resize(len_a, len_b, lowest);  // ins to match matrix
    work.ins_ins.resize(len_a, len_b, lowest);  // ins to ins matrix
    work.ins_del.resize(len_a, len_b, lowest);  // ins to del matrix

    // initialize the margins of the matrices
    work.mch(0, 0) = 0.0;

    for(size_t i = 1; i < len_a; ++i) {
        work.del(i, 0) = work.del_del(i, 0) =
            no_gap + gap_open + gap_extend * (i - 1);
    }
    for(size_t j = 1; j < len_b; ++j) {
        work.ins(0, j) = work.ins_ins(0, j) = gap_open + gap_extend * (j - 1);
    }

    // fill the body of the matrices
    for(size_t i = 1; i < len_a; ++i) {
        for(size_t j = 1; j < len_b; ++j) {
            //  from match, ins, or del to match
            auto mch = match(a[i - 1], b[j - 1]);
            work.mch_mch(i, j) = work.mch(i - 1, j - 1) + 2 * no_gap + mch;
            work.del_mch(i, j) = work.del(i - 1, j - 1) + gap_stop + mch;
            work.ins_mch(i, j) =
                work.ins(i - 1, j - 1) + gap_stop + no_gap + mch;

            // from match or del to del
            work.mch_del(i, j) = work.mch(i - 1, j) + no_gap + gap_open;
            work.del_del(i, j) = work.del(i - 1, j) + gap_extend;
            work.ins_del(i, j) = work.ins(i - 1, j) + gap_stop + gap_open;

            // from match, del, or ins to ins
            work.mch_ins(i, j) = work.mch(i, j - 1) + gap_open;
            work.ins_ins(i, j) = work.ins(i, j - 1) + gap_extend;

            // save score
            work.mch(i, j) = maximum(work.mch_mch(i, j), work.del_mch(i, j),
                                     work.ins_mch(i, j));
            work.del(i, j) = maximum(work.mch_del(i, j), work.del_del(i, j),
                                     work.ins_del(i, j));
            work.ins(i, j) = maximum(work.mch_ins(i, j), work.ins_ins(i, j));
        }
    }
    {
        // adjust the terminal state
        work.mch(len_a - 1, len_b - 1) += no_gap;
        work.ins(len_a - 1, len_b - 1) += gap_stop;
    }
}

void traceback(const align_pair_work_t &work, const std::string &a,
               const std::string &b, alignment_t &aln) {
    // sequence_pair_t ret(2);
    size_t i = work.mch_mch.rows() - 1;
    size_t j = work.mch_mch.cols() - 1;

    aln.f.seq_data.resize(2);
    aln.f.seq_data[0].reserve(i + j);
    aln.f.seq_data[1].reserve(i + j);

    aln.weight = work.mch(i, j);
    int m = 0;
    float_t score;
    if(work.del(i, j) > aln.weight) {
        aln.weight = work.del(i, j);
        m = 1;
    }
    if(work.ins(i, j) > aln.weight) {
        aln.weight = work.ins(i, j);
        m = 2;
    }

    while((j > 0) || (i > 0)) {
        switch(m) {
        case 0:  // match
            aln.f.seq_data[0].push_back(a[i - 1]);
            aln.f.seq_data[1].push_back(b[j - 1]);
            score = work.mch_mch(i, j);
            m = 0;
            if(work.del_mch(i, j) > score) {
                score = work.del_mch(i, j);
                m = 1;
            }
            if(work.ins_mch(i, j) > score) {
                m = 2;
            }
            i--;
            j--;
            break;
        case 1:  // deletion
            aln.f.seq_data[0].push_back(a[i - 1]);
            aln.f.seq_data[1].push_back('-');
            score = work.mch_del(i, j);
            m = 0;
            if(work.del_del(i, j) > score) {
                score = work.del_del(i, j);
                m = 1;
            }
            if(work.ins_del(i, j) > score) {
                m = 2;
            }
            i--;
            break;
        case 2:  // insertion
            aln.f.seq_data[0].push_back('-');
            aln.f.seq_data[1].push_back(b[j - 1]);
            score = work.mch_ins(i, j);
            m = 0;
            if(work.ins_ins(i, j) > score) {
                m = 2;
            }
            j--;
            break;
        }
    }

    std::reverse(aln.f.seq_data[0].begin(), aln.f.seq_data[0].end());
    std::reverse(aln.f.seq_data[1].begin(), aln.f.seq_data[1].end());
}

}  // namespace coati
