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
                const seq_view_t &b, const Matrixf &match, input_t &in_data) {
    // calculate log(1-g) log(1-e) log(g) log(e)
    float_t no_gap = std::log1pf(-in_data.gapo);
    float_t gap_stop = std::log1pf(-in_data.gape);
    float_t gap_open = ::logf(in_data.gapo);
    float_t gap_extend = ::logf(in_data.gape);
    size_t look_back = in_data.g_len;

    const float_t lowest = std::numeric_limits<float_t>::lowest();

    // create matrices
    size_t len_a = a.length() + 1;      // length of ancestor
    size_t len_b = b.length() + 1;      // length of descendant
    work.resize(len_a, len_b, lowest);  // resize work matrices

    // initialize the margins of the matrices
    work.mch(0, 0) = 0.0;

    for(size_t i = look_back; i < len_a; i += look_back) {
        work.del(i, 0) = work.del_del(i, 0) =
            no_gap + gap_open + gap_extend * static_cast<float_t>(i - 1);
    }
    for(size_t j = look_back; j < len_b; j += look_back) {
        work.ins(0, j) = work.ins_ins(0, j) =
            gap_open + gap_extend * static_cast<float_t>(j - 1);
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

            if(i > look_back - 1) {
                // from match or del to del
                work.mch_del(i, j) =
                    work.mch(i - look_back, j) + no_gap + gap_open +
                    gap_extend * static_cast<float_t>(look_back - 1);
                work.ins_del(i, j) =
                    work.ins(i - look_back, j) + gap_stop + gap_open +
                    gap_extend * static_cast<float_t>(look_back - 1);
                work.del_del(i, j) =
                    work.del(i - look_back, j) + gap_extend * look_back;
            }
            if(j > look_back - 1) {
                // from match, del, or ins to ins
                work.mch_ins(i, j) =
                    work.mch(i, j - look_back) + gap_open +
                    gap_extend * static_cast<float_t>(look_back - 1);
                work.ins_ins(i, j) =
                    work.ins(i, j - look_back) +
                    gap_extend * static_cast<float_t>(look_back - 1);
            }

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
               const std::string &b, alignment_t &aln, size_t look_back) {
    size_t i = static_cast<int>(work.mch_mch.rows() - 1);
    size_t j = static_cast<int>(work.mch_mch.cols() - 1);

    aln.f.seq_data.resize(2);
    aln.f.seq_data[0].reserve(i + j);
    aln.f.seq_data[1].reserve(i + j);

    aln.weight = maximum(work.mch(i, j), work.del(i, j), work.ins(i, j));
    int m = max_index(work.mch(i, j), work.del(i, j), work.ins(i, j));

    while((j > 0) || (i > 0)) {
        switch(m) {
        case 0:  // match
            aln.f.seq_data[0].push_back(a[i - 1]);
            aln.f.seq_data[1].push_back(b[j - 1]);
            m = max_index(work.mch_mch(i, j), work.del_mch(i, j),
                          work.ins_mch(i, j));
            i--;
            j--;
            break;
        case 1:  // deletion
            for(size_t k = i; k > (i - look_back); k--) {
                aln.f.seq_data[0].push_back(a[k - 1]);
                aln.f.seq_data[1].push_back('-');
            }
            m = max_index(work.mch_del(i, j), work.del_del(i, j),
                          work.ins_del(i, j));
            i -= look_back;
            break;
        case 2:  // insertion
            for(size_t k = j; k > (j - look_back); k--) {
                aln.f.seq_data[0].push_back('-');
                aln.f.seq_data[1].push_back(b[k - 1]);
            }
            m = work.mch_ins(i, j) > work.ins_ins(i, j) ? 0 : 2;
            j -= look_back;
            break;
        }
    }

    std::reverse(aln.f.seq_data[0].begin(), aln.f.seq_data[0].end());
    std::reverse(aln.f.seq_data[1].begin(), aln.f.seq_data[1].end());
}

}  // namespace coati
