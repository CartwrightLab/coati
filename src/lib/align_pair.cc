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
    size_t len_a = a.length() + 1;           // length of ancestor
    size_t len_b = b.length() + 1;           // length of descendant
    work.mch.resize(len_a, len_b, lowest);   // match matrix
    work.del.resize(len_a, len_b, lowest);   // deletion matrix
    work.ins.resize(len_a, len_b, lowest);   // insertion matrix
    work.path_mch.resize(len_a, len_b, -1);  // path match matrix
    work.path_del.resize(len_a, len_b, -1);  // path deletion matrix
    work.path_ins.resize(len_a, len_b, -1);  // path insertion matrix

    // initialize the margins of the matrices
    work.mch(0, 0) = 0.0;
    work.path_mch(0, 0) = 0;

    for(size_t i = 1; i < len_a; ++i) {
        work.mch(i, 0) = no_gap + gap_open + gap_extend * (i - 1);
        work.del(i, 0) = no_gap + gap_open + gap_extend * (i - 1);
        work.path_mch(i, 0) = 1;  // del
        work.path_del(i, 0) = 2;  // gap extension
    }
    for(size_t j = 2; j < len_b; ++j) {
        work.mch(0, j) = gap_open + gap_extend * (j - 1);
        work.ins(0, j) = gap_open + gap_extend * (j - 1);
        work.path_mch(0, j) = 2;  // ins
        work.path_ins(0, j) = 2;  // gap extension
    }

    work.path_del(1, 0) = 1;
    work.path_ins(0, 1) = 1;

    // fill the body of the matrices
    for(size_t i = 1; i < len_a; ++i) {
        for(size_t j = 1; j < len_b; ++j) {
            // matches can follow matches, ins, or del
            switch(work.path_mch(i - 1, j - 1)) {
            case 0:
                work.mch(i, j) = 2 * no_gap + work.mch(i - 1, j - 1) +
                                 match(a[i - 1], b[j - 1]);
                break;
            case 1:
                work.mch(i, j) = gap_stop + work.mch(i - 1, j - 1) +
                                 match(a[i - 1], b[j - 1]);
                break;
            case 2:
                work.mch(i, j) = gap_stop + no_gap + work.mch(i - 1, j - 1) +
                                 match(a[i - 1], b[j - 1]);
                break;
            }

            // deletions can follow matches, ins, or del
            auto del_del = work.del(i - 1, j) + gap_extend;
            auto del_mch = work.mch(i - 1, j) + no_gap + gap_open;
            auto del_ins = work.ins(i - 1, j) + gap_stop + gap_open;

            // insertions can follow matches or ins.
            auto ins_mch = work.mch(i, j - 1) + gap_open;
            auto ins_ins = work.ins(i, j - 1) + gap_extend;

            // save score
            work.del(i, j) = maximum(del_mch, del_del, del_ins);
            work.ins(i, j) = maximum(ins_mch, ins_ins);

            // save ins & del path: (1) gap open, (2) gap extend
            work.path_del(i, j) =
                ((del_del > del_mch) && (del_del > del_ins)) ? 2 : 1;
            work.path_ins(i, j) = (ins_mch > ins_ins) ? 1 : 2;

            // save mch path: (0) match, (1) deletion, (2) insertion
            int path = 0;
            auto score = work.mch(i, j);
            if(work.del(i, j) > score) {
                score = work.del(i, j);
                path = 1;
            }
            if(work.ins(i, j) > score) {
                score = work.ins(i, j);
                path = 2;
            }
            work.mch(i, j) = score;
            work.path_mch(i, j) = path;
        }
    }
    {
        // adjust the terminal state
        work.mch(len_a - 1, len_b - 1) += no_gap;
        work.ins(len_a - 1, len_b - 1) += gap_stop;
        int path = 0;
        auto score = work.mch(len_a - 1, len_b - 1);
        if(work.del(len_a - 1, len_b - 1) > score) {
            score = work.del(len_a - 1, len_b - 1);
            path = 1;
        }
        if(work.ins(len_a - 1, len_b - 1) > score) {
            score = work.ins(len_a - 1, len_b - 1);
            path = 2;
        }
        work.path_mch(len_a - 1, len_b - 1) = path;
        work.mch(len_a - 1, len_b - 1) = score;
        work.path_mch(0, 0) = -1;
    }
}

void traceback(const align_pair_work_t &work, const std::string &a,
               const std::string &b, alignment_t &aln) {
    // sequence_pair_t ret(2);
    size_t i = work.path_mch.rows() - 1;
    size_t j = work.path_mch.cols() - 1;

    aln.weight = work.mch(i, j);

    aln.f.seq_data.resize(2);
    aln.f.seq_data[0].reserve(i + j);
    aln.f.seq_data[1].reserve(i + j);

    for(int m = work.path_mch(i, j); m != -1; m = work.path_mch(i, j)) {
        switch(m) {
        case 0:  // match
            aln.f.seq_data[0].push_back(a[i - 1]);
            aln.f.seq_data[1].push_back(b[j - 1]);
            i--;
            j--;
            break;
        case 1:  // deletion
            do {
                aln.f.seq_data[0].push_back(a[i - 1]);
                aln.f.seq_data[1].push_back('-');
                i--;
            } while(work.path_del(i + 1, j) == 2);
            break;
        case 2:  // insertion
            do {
                aln.f.seq_data[0].push_back('-');
                aln.f.seq_data[1].push_back(b[j - 1]);
                j--;
            } while(work.path_ins(i, j + 1) == 2);
            break;
        }
    }

    std::reverse(aln.f.seq_data[0].begin(), aln.f.seq_data[0].end());
    std::reverse(aln.f.seq_data[1].begin(), aln.f.seq_data[1].end());
}

}  // namespace coati
