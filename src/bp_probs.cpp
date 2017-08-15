#include "bp_probs.h"

#include <string>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

extern "C" {
#include "ViennaRNA/model.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/part_func.h"
}

size_t get_pos_from_vc(vrna_fold_compound_t  *vc, size_t i,size_t j)
{
    return vc->iindx[i]-j;
}

float get_prob_from_vc(vrna_fold_compound_t *vc, size_t i, size_t j)
{
    return vc->exp_matrices->probs[get_pos_from_vc(vc, i,j)];
}

int get_bp_prob(const std::string &seq, vrna_md_t *md, struct bp_prob *bp)
{
    float f;
    unsigned int i, j;
    vrna_fold_compound_t *vc;

    vc  = vrna_fold_compound(seq.c_str(), md, VRNA_OPTION_PF);
    vrna_pf(vc, NULL);
    
    for (i = 1; i <= seq.length(); ++i)
    {
        for (j = 1; j <= seq.length(); ++j)
        {
            f = get_prob_from_vc(vc, i, j);
            //if (f)
            //printf("%d %d %1.9f %1.9f\n", i, j, f, sqrt(f));
            //printf("m[%d][%d] = %1.9f;\n", i, j, f);
            if (f > 0.01)
                bp->m[i][j] = f * 20;
            else
                bp->m[i][j] = 0;
        }
    }
    vrna_fold_compound_free(vc);
    return 0;
}

void set_md_config(vrna_md_t *md)
{
    vrna_md_set_default(md);

    md->backtrack = 0;
    md->noLP = 0;
    md->compute_bpp = 1;
}

int get_bp_prob(const std::string &seq, struct bp_prob *bp)
{
    vrna_md_t md;
    set_md_config(&md);
    get_bp_prob(seq.c_str(), &md, bp);
    return 0;
}

int _main(void)
{
    vrna_md_t md;
    struct bp_prob *bp1, *bp2;
    std::string seq1 = "CGCAGGGAUACCCGCG";
    std::string seq2 = "GCGCCCAUAGGGACGC";

    bp1 = new struct bp_prob();
    bp2 = new struct bp_prob();

    set_md_config(&md);
    get_bp_prob(seq1.c_str(), &md, bp1);
    get_bp_prob(seq2.c_str(), &md, bp2);

    delete bp1;
    delete bp2;
    return 0;
}
