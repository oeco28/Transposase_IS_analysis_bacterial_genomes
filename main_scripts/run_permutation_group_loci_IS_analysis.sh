mkdir -p grid
for GAP in 2000 5000; do for F in 0.5 0.8; do for M in 3 10 20; do
  python3 group_by_loci.py \
    --cds inputs/cds.tsv --is_table inputs/is.tsv \
    --genes inputs/genes.tsv --groups inputs/groups_norm.tsv \
    --focal 12 --locus_gaps $GAP --min_frac $F --min_pair_genomes $M \
    --nperm 5000 --loo_nperm 500 --out grid/g${GAP}_f${F}_m${M}
done; done; done

python3 - << 'EOF'
import glob, os, pandas as pd
rows = []
for d in sorted(glob.glob("grid/*")):
    t = pd.read_csv(f"{d}/group12_locus_tests.tsv", sep="\t").iloc[0]
    L = pd.read_csv(f"{d}/group12_loci.tsv", sep="\t")
    rows.append(dict(setting=os.path.basename(d), n_loci=int(t.n_loci),
                     max_locus_genes=int(L.n_genes.max()),
                     locus_diff=round(t.locus_level_diff, 4), locus_p=t.locus_level_p))
print(pd.DataFrame(rows).to_string(index=False))
EOF
