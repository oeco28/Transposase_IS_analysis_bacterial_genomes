# Identification of Insertion sequences (IS) from transposase activity in bacterial genomes

This repository contains the main scripts necessary to identify (and catalogue) Insertion sequences (ISs) in all bacterial genomes using ISEScan. We further limit the count of IS elements to intact elements, which we refer to as "complete elements" or "complete' ISs. 

We investigate the evidence for recent activity of transposase activity and further test if this activity is more pronounced in a focal group (gruop 12 in our case, but this can be modified while running the pipeline using the --focal-group flag).

The scripts included in main are required for producing the results presented in our manuscript, and the additional script folder allows additional analyses, not included in our work, of the IS distribution in the genomes of *Streptococcus pyogenes*, including the analysis of the recombination map.

## Installation of environment

We provide with an installation script that creates a conda environment *ismobile*. A configuration file (env_core.yml) needs to be downloaded and to create the environment and install the necessary dependencies, the user should run:

```bash
./install_tools.sh
```

(the env_core.yml file needs to be in the same directory as the install_tools.sh script)

We also provide with an additional installation option that includes optional tools not included in the analysis of the manuscript but available to any researcher interested in performing additional analyses. One of the potential analyses is the association between the recombination map recovered from gubbins and the distribution of IS sequences in the genomes of focal groups of interest.

In order to install the additional dependencies, the file env_optional.yml should be saved in the same directory of the install_tools.sh and run

```bash
./install_tools.sh --optional
```
After the tool runs, the installation can be checked with

```bash
./install_tools.sh --verify-only
```
NOTE: the installation step also downloads the necessary databases for the identification of insertion sequences (IS). If for some reason, the process of download of databases fail, you can run the following command to download databases and prepare them for analyses:

```bash
./install_tools.sh --dbs
```
or
```bash
./install_tools --optional --dbs
```

## IS distribution analyses

In order to run analyses, first activate the conda environment

```bash
conda activate ismobile
```

### Step 1: prepare inputs

```bash
bash prep_inputs.sh genomes fastbaps_clusters.csv
```

Creates `assemblies/` (symlinks named `<genome>.fna`, since script build_panel.py globs `*.fna`) and `metadata.tsv` with genome, fastbaps_group, n_contigs, N50 and total_len. Gzipped inputs are decompressed rather than linked, because ISEScan will not read them.

Read the group-size summary it prints and note the label of your focal group. Everything below assumes `12`; substitute yours if interested in running our pipeline.

### Step 2: build the genome panels

```bash
python build_panel.py --meta metadata.tsv --tree core_gene_tree.nwk --focal 12
```

Writes `panel_discovery.txt` (all genomes), `panel_scoring.txt`, and `panel_weights.tsv`. Check the clonal-redundancy warning before going further. If clonal cluster sizes are uneven across groups, an apparent IS excess in the focal group can be produced entirely by that.

## Step 3: annotate IS elements

```bash
bash annotate_IS.sh assemblies is_out 16     # dir, outdir, threads
```

This script is resumable. Rerun it after any timeout or failure and it skips completed genomes, completed stages, and completed BLAST chunks. `--status` shows what is left without running anything, `--force` redoes everything, `--stage <annotate|collate|extract|blast>` runs one stage.

For ~1900 genomes, you can run this step as a SLURM array instead of one long job:

```bash
bash annotate_IS.sh assemblies is_out 1 --status   # writes todo.txt + array_job.sh
sbatch --array=1-$(wc -l < is_out/todo.txt)%50 array_job.sh
bash annotate_IS.sh assemblies is_out 16           # collate/extract/blast after
```
the tool is easy to port to a system managed with SLURM

This step runs ISEScan on every genome, collates `is_out/IS_calls.tsv`, extracts `is_out/IS_seqs.fna`, and does the all-vs-all identity search.

Run this on ALL genomes, not just the focal group. Restricting it creates ascertainment bias you cannot undo later. Budget a few CPU minutes per genome; it parallelizes across assemblies.

Sanity check before moving on: the family counts printed at the end should be dominated by a handful of IS families, and the complete (TIR-intact) fraction should be somewhere in the tens of percent. If nearly everything is partial, the assemblies are too fragmented for the insertion-site analysis.

### Step 4: IS load and recency

```bash
Rscript 02_IS_load_and_recency.R --tree strep_tree.nwk --meta metadata.tsv --focal 12 2>&1 | tee is_out/02_log_new.txt
```

Edit the `FOCAL` variable at the top first. Reads `is_out/IS_calls.tsv`, `metadata.tsv`, `core_gene_tree.nwk`. Produces the quasibinomial models, a phylogenetic ANOVA (if you are interested in checking how phylogenetic inertia plays a role in what you are detecting), and two exploratory PDFs for your own consumption. It allows for visualization of the number of IS identified as a function of the number of genomes analyzed.

Read the models in order. If the focal effect disappears once `n_contigs` is in the model, that is the answer, and it is worth knowing before you build anything else on top of it.

### Step 4a: insertion-site occupancy matrix

```bash
python insertion_site_polymorphism.py \
    --calls is_out/IS_calls.tsv \
    --asm assemblies \
    --meta metadata.tsv \
    --threads 16
```

This is the slowest step and the one that matters most. It BLASTs every probe against every genome, so runtime scales as roughly (IS copies) x (genomes). For a large collection, restrict scoring to `panel_scoring.txt` by passing a filtered metadata file.

Outputs `is_out/insertion_matrix.tsv` (occupied/empty/unknown per genome per site) and `is_out/insertion_frequencies.tsv`.

Set the focal group label inside the script (it is hardcoded as `12` near the bottom) before running.

### Step 5: selection and recombination

```bash
python IS_selection_and_recombination.py \
    --focal 12 \
    --calls is_out/IS_calls.tsv \
    --gff /path_to_converted_gffs \
    --pa /path_to_gene_presence_absence.csv \
    --selected /path_to_list_selected_genes.txt \
    --meta metadata.tsv
```

As written, `__main__` runs only test A (physical association between selected accessory genes and IS elements). Tests B and C are importable functions, because they need inputs specific to your setup:

- `test_B(is_ref_coords, gubbins_path, ref_len)` needs IS positions
  projected onto whichever reference genome gubbins was run against.
  Get those by BLASTing `is_out/IS_seqs.fna` against that reference.

- `test_C(matrix_path, tree_path)` needs step 03 to have finished.

Although the recombination analysis is implemented, we mainly explored the impact of transposase activity on the movement of genes under selection in group A Strep.

