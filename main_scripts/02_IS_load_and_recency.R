## ------------------------------------------------------------------
## 02_IS_load_and_recency.R
##
## Q1: does the focal fastbaps group really carry more IS elements,
##     after controlling for assembly fragmentation and phylogeny?
## Q2: is the excess made of recent copies (active transposition) or
##     old copies (an ancient expansion inherited vertically)?
## ------------------------------------------------------------------

library(data.table); library(MASS); library(ape); library(phytools)
library(ggplot2)

is_calls <- fread("is_out/IS_calls.tsv")
meta     <- fread("metadata.tsv")          # genome, fastbaps_group, n_contigs, N50, total_len
tree     <- read.tree("core_gene_tree.nwk")
FOCAL    <- "G4"                            # <- your interspersed group

## ---------- Q1: IS load ------------------------------------------
load_tab <- is_calls[, .(n_IS       = .N,
                         n_complete = sum(complete),
                         n_partial  = sum(!complete)), by = genome]
d <- merge(meta, load_tab, by = "genome", all.x = TRUE)
d[is.na(n_IS), c("n_IS","n_complete","n_partial") := 0]
d[, focal := factor(fifelse(fastbaps_group == FOCAL, "focal", "other"),
                    levels = c("other","focal"))]

## Short-read assemblies collapse repeated IS copies and break contigs
## at them, so n_contigs is both a confounder and a weak proxy for IS
## load. Put it in the model and report the effect with and without it;
## if the focal effect survives, that is a much harder result.
m0 <- glm.nb(n_IS ~ focal + offset(log(total_len)), data = d)
m1 <- glm.nb(n_IS ~ focal + scale(log(n_contigs)) + scale(log(N50)) +
                    offset(log(total_len)), data = d)
print(summary(m1)); print(anova(m0, m1))

## Complete (TIR-intact) copies only: this is the fraction that can
## still move, and it is the number your argument actually needs.
m2 <- glm.nb(n_complete ~ focal + scale(log(n_contigs)) +
                          offset(log(total_len)), data = d)
print(summary(m2))

## Group membership is confounded with phylogeny: the focal group is a
## clade (mostly), so ordinary GLM p-values are anticonservative.
## Phylogenetic ANOVA with tip-label simulation on the core tree.
x <- setNames(d$focal, d$genome)[tree$tip.label]
y <- setNames(log1p(d$n_IS / d$total_len * 1e6), d$genome)[tree$tip.label]
print(phylANOVA(tree, x, y, nsim = 10000))

## Family-level: which IS families drive the excess? A single expanding
## family is a very different story from a general permissiveness to MGE.
fam <- is_calls[, .N, by = .(genome, family)]
fam <- merge(fam, d[, .(genome, focal, total_len, n_contigs)], by = "genome")
for (f in unique(fam$family)) {
  sub <- fam[family == f]
  if (nrow(sub) < 20) next
  mf <- try(glm.nb(N ~ focal + scale(log(n_contigs)) + offset(log(total_len)),
                   data = sub), silent = TRUE)
  if (!inherits(mf, "try-error"))
    cat(f, coef(summary(mf))["focalfocal", c(1,4)], "\n")
}

## ---------- Q2: how recent are the copies? -----------------------
ident <- fread("is_out/IS_identity.tsv")
ident[, `:=`(qg = tstrsplit(q, "\\|")[[1]], sg = tstrsplit(s, "\\|")[[1]],
             fam = tstrsplit(q, "\\|")[[5]])]
ident <- merge(ident, d[, .(qg = genome, qfocal = focal)], by = "qg")
ident <- merge(ident, d[, .(sg = genome, sfocal = focal)], by = "sg")

## Within-genome comparisons are the cleanest: copies inside one genome
## that are ~100% identical to each other cannot have been sitting there
## since the common ancestor of the genome's core genes.
within <- ident[qg == sg]
ggplot(within, aes(pident, fill = qfocal)) +
  geom_density(alpha = .4) + facet_wrap(~fam, scales = "free_y") +
  labs(x = "pairwise identity between IS copies within the same genome")
ggsave("IS_within_genome_identity.pdf", width = 10, height = 7)

## Quantify: proportion of within-genome pairs above 99% identity
recent <- within[, .(p_recent = mean(pident > 99), n = .N), by = .(qg, qfocal)]
print(wilcox.test(p_recent ~ qfocal, data = recent))

## Calibration: compare IS divergence to core-genome divergence for the
## same pair of genomes. IS pairs much less diverged than the core imply
## either recent transposition or horizontal acquisition of the IS itself.
core_d <- cophenetic(tree)
btw <- ident[qg != sg]
btw[, core_dist := core_d[cbind(qg, sg)]]
ggplot(btw, aes(core_dist, 100 - pident, colour = qfocal)) +
  geom_point(alpha = .1) + geom_smooth(method = "gam") +
  labs(x = "core-genome distance", y = "IS divergence (%)")
ggsave("IS_vs_core_divergence.pdf", width = 7, height = 5)

## ---------- Are the transposases themselves under constraint? -----
## An actively transposing family should show purifying selection on the
## transposase ORF alongside a tail of pseudogenised copies. A family
## that is uniformly degenerate is a fossil and cannot be your driver.
## (feed IS_seqs.fna per family to your existing dN/dS pipeline)
