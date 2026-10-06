# Evidence carried into downstream AMR analysis

metaMARS prepares assemblies, genome bins, quality records and taxonomy. These
outputs support subsequent AMR annotation and pathway analysis; preparation alone
does not establish resistance phenotypes, gene transfer, or complete pathways.

## Canonical contigs and genome membership

Keep one canonical contig collection for each assembly, including contigs shorter
than the binning threshold. A membership table links exact contig identifiers to
their genome bins. Sequence identifiers and their sample or assembly provenance
must remain stable across coverage, membership and annotation tables.

Annotate the canonical contigs once in later evidence layers, then aggregate
annotations using membership. A gene on a binned contig is one observation with
both contig and genome context, not two independent AMR observations. In an isolate
assembly, the assembly represents the isolate genome; binning is unnecessary.

Unbinned contigs, short contigs and contigs in excluded bins remain available for
later annotation. An unbinned contig has unresolved genome membership. Do not
concatenate all unbinned contigs into a synthetic genome for CheckM2 or GTDB-Tk.

## Quality determines eligibility, not retention

CheckM2 reports genome completeness and contamination. The module's quality gates
determine which genomes are eligible for population comparison. Excluded genomes
and their sequences remain in the evidence collection with quality values and
exclusion status. A high CheckM2 score does not establish the correct host of
every accessory contig in a bin.

MetaQUAST evaluates assembly characteristics. Reference-dependent metrics require
suitable references and should be interpreted with their reference provenance.
Assembly continuity and estimated genome completeness answer different questions.

An AMR gene missing from an incomplete genome is **not detected**, rather than
confirmed absent. Contamination can also create apparent gene combinations from
different organisms. Downstream pathway comparisons should retain these quality
limits and distinguish a single-genome pathway from a community-level combination.

## Taxonomy and comparison populations

GTDB-Tk classifies genomes, using MAG or isolate assembly FASTA inputs. Its
documentation reports validation on genomes estimated to be at least 50% complete
and at most 10% contaminated. These are useful eligibility criteria, not proof
that every qualifying genome is accurately classified. Keep classification
warnings, the complete lineage, and the GTDB reference release.

The `--tax-rank` setting chooses domain, phylum, class, order, family, genus or
species for grouping. A stable group identifier is derived from the complete
lineage through that rank, so identically named taxa under different ancestors
remain distinct. Genomes missing the selected rank receive `unresolved` status and
an empty group identifier. They do not form a shared unknown population. Quality
exclusions receive `quality_excluded` status even when taxonomy is available.

GTDB-Tk taxonomy does not measure abundance. Counts of recovered MAGs are not
community relative abundance, because assembly, coverage and genome recovery
affect those counts. Abundance estimation requires a separately defined mapping
and normalization procedure. Likewise, sharing a species or genus assignment
does not establish strain identity.

Retain biological sample metadata separately from taxonomy: condition or cohort,
subject, site, time point, and replicate relationships where available. Comparisons
should respect these experimental units. Technical replicates and repeated genome
recovery must not silently become independent biological observations. If later
genome dereplication creates representatives, preserve every original genome's
sample association and representative membership.

## Host attribution and information loss

Short-read assembly and binning can miss mobile elements and AMR loci. Preserve
cleaned reads, assembly graphs, contig coverage and mapping provenance for later
evidence layers. Additional read-level AMR searches can capture signals absent
from assemblies, but read matches alone do not assign a host or genomic context.

Host associations should retain the evidence and uncertainty behind each link:
bin membership, an assembly connection, read-pair support, co-abundance, or another
validated method. Mobile-element copy number and sequence composition can differ
from the host chromosome. Co-abundance or a shared AMR sequence alone must not be
treated as proof of a host, physical linkage, or horizontal transfer. Keep candidate
associations separate from established membership and retain unresolved links.

References: [CheckM2 documentation](https://github.com/chklovski/CheckM2),
[GTDB-Tk workflow documentation](https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html),
[QUAST/MetaQUAST manual](https://quast.sourceforge.net/docs/manual.html), and
[Maguire et al. on mobile-element losses in MAG workflows](https://pmc.ncbi.nlm.nih.gov/articles/PMC7660262/).
