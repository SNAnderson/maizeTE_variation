Processed data files for TE variation manuscript 
Compiled by Sarah Anderson (permanent email: asarahanderson@gmail.com)
____________________________________________

FilteredTE annotations
PH207.structuralTEv2.2018.10.30.filteredTE.gff3Mo17.structuralTEv2.2018.12.28.filteredTE.gff3W22.structuralTEv2.10.08.2018.filteredTE.gff3B73.structuralTEv2.2018.12.20.filteredTE.gff3


The TE annotation files include metadata including LTR similarity scores for full-length LTR retrotransposons as output by LTRharvest, terminal inverted repeat (TIR) lengths and number of mismatches for TIR elements, and target site duplication (TSD) lengths for TIRs, LINEs, SINEs, and soloLTRs. TIR elements were filtered to remove annotations which were nested within another TIR element when the disjoined length of the outer element was less than twice the outer element’s TSD length since in many cases this pattern resulted from ambiguity between 2 bp TIR/TSD patterns. The TE annotation was further filtered to remove ~ 1800 TEs from each genome where the TE coordinates completely overlapped a gene with synteny to rice or sorghum, suggesting a false positive annotation. 

____________________________________________

Disjoined TE annotation
PH207.structuralTEv2.10.30.2018.filteredTE.disjoined.gff3Mo17.structuralTEv2.2018.12.28.filteredTE.disjoined.gff3W22.structuralTEv2.10.08.2018.filteredTE.disjoined.gff3B73.structuralTEv2.2018.12.20.filteredTE.disjoined.gff3

The TE annotation file was disjoined to resolve nested TE insertions and to create a file where each bp of the genome is assigned to only the element contributing the DNA sequence of that region. 

____________________________________________

Gene keys
2018-08-21_w22-ph207_key-beta.txt2018-08-21_b73-ph207_key-beta.txt2018-08-21_b73-mo17_key-beta.txt2018-08-21_b73-w22_key-beta.txt2018-08-21_w22-mo17_key-beta.txt2018-08-21_ph207-mo17_key-beta.txt

These files were created by Alex Brohammer. A cross-reference of homologous maize genes between assemblies was produced using multiple complementary approaches used in iterative fashion. First, a local version of the SynMap pipeline (Lyons et al. 2008), was used in order to identify stretches of collinear genes in pairwise comparisons of genomes. This pipeline uses the LAST aligner version 963 (Kiełbasa et al. 2011) to identify hits between CDS sequence from each genome and then incorporates DAGchainer to identify “chains” of collinear hits. The Lastal algorithm was run using default parameters, however hits were later filtered to have a c-score of 0.10 before being supplied to DAGchainer (Haas et al. 2004). DAGchainer was ran using an e-value cutoff of 0.05, allowing a maximum distance of 10 genes between matches (-D), requiring a minimum of 12 aligned pairs per chain (-A), and a gap distance of 7 (-g).

The Nucmer alignment algorithm within MUMmer version 3.32 (Kurtz et al. 2004) was then used to perform whole-genome alignments between homologous chromosomes of each assembly in pairwise fashion (-c 5000). The ‘show-coords’ command was used to filter any alignments not included in the longest ascending subset (-g flag) and structural gene annotation files are used to identify genic positions in the alignments. The Nucmer-based gene assignments were cross-referenced with the assignments obtained using the SynMap pipeline. Any gene assignment unique to Nucmer was required to be no further than 3 genes upstream and downstream from the nearest SynMap based gene assignment. This allowed genes that are split across multiple gene models, genes affected by local rearrangements, and genes missed by the SynMap approach to be recovered.

A third approach using the OrthoFinder clustering algorithm version 2.2.7 (Emms and Kelly 2015) was also used to identify homologous genes. OrthoFinder was ran in manual mode by first performing blastp (Altschul et al. 1997), searches across genomes requiring an e-value of 1e-3. The clustering algorithm was then run with the ‘-og’ option to output groups. Each orthogroup cluster was scanned to identify genes from different genotypes that were located on homologous chromosomes. Similar to the nucmer-based gene assignments, collinearity between genes meeting these criteria was assessed using the SynMap and Nucmer assignments. Any assignment unique to this method was required to be no further than 8 genes upstream and downstream from the nearest assignment.

____________________________________________

All contrasts
PH207_filteredTE_allContrasts_1Feb19.txtMo17_filteredTE_allContrasts_1Feb19.txtW22_filteredTE_allContrasts_1Feb19.txtB73_filteredTE_allContrasts_1Feb19.txt

These files summarize the shared, non shared, and unresolved calls for all TEs compared with all other genomes. Inferred coordinates in the target genome are listed where possible, along with the number of genes where the TE is present, absent, or unresolved.

____________________________________________

TE attributes
TE_attributes_4genomes_1Feb19.txt

TE attributes including (where applicable) TSD length, TIR length, mismatches in TIR, set (TE type), LTR similarity,  family, superfamily, order, disjoined length. These attributes are extracted from the annotation files. 

____________________________________________

Resolved non-redundant TE set
non-redundant_TEs_4genomes_1Feb19.txt

A non-redundant TE set was created for all elements that could be resolved as shared or non-shared compared to all other genome assemblies. This set was created in an iterative manner, starting with all B73 TEs that could be resolved compared to W22, Mo17, and PH207. Resolved W22 TEs that were non-shared with B73 and those that were defined as missing annotations in B73 were then added. Mo17 TEs called non-shared with both B73 and W22 and those that were shared with either B73 or W22 but were defined as missing annotations where shared were added. Finally, fully resolved TEs unique to PH207 and those that were shared with any combination of B73, W22, and Mo17 but were defined as missing annotations where present were added. This resulted in a set of 510,533 TEs present in at least one genome (Supplemental dataset X). The source of each TE annotation can be determined by the genome identifier in the TE name, labeled as Zm00001d for B73, Zm00004b for W22, Zm00014a for Mo17, and Zm00008a for PH207. Disjoined TE length and LTR similarity for the non-redundant TE set were extracted from the genome which was the source of the annotation. TEs were defined as shared when present in all four genomes and as variable when called non-shared with at least one other genome. For family and superfamily-level analyses, the non-redundant TE set was used to summarize TE variability across genomes. This means that in many cases the number of members analyzed for a given family is often less than the total number of annotated elements. Families unique to one genome were required to have annotated copies in only one genome and resolved members in only one genome. Resolved and annotated members in these families are listed in Table S3. 

____________________________________________

TE expression across development
Walley_expresslion_18Jan19.txt

TE family RPM values for Walley et al. 2016 developmental atlas (B73). 








