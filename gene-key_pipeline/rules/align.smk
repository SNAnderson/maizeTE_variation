rule get_rep_cds:
    input:
        lambda w: config['samples']["{}".format(w.query)]['cds']
    output:
        "{run}/{query}.repcds.fa"
    params:
        feature = lambda w: config["repcds"]["{}".format(w.query)],
        genotype = "{query}"
    run:
        shell("~/miniconda3/envs/py27/bin/python munge/get-rep-cds.py --genotype {params.genotype} --outfile {output} --extract_id {params.feature} --mode fasta")

rule get_bed:
    input:
        "{run}/{query}.repcds.fa"
    output:
        "{run}/{query}.repcds.bed"
    params:
        feat=lambda w: config["repcds"]["{}".format(w.query)],
        genotype="{query}"
    run:
        shell("~/miniconda3/envs/py27/bin/python munge/get-rep-cds.py --genotype {params.genotype} --outfile {output} --extract_id {params.feat} --mode bed")

rule build_last_index:
    input:
        lambda w: config['samples']["{}".format(w.subj)]['cds']
        #config['samples']["{subj}"]['bed']
    output:
        "{run}/lastdb/{subj}"
    run:
        shell("lastdb {output} {input}")
        shell("touch {output}")

rule last:
    input:
        "{run}/{query}.repcds.fa",
        "{run}/lastdb/{subj}"
    output:
        "{run}/{query}-{subj}.last"
    params:
        extra=config["last"]["extra"],
        threads=config["last"]["threads"]
    log:
        "logs/last/{run}.log"
    threads:
        config["last"]["threads"]
    run:
        shell("lastal -P {params.threads} {params.extra} {input[1]} {input[0]} > {output}")

# Move bed params to input if create_bed rule implemented
rule blast2raw:
    input:
        "{run}/{query}-{subj}.last",
        "{run}/{query}.repcds.bed",
        "{run}/{subj}.repcds.bed"
    output:
        "{run}/{query}-{subj}.filtered"
    params:
        tandemn=0,
        cscore=0.1
    log:
        "logs/last/{run}.log"
    shell:
        "~/miniconda3/envs/py27/bin/python quota-alignment/scripts/blast_to_raw.py {input[0]} --qbed {input[1]} --sbed {input[2]} --tandem_Nmax {params.tandemn} --cscore {params.cscore} --write-filtered-blast > {output} 2>>{log}"

rule dag_format:
    input:
        "{run}/{query}-{subj}.filtered"
    output:
        "{run}/{query}-{subj}.filtered.dag"
    params:
        q="Zmays-{query}",
        s="Zmays-{subj}",
        qbed="{run}/{query}.repcds.bed",
        sbed="{run}/{subj}.repcds.bed"
    log:
    shell:
        "~/miniconda3/envs/py27/bin/python munge/dag_tools.py -q {params.q} -s {params.s} -b {input} -c True --query-bed {params.qbed} --subj-bed {params.sbed} > {output}"

rule geneorder_conversion:
    input:
        "{run}/{query}-{subj}.filtered.dag"
    output:
        "{run}/{query}-{subj}.filtered.dag.go"
    params:
        qbed="{run}/{query}.repcds.bed",
        sbed="{run}/{subj}.repcds.bed"
    shell:
        "~/miniconda3/envs/py27/bin/python munge/gene_order.py {input} {output} --query-bed {params.qbed} --subj-bed {params.sbed}"

rule dagchainer:
    input:
        "{run}/{query}-{subj}.filtered.dag.go"
    output:
        "{run}/{query}-{subj}.filtered.dag.go.aligncoords"
    params:
        evalue=0.05,
        maxMatchDist=10,
        minAlignedPairs=12,
        gapDistance=7
    shell:
        "~/miniconda3/envs/py27/bin/python dagtools_bp/dag_chainer.py -i {input} -E {params.evalue} -D {params.maxMatchDist} -A {params.minAlignedPairs} -g {params.gapDistance} > {output}"

rule quota_align:
    input:
        "{run}/{query}-{subj}.filtered.dag.go.aligncoords"
    output:
        "{run}/{query}-{subj}.filtered.dag.go.aligncoords.qa.merged"
    params:
        maxDistance=0
    shell:
        "perl munge/quota_align_merge.pl --infile {input} --outfile {output} --max_distance {params.maxDistance}"

rule syntenic_depth:
    input:
        "{run}/{query}-{subj}.filtered.dag.go.aligncoords.qa.merged"
    output:
        "{run}/{query}-{subj}.filtered.dag.go.aligncoords.qa.merged.qa.filtered.cov"
    params:
        depth1=1,
        depth2=1,
        depthoverlap=40,
        prefix="{run}/{query}-{subj}.filtered.dag.go.aligncoords.qa.merged.qa"
    shell:
        "perl munge/quota_align_coverage.pl --infile {input} --depth_ratio_org1 {params.depth1} --depth_ratio_org2 {params.depth2} --depth_overlap {params.depthoverlap} --prefix {params.prefix} --outfile {output}"

rule coordinate_conversion:
    input:
        "{run}/{query}-{subj}.filtered.dag.go.aligncoords.qa.merged.qa.filtered.cov"
    output:
        "{run}/{query}-{subj}.filtered.dag.go.aligncoords.qa.merged.qa.filtered.cov.coor"
    params:
        qbed="{run}/{query}.repcds.bed",
        sbed="{run}/{subj}.repcds.bed"
    shell:
        "python munge/gene_order.py {input} {output} --query-bed {params.qbed} --subj-bed {params.sbed} --positional"

rule make_key:
    input:
        "{run}/{query}-{subj}.filtered.dag.go.aligncoords.qa.merged.qa.filtered.cov"
    output:
        "{run}/{query}-{subj}_key.txt"
    shell:
        """
        python3 munge/synmap-to-key.py --infile {input} --outfile {output}.tmp
        sort -k 1,1V -k 2,2n -k3,3n {output}.tmp > {output} && rm {output}.tmp
        """
