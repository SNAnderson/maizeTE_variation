rule nucmer:
    input:
        lambda w: config['samples']["{}".format(w.query)]['fasta'],
        lambda w: config['samples']["{}".format(w.subj)]['fasta']
    output:
        "{run}/{query}_{subj}-nucmer.coords"
    run:
        shell("/panfs/roc/msisoft/mummer/3.23/nucmer -c 5000 -p {wildcards.run}/{wildcards.query}_{wildcards.subj}-nucmer {input[0]} {input[1]}")
        shell("/panfs/roc/msisoft/mummer/3.23/show-coords -r -c -l {wildcards.run}/{wildcards.query}_{wildcards.subj}-nucmer.delta > {output}")

rule parse_nucmer:
    input:
        "{run}/{query}_{subj}-nucmer.coords"
    output:
        "{run}/{query}_{subj}_query-nucmer.coords_parsed.txt",
        "{run}/{query}_{subj}_subj-nucmer.coords_parsed.txt"
    run:
        shell("perl munge/parse-nucmer.pl {input} {output[0]} {output[1]}")

rule add_nucmer:
    input:
        "{run}/{query}-{subj}_key-wOrtho.txt",
        "{run}/{query}_{subj}_query-nucmer.coords_parsed.txt",
        "{run}/{query}_{subj}_subj-nucmer.coords_parsed.txt",
        "{run}/{query}-{subj}.last"
    output:
        "{run}/{query}-{subj}_key-beta.txt"
    run:
        shell("python3 munge/incorporate_nucmer.py --key {input[0]} --last {input[3]} --query {wildcards.query} --subj {wildcards.subj} --mummer_query {input[1]} --mummer_subj {input[2]} --output {output}")
