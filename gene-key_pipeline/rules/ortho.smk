rule prepare_ortho:
    input:
        "{query}-{subj}/{query}-{subj}_key.txt",
        "{query}-{subj}/{query}.repcds.bed",
        "{query}-{subj}/{subj}.repcds.bed"
    output:
        "{query}-{subj}/orthofinder/Species0.fa",
        "{query}-{subj}/orthofinder/Species1.fa",
        "{query}-{subj}/orthofinder/SequenceIDs.txt",
        "{query}-{subj}/orthofinder/SpeciesIDs.txt"
    params:
        dir = "{query}-{subj}/orthofinder"
    shell:
        """
        python3 munge/prepare-orthofinder.py --query {wildcards.query} --querybed {input[1]} --subj {wildcards.subj} --subjbed {input[2]} --dir {params.dir}
        touch {output[0]}
        touch {output[1]}
        touch {output[2]}
        touch {output[3]}
        """

rule make_blast_db:
    input:
        "{run}/orthofinder/Species0.fa",
        "{run}/orthofinder/Species1.fa",
        "{run}/orthofinder/SequenceIDs.txt",
        "{run}/orthofinder/SpeciesIDs.txt"
    output:
        "{run}/orthofinder/{run}"
    shell:
        """
        makeblastdb -dbtype prot -parse_seqids -in {input[0]}
        makeblastdb -dbtype prot -parse_seqids -in {input[1]}
        touch {output}
        """

rule prepare_blast:
    input:
        "{run}/orthofinder/{query}-{subj}"
    output:
        "{run}/orthofinder/{query}-{subj}-commands.txt"
    params:
        querydb="{run}/orthofinder/{query}.fa",
        subjdb="{run}/orthofinder/{subj}.fa"
    shell:
        """
        echo "blastp -outfmt 6 -evalue 0.0005 -query {wildcards.query}-{wildcards.subj}/orthofinder/{wildcards.query}.fa -db {params.querydb} -out {wildcards.query}-{wildcards.subj}/orthofinder/Blast0_0.txt" >> {output}
        echo "blastp -outfmt 6 -evalue 0.0005 -query {wildcards.query}-{wildcards.subj}/orthofinder/{wildcards.query}.fa -db {params.subjdb} -out {wildcards.query}-{wildcards.subj}/orthofinder/Blast0_1.txt" >> {output}
        echo "blastp -outfmt 6 -evalue 0.0005 -query {wildcards.query}-{wildcards.subj}/orthofinder/{wildcards.subj}.fa -db {params.querydb} -out {wildcards.query}-{wildcards.subj}/orthofinder/Blast1_0.txt" >> {output}
        echo "blastp -outfmt 6 -evalue 0.0005 -query {wildcards.query}-{wildcards.subj}/orthofinder/{wildcards.subj}.fa -db {params.subjdb} -out {wildcards.query}-{wildcards.subj}/orthofinder/Blast1_1.txt" >> {output}
        """

rule run_orthofinder_blast:
    input:
        "{run}/orthofinder/{query}-{subj}-commands.txt"
    output:
        "{run}/orthofinder/Blast_{query}_{subj}-orthofinder.txt"
    shell:
        """
        module load parallel
        module load orthofinder
        #export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY"
        env_parallel --jobs 4 \
        --sshloginfile $PBS_NODEFILE \
        --workdir $PWD \
        --block 1k -k --recstart '>' \
        < {input}
        touch {output}
        """

rule run_orthofinder_cluster:
    input:
        "{run}/orthofinder/Blast_{query}_{subj}-orthofinder.txt"
    output:
        "{run}/orthofinder/Blast_{query}_{subj}-orthofinder-groups.txt"
    params:
        ortho_dir = "{run}/orthofinder"
    run:
        shell("orthofinder -b {params.ortho_dir} -t 8 -og")
        shell("touch {output}")

rule add_orthofinder:
    input:
        "{run}/orthofinder/Blast_{query}_{subj}-orthofinder-groups.txt",
        "{run}/{query}-{subj}_key.txt"
    params:
        ortho_csv="{run}/orthofinder/Orthogroups.csv"
    output:
        "{run}/{query}-{subj}_key-wOrtho.txt"
    run:
        shell("python3 munge/incorporate-orthofinder.py --infile {params.ortho_csv} --query {wildcards.query} --subj {wildcards.subj} --key {input[1]} --outfile {output}")
