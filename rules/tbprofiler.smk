
rule tbprofiler:
    input:
        unpack(get_trimmed_reads_dict) 
    output: 
        "output/tbprofiler_" + tbprofiler_caller + "/results/{sample}.results.json",
        "output/tbprofiler_" + tbprofiler_caller + "/results/{sample}.results.txt"
    params: 
        sample_name="{sample}",
        working_dir="output/tbprofiler_" + tbprofiler_caller,
        version=config["tool_params"]["tbprofiler"]["version"],
        caller=config["tool_params"]["tbprofiler"]["caller"]
    log:
        "log/tbprofiler_" + tbprofiler_caller + "/{sample}.log"
    threads: 8
    # conda:
    #     "../envs/tbprofiler.yaml"
    shell:"""
        module purge
        module load tbprofiler/4.3.0_python3
        tb-profiler profile --threads {threads} \
            --no_trim \
            --caller {params.caller} \
            --read1 {input.r1} --read2 {input.r2} \
            --prefix {params.sample_name} \
            --dir {params.working_dir} --txt 2> {log}
    """

rule profiler_merge_list:
    input: 
        expand("output/tbprofiler_" + tbprofiler_caller + "/results/{sample}.results.json", sample=samples.index)
    output: 
        f"output/tbprofiler_{tbprofiler_caller}/merge/{merged_sample_name}_samples_to_merge.txt"
    params: 
        samples=samples.index
    run:
        with open(output[0], "w") as f:
            for sample in params.samples:
                f.write("{}\n".format(sample))

rule merge_tbprofiler:
    input:
        profiler_res=get_profiler_res,
        samples_to_merge=rules.profiler_merge_list.output
    output: 
        f"output/tbprofiler_{tbprofiler_caller}/merge/{merged_sample_name}.json",
        f"output/tbprofiler_{tbprofiler_caller}/merge/{merged_sample_name}.variants.txt",
        f"output/tbprofiler_{tbprofiler_caller}/merge/{merged_sample_name}.txt"
    params: 
        working_dir=f"output/tbprofiler_{tbprofiler_caller}/results",
        out_name=f"output/tbprofiler_{tbprofiler_caller}/merge/{merged_sample_name}",
        version=config["tool_params"]["tbprofiler"]["version"]
    log:
        f"log/tbprofiler_merge/{merged_sample_name}_tbprofiler_{tbprofiler_caller}.log"
    threads: 8
    # conda:
    #     "../envs/tbprofiler.yaml"
    shell:"""
        module purge
        module load tbprofiler/4.3.0_python3
        tb-profiler collate --dir {params.working_dir} --samples {input.samples_to_merge} --prefix {params.out_name} 2> {log}
    """

rule update_tbprofiler:
    input:
        json=f"output/tbprofiler_{tbprofiler_caller}/merge/{merged_sample_name}.json",
        variants=f"output/tbprofiler_{tbprofiler_caller}/merge/{merged_sample_name}.variants.txt",
        summary=f"output/tbprofiler_{tbprofiler_caller}/merge/{merged_sample_name}.txt"
    output: 
        json=f"results/{merged_sample_name}/{merged_sample_name}_tbprofiler.json",
        variants=f"results/{merged_sample_name}/{merged_sample_name}_tbprofiler.variants.txt",
        summary=f"results/{merged_sample_name}/{merged_sample_name}_tbprofiler.summary.txt"
    params: 
        isodate=isodate,
        version=config["tool_params"]["tbprofiler"]["version"],
        caller=config["tool_params"]["tbprofiler"]["caller"]
    log:
        f"log/tbprofiler_merge/{merged_sample_name}_tbprofiler_{tbprofiler_caller}_update.log"
    # conda:
    #     "../envs/tbprofiler.yaml"
    shell:"""
        cp {input.json} {output.json}
        cp {input.variants} {output.variants}
        cp {input.summary} {output.summary}
    """


