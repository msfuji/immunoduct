bin_dir=config["env_dir"]+"/bin/"

rule all:
    input:
        "output/immunoduct.gct",
        "cluster/cluster.txt"

rule merge_input:
    input:
        config["expression"]
    output:
        file="input/expression.gct"
    log:
        "log/merge_input/"
    shell:
        config["rscript_path"]+" scripts/merge_gct.R col {output.file} {input}"

rule merge_gene_sets:
    input:
        config["gene_set"]
    output:
        file="input/gene_sets.gmt"
    log:
        "log/merge_gene_sets/"
    shell:
        "cat {input} > {output.file}"

################################################################################

rule cyt:
    input:
        "input/expression.gct"
    output:
        file="signature/cyt.gct"
    log:
        "log/cyt/"
    shell:
        config["rscript_path"]+" scripts/cyt.R {input} {output.file}"

rule goi:
    input:
        expr="input/expression.gct",
        name=config["goi"]
    output:
        file="goi/goi.gct"
    log:
        "log/goi/"
    shell:
        config["rscript_path"]+" scripts/filter_gct.R {input.expr} {input.name} {output.file}"

rule ssgsea:
    input:
        expr="input/expression.gct",
        gmt="input/gene_sets.gmt"
    output:
        file="signature/ssgsea.gct"
    log:
        "log/ssgsea/"
    threads: 8
    shell:
        config["rscript_path"]+" scripts/ssgsea.R {threads} {input.expr} {input.gmt} {output.file}"

rule estimate:
    input:
        "input/expression.gct"
    output:
        file="signature/estimate.gct"
    log:
        "log/estimate/"
    shell:
        config["rscript_path"]+" scripts/estimate.R {input} {output.file}"

rule epic:
    input:
        "input/expression.gct"
    output:
        file="cell/epic.gct"
    log:
        "log/epic/"
    shell:
        config["rscript_path"]+" scripts/epic.R {input} {output.file}"

rule mcp_counter:
    input:
        "input/expression.gct"
    output:
        file="cell/mcp_counter.gct"
    log:
        "log/mcp_counter/"
    shell:
        config["rscript_path"]+" scripts/mcp_counter.R {input} {output.file}"

rule xcell:
    input:
        "input/expression.gct"
    output:
        file="cell/xcell.gct"
    log:
        "log/xcell/"
    threads: 4
    shell:
        config["rscript_path"]+" scripts/xcell.R {threads} {input} {output.file}"

rule cibersort:
    input:
        "input/expression.gct"
    output:
        file="cell/cibersort.gct"
    log:
        "log/cibersort/"
    threads: 3
    shell:
        config["rscript_path"]+" scripts/cibersort.R {input} {output.file}"

rule impres:
    input:
        "input/expression.gct"
    output:
        file="signature/impres.gct"
    log:
        "log/impres/"
    shell:
        config["rscript_path"]+" scripts/impres.R {input} {output.file}"

rule log_average:
    input:
        expr="input/expression.gct",
        gmt="data/Ayers2017.gmt"
    output:
        file="signature/log_average.gct"
    log:
        "log/log_average/"
    shell:
        config["rscript_path"]+" scripts/log_average.R --qn {input.expr} {input.gmt} {output.file}"

################################################################################

def input_of_merge_output(wildcards):
    inputs=[
    "signature/cyt.gct",
    "goi/goi.gct",
    "signature/ssgsea.gct",
    "signature/estimate.gct",
    "signature/impres.gct",
    "signature/log_average.gct",
    "cell/epic.gct",
    "cell/mcp_counter.gct",
    "cell/xcell.gct",
    ]

    # check files necessary to run CIBERSORT
    has_cibersort=os.path.exists("scripts/CIBERSORT.R")
    has_lm22=os.path.exists("scripts/LM22.txt")
    if has_cibersort and has_lm22:
        inputs.append("cell/cibersort.gct")

    return inputs

rule merge_output:
    input:
        input_of_merge_output
    output:
        file="output/immunoduct.gct"
    log:
        "log/merge_output/"
    shell:
        config["rscript_path"]+" scripts/merge_gct.R row {output.file} {input}"

rule make_cluster_input:
    input:
        imm="output/immunoduct.gct",
        ann=config["annotation"]
    output:
        file="cluster/cluster.txt"
    log:
        "log/make_cluster_input/"
    script:
        "scripts/make_cluster_input.py"
