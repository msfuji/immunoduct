bin_dir=config["env_dir"]+"/bin/"

rule all:
    input:
        "output/immunoduct.gct",
        "cluster/cluster.txt"

rule merge_input:
    input:
        config["expression"]
    output:
        file="input/expression.gct",
        dir="input/"
    log:
        "log/merge_input/"
    shell:
        bin_dir+"Rscript scripts/merge_gct.R col {output.file} {input}"

rule merge_gene_sets:
    input:
        config["gene_set"]
    output:
        file="input/gene_sets.gmt",
        dir="input/"
    log:
        "log/merge_gene_sets/"
    shell:
        "cat {input} > {output.file}"

################################################################################

rule cyt:
    input:
        "input/expression.gct"
    output:
        file="signature/cyt.gct",
        dir="signature/"
    log:
        "log/cyt/"
    shell:
        bin_dir+"Rscript scripts/cyt.R {input} {output.file}"

rule goi:
    input:
        expr="input/expression.gct",
        name=config["goi"]
    output:
        file="goi/goi.gct",
        dir="goi"
    log:
        "log/goi/"
    shell:
        bin_dir+"Rscript scripts/filter_gct.R {input.expr} {input.name} {output.file}"

rule ssgsea:
    input:
        expr="input/expression.gct",
        gmt="input/gene_sets.gmt"
    output:
        file="signature/ssgsea.gct",
        dir="signature/"
    log:
        "log/ssgsea/"
    threads: 8
    shell:
        bin_dir+"Rscript scripts/ssgsea.R {threads} {input.expr} {input.gmt} {output.file}"

rule estimate:
    input:
        "input/expression.gct"
    output:
        file="signature/estimate.gct",
        dir="signature/"
    log:
        "log/estimate/"
    shell:
        bin_dir+"Rscript scripts/estimate.R {input} {output.file}"

rule epic:
    input:
        "input/expression.gct"
    output:
        file="cell/epic.gct",
        dir="cell/"
    log:
        "log/epic/"
    shell:
        bin_dir+"Rscript scripts/epic.R {input} {output.file}"

rule mcp_counter:
    input:
        "input/expression.gct"
    output:
        file="cell/mcp_counter.gct",
        dir="cell/"
    log:
        "log/mcp_counter/"
    shell:
        bin_dir+"Rscript scripts/mcp_counter.R {input} {output.file}"

rule xcell:
    input:
        "input/expression.gct"
    output:
        file="cell/xcell.gct",
        dir="cell/"
    log:
        "log/xcell/"
    shell:
        bin_dir+"Rscript scripts/xcell.R {input} {output.file}"

rule cibersort:
    input:
        "input/expression.gct"
    output:
        file="cell/cibersort.gct",
        dir="cell/"
    log:
        "log/cibersort/"
    shell:
        bin_dir+"Rscript scripts/cibersort.R {input} {output.file}"

################################################################################
rule merge_output:
    input:
        "signature/cyt.gct",
        "goi/goi.gct",
        "signature/ssgsea.gct",
        "signature/estimate.gct",
        "cell/epic.gct",
        "cell/mcp_counter.gct",
        "cell/xcell.gct",
        "cell/cibersort.gct"
    output:
        file="output/immunoduct.gct",
        dir="output/"
    log:
        "log/merge_output/"
    shell:
        bin_dir+"Rscript scripts/merge_gct.R row {output.file} {input}"

rule make_cluster_input:
    input:
        imm="output/immunoduct.gct",
        ann=config["annotation"]
    output:
        file="cluster/cluster.txt",
        dir="cluster/"
    log:
        "log/make_cluster_input"
    script:
        "scripts/make_cluster_input.py"
