bin_dir=config["env_dir"]+"/bin/"

rule all:
    input:
        "output/immunoduct.gct"

################################################################################
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

################################################################################
rule merge_output:
    input:
        cyt="signature/cyt.gct"
    output:
        file="output/immunoduct.gct",
        dir="output/"
    log:
        "log/merge_output/"
    shell:
        bin_dir+"Rscript scripts/merge_gct.R row {output.file} {input}"
