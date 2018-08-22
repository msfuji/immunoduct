bin_dir=config["env_dir"]+"/bin/"

rule all:
    input:
        "output/immunoduct.gct"

################################################################################
rule merge_input:
    input:
        config["expression"]
    output:
        "input/expression.gct",
        dir="input/"
    log:
        "log/merge_input/"
    shell:
        bin_dir+"Rscript scripts/merge_input.R {output.dir} {input}"

rule cyt:
    input:
        "input/expression.gct"
    output:
        "signature/cyt.gct"
        dir="signature/"
    log:
        "log/cyt/"
    shell:
        bin_dir+"Rscript scripts/cyt.R {input} {output.dir}"

################################################################################
rule merge_output:
    input:
        cyt="signature/cyt.gct"
    output:
        "output/immunoduct.gct"
        dir="output/"
    log:
        "log/merge_output/"
    shell:
        bin_dir+"Rscript scripts/merge_output.R {output.dir} {input}"
