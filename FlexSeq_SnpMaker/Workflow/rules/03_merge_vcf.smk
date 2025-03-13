# This rule will merge the bcf files into one big file
rule index_bcf:
    input:
        bcf_file=vcfs + "{sample}.bcf"
    output:
        csi_file=vcfs + "{sample}.bcf.csi"
    conda:
        "../Conda_Envs/bcftools.yaml"
    params:
        bcf_path= vcfs  # Path to the directory containing BCF file
    shell:
        """
        bcftools index --csi {input.bcf_file}
         """

# Can not get the merge to work
# rule merge:
#     input:
#         calls=expand(vcfs + "{sample}.bcf", sample=SAMPLES),
#         index=expand(vcfs + "{sample}.bcf.csi", sample=SAMPLES)
#     output:
#         output
#     log:
#         "logs/merge.log"
#     conda:
#         "../Conda_Envs/bcftools.yaml"
#     shell:
#         """
#         bcftools merge {input.calls} -o {output}
#         """

