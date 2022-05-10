def get_fast5_dir(wildcards):
    return samplesheet.at[wildcards.sample, "fast5dir"]
