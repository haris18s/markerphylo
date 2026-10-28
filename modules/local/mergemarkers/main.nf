
process MERGEMARKERS {
    tag "all_genomes"
    label 'process_single'

    conda     "conda-forge::python=3.14.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'oras://community.wave.seqera.io/library/python:3.14.1--595559071f0b240b' :
    'community.wave.seqera.io/library/python:3.14.1--6e31340b9584a1fe' }"


    input:
        path(marker_files)


    output:
        path("markers/*.faa"), emit: markers
        tuple val("${task.process}"), val('mergemarkers'), eval("python --version sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    merge_markers.py \\
        --inputs ${marker_files} \\
        --outdir markers
        ${args}

    """

}
