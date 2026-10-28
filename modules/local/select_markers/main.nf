process SELECT_MARKERS {
    tag "${meta.id}"
    label 'process_single'

    conda     "conda-forge::python=3.14.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'oras://community.wave.seqera.io/library/python:3.14.1--595559071f0b240b' :
    'community.wave.seqera.io/library/python:3.14.1--6e31340b9584a1fe' }"

    input:
        tuple val(meta), path(proteins), path(hits) // include both protin file and hit table to match ids with the protein seqs

    output:
        tuple val(meta), path("${prefix}_markers.faa"), emit: markers
        tuple val(meta), path("${prefix}_marker_summary.tsv"), emit: summary
        tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), emit:versions_python, topic:versions


    when:
    task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"

        """
        select_markers.py \\
            --proteins ${proteins} \\
            --hits ${hits} \\
            --genome ${meta.id} \\
            --outfile ${prefix}_markers.faa \\
            --summary ${prefix}_marker_summary.tsv \\
            ${args}
        """

}