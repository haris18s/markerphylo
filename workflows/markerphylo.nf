/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_markerphylo_pipeline'
include { NCBIGENOMEDOWNLOAD     } from '../modules/nf-core/ncbigenomedownload/main'
include { CHECKM2_PREDICT        } from '../modules/nf-core/checkm2/predict/main'
include { PYRODIGAL              } from '../modules/nf-core/pyrodigal/main'
include { HMMER_HMMSEARCH        } from '../modules/nf-core/hmmer/hmmsearch/main'
include { GUNZIP                 } from '../modules/nf-core/gunzip/main'
include { SELECT_MARKERS         } from '../modules/local/select_markers/main'
include { MERGEMARKERS           } from '../modules/local/mergemarkers/main'
include { MAFFT_ALIGN            } from '../modules/nf-core/mafft/align/main'
include { CLIPKIT                } from '../modules/nf-core/clipkit/main'
include { SEQKIT_CONCAT          } from '../modules/nf-core/seqkit/concat/main'
include { IQTREE                 } from '../modules/nf-core/iqtree/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MARKERPHYLO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // read file 1 Parse file

    ch_samplesheet.branch { meta, accession, fasta ->
        download: accession
        local: fasta
        }
        .set { input_branch } // equal to setting a variable
    // input_branch.local.view {"local: $it"}
    // input_branch.download.view {"download: $it"}

    // parse accession file to fownload genomes in ncbi
    acc_file = input_branch.download
        .map { _meta, accession, _fasta -> accession}
        .collectFile(name: 'accession_file.txt', newLine:true)

    // acc_file.view { f -> "Accession file ${f.text}"  }

    // Step 2: Download genomes

    NCBIGENOMEDOWNLOAD(
        [id: 'ncbi'],
        acc_file,
        [],
        'bacteria'
    )
    NCBIGENOMEDOWNLOAD.out.fna.view { genome -> "Downloded genome: genome ${genome}" }

    ch_downloaded_files = NCBIGENOMEDOWNLOAD.out.fna
        .flatMap { _meta, files ->
            (files instanceof List ? files : [files]).collect { f ->
            def acc = f.name.tokenize('_')[0..1].join('_') // glues the downloaded accession name
            [acc, f]
            }
        }

    ch_acc_meta = input_branch.download.map {meta, accession, _fasta -> [accession, meta]}

    donwload_ch = ch_acc_meta
        .join(ch_downloaded_files)
        .map { _accession, meta, fasta -> [meta, fasta] }


    // Format local genomes as [meta, fasta]
    input_branch.local
        .map{ meta, _accession, fasta ->
            [meta, file(fasta)]

    }.set{ local_ch }


    // local_ch.view{"here ${it}"}

    //Combine donwload & local genomes
    genomes_ch = local_ch.mix(donwload_ch)
    // genomes_ch.view()

    // Step 3: Genome QC, check if not skip_qc is enabled
    if (!params.skip_qc) {
        // channel meta for checkm2 db
        checkm2_db_ch = channel.value([
            [id: 'checkm2_db'],
            file(params.checkm2_db, checkIfExists:true)
        ])

        CHECKM2_PREDICT(genomes_ch,checkm2_db_ch)
    }

    // Step 4: Gene calling using Pyrodical
    PYRODIGAL(genomes_ch, 'gff')

    GUNZIP(PYRODIGAL.out.faa)

    proteins_ch = GUNZIP.out.gunzip

    // Step 5: HMM-search of marker profiles against each genome's proteins
    // create channel for input
    hmm_input_ch = proteins_ch.map {meta, faa ->
        [
            meta, //sample info
            file(params.markers_hmm, checkIfExists: true), //Marker HMM file
            faa, //genome's protein
            false, //write align: not need , can eb large
            true, //write targer: per protein hit tablee
            true // write_domain: per-domain table
        ]
    }

    // Run HMMEer_searcch
    HMMER_HMMSEARCH(hmm_input_ch)

    HMMER_HMMSEARCH.out.target_summary.view { t -> "Hits: ${t}" }

    //Select_markers.py need 2 files per genome, come from 2 channels,
    // 1) proteins_ch :[meta, proteins.faa] 2) HMMER_HMMSEARCH.out.targer_summary ->. [meta, hits.tbl.gz]

    select_input_ch = proteins_ch.join(HMMER_HMMSEARCH.out.target_summary)

    SELECT_MARKERS(select_input_ch)

    SELECT_MARKERS.out.markers.view {m -> "markers: ${m}"}

    // merge script to merge  a per-genome file FASTA into file per marker with all genomes
    MERGEMARKERS(SELECT_MARKERS.out.markers.map {meta, f -> f}.collect() )

    markers_ch = MERGEMARKERS.out.markers.flatten() // flatten turns one list of files into seperate items
        .map{ f -> [ [id:f.baseName], f ]} // f.baseName is the file without its extension, so each alignemnt is named after each marker

    MAFFT_ALIGN(
        markers_ch,
        [[:], []], // first [:] is the empty meta and second is empty file
        [[:], []],
        [[:], []],
        [[:], []],
        [[:], []],
        false //keep output unccompressed for trimming
    )

    // trimming with clipkit
    CLIPKIT(
        MAFFT_ALIGN.out.fas,
        'fasta', //output format: plain FASTa for the NEXT concat step
        [] // auxilalry file: not needed
    )

    //Collect all trimmed alignemntd into one item for concatenation
    concat_input_ch = CLIPKIT.out.clipkit
        .map{ meta, f -> f } // keep only the files
        .collect() // all markers together
        .map {files -> [ [id: 'concatenated'], files] }// one meta for the batch

    // Concat multiple uncompressed sequence files togethersequences with seqkit/concat
    SEQKIT_CONCAT(
        concat_input_ch

    )

    //IQ tree needs at least 3 genomes
    IQTREE(
        // IQTREE expects 3 values (meta, alignment, starting tree), but SEQKIT_CONCAT outputs 2 (meta, alignment),
         SEQKIT_CONCAT.out.fastx.map{ meta, concatenated_fasta -> [meta, concatenated_fasta, [] ] },
         [], [], [], [], [] , [], [], [], [], [], [], [] // 12 optional inputs
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'markerphylo_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    // MULTIQC(
    //     ch_multiqc_files.flatten().collect().map { files ->
    //         [
    //             [id: 'markerphylo'],
    //             files,
    //             multiqc_config
    //                 ? file(multiqc_config, checkIfExists: true)
    //                 : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
    //             multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
    //             [],
    //             [],
    //         ]
    //     }
    // )
    emit:multiqc_report = channel.empty() // MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
