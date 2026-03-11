//
// Check input samplesheet and get read channels
//

include { fromSamplesheet } from 'plugin/nf-validation'

workflow INPUT_CHECK {
    take:
    _samplesheet // unused: nf-validation reads params.input directly

    main:
    channel.fromSamplesheet('input')
        .map { meta, fastq_1, fastq_2 ->
            // meta already contains: id, nucl_acid_conc (from schema meta annotations)
            def new_meta = meta + [ single_end: !fastq_2 ]
            def reads = fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ]
            [ new_meta, reads ]
        }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}
