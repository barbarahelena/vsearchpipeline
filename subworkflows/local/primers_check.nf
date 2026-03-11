//
// Check primersheet and emit primer sequences
//

workflow PRIMERS_CHECK {
    take:
    primersheet // path: /path/to/primers.csv

    main:
    primersheet
        .map { csv ->
            def rows = csv.splitCsv(header: true, strip: true)
            if (rows.size() != 1) {
                error "Primersheet must contain exactly one row of primers, found ${rows.size()} rows."
            }
            def row = rows[0]
            // Accept both 'forward_primer'/'reverse_primer' and 'fwd_primer'/'rev_primer'
            def fwd = row['forward_primer'] ?: row['fwd_primer']
            def rev = row['reverse_primer'] ?: row['rev_primer']
            if (!fwd || !rev) {
                error "Primersheet must have columns 'forward_primer' and 'reverse_primer' (or 'fwd_primer'/'rev_primer')."
            }
            if (!fwd.matches('[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+')) {
                error "Forward primer '${fwd}' is not a valid IUPAC DNA sequence."
            }
            if (!rev.matches('[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+')) {
                error "Reverse primer '${rev}' is not a valid IUPAC DNA sequence."
            }
            return [ forward: fwd.trim(), reverse: rev.trim() ]
        }
        .set { primers }

    emit:
    primers // channel: [ val(map with forward/reverse keys) ]
}

