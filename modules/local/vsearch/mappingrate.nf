process VSEARCH_MAPPINGRATE {
    label 'process_single_low'

    input:
    path(count_table)
    path(mapping_stats)
    path(filter_stats)   // list of per-sample *.filter_stats.txt files

    output:
    path "mapping_rate_summary.tsv"  , emit: summary
    path "mapping_rate_overall.txt"  , emit: overall

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    import re
    import glob
    from collections import defaultdict

    # -------------------------------------------------------
    # 1. Parse overall mapping rate from vsearch stderr
    # -------------------------------------------------------
    overall_line = ""
    with open("${mapping_stats}") as fh:
        for line in fh:
            if "Matching unique query sequences" in line:
                overall_line = line.strip()
                break

    with open("mapping_rate_overall.txt", "w") as fh:
        if overall_line:
            fh.write(overall_line + "\\n")
            fh.write("(Note: this reflects mapping of unique dereplicated sequences, not raw reads)\\n")
        else:
            fh.write("Overall mapping stats not found.\\n")

    # -------------------------------------------------------
    # 2. Parse total reads per sample from fastqfilter stderr
    #    vsearch prints: "X sequences kept (of Y)"
    #    File names are: <sample>.filter_stats.txt
    # -------------------------------------------------------
    total_reads = {}
    for stats_file in glob.glob("*.filter_stats.txt"):
        sample = stats_file.replace(".filter_stats.txt", "")
        kept = 0
        with open(stats_file) as fh:
            for line in fh:
                m = re.search(r"(\\d+) sequences kept", line)
                if m:
                    kept = int(m.group(1))
                    break
        total_reads[sample] = kept

    # -------------------------------------------------------
    # 3. Sum mapped reads per sample from count_table
    #    count_table is tab-separated: first col = #OTU ID (ASV),
    #    remaining cols = sample names
    # -------------------------------------------------------
    mapped_reads = defaultdict(int)
    with open("${count_table}") as fh:
        header_cols = fh.readline().rstrip("\\n").split("\\t")[1:]  # sample names
        for line in fh:
            parts = line.rstrip("\\n").split("\\t")
            for i, val in enumerate(parts[1:]):
                try:
                    mapped_reads[header_cols[i]] += int(float(val))
                except ValueError:
                    pass

    # -------------------------------------------------------
    # 4. Write per-sample summary
    # -------------------------------------------------------
    all_samples = sorted(set(list(total_reads.keys()) + list(mapped_reads.keys())))
    with open("mapping_rate_summary.tsv", "w") as fh:
        fh.write("sample\\ttotal_reads_post_filter\\tmapped_reads\\tmapping_rate_pct\\n")
        for sample in all_samples:
            total = total_reads.get(sample, 0)
            mapped = mapped_reads.get(sample, 0)
            rate = (mapped / total * 100) if total > 0 else 0.0
            fh.write(f"{sample}\\t{total}\\t{mapped}\\t{rate:.2f}\\n")

    print(f"Mapping rate summary written for {len(all_samples)} samples.")
    """

    stub:
    """
    touch mapping_rate_summary.tsv
    touch mapping_rate_overall.txt
    """
}
