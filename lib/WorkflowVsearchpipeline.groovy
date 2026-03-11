//
// This file holds several functions specific to the workflow/vsearchpipeline.nf in the vsearchpipeline pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowVsearchpipeline {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        // here you can add checks for mandatory parameters

    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        def citation_text = [
                "Tools used in the workflow included:",
                "FastQC (Andrews 2010),",
                "MultiQC (Ewels et al. 2016),",
                "VSEARCH (Rognes et al. 2016),",
                "seqtk (Li 2023),",
                "MAFFT (Katoh & Standley 2013),",
                "FastTree (Price et al. 2010),",
                "SILVA (Quast et al. 2013),",
                "DADA2 (Callahan et al. 2016),",
                "phyloseq (McMurdie & Holmes 2013),",
                params.run_decontam ? "decontam (Davis et al. 2018)," : "",
                params.skip_picrust ? "" : "PICRUSt2 (Douglas et al. 2020),",
                "."
            ].findAll { it }.join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.</li>",
                "<li>Ewels P, Magnusson M, Lundin S, Käller M. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19):3047-3048. doi: 10.1093/bioinformatics/btw354</li>",
                "<li>Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4:e2584. doi: 10.7717/peerj.2584</li>",
                "<li>Li H. (2023) seqtk: Toolkit for processing sequences in FASTA/Q formats. URL: https://github.com/lh3/seqtk</li>",
                "<li>Katoh K, Standley DM. (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol, 30(4):772-80. doi: 10.1093/molbev/mst010</li>",
                "<li>Price MN, Dehal PS, Arkin AP. (2010) FastTree 2 – approximately maximum-likelihood trees for large alignments. PLoS ONE, 5(3):e9490. doi: 10.1371/journal.pone.0009490</li>",
                "<li>Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO. (2013) The SILVA ribosomal RNA gene database project. Nucleic Acids Res, 41(D1):D590-D596. doi: 10.1093/nar/gks1219</li>",
                "<li>Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. (2016) DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods, 13(7):581-583. doi: 10.1038/nmeth.3869</li>",
                "<li>McMurdie PJ, Holmes S. (2013) phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE, 8(4):e61217. doi: 10.1371/journal.pone.0061217</li>",
                params.run_decontam ? "<li>Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. (2018) Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome, 6(1):226. doi: 10.1186/s40168-018-0605-2</li>" : "",
                params.skip_picrust ? "" : "<li>Douglas GM et al. (2020) PICRUSt2 for prediction of metagenome functions. Nat Biotechnol, 38:685-688. doi: 10.1038/s41587-020-0548-6</li>",
            ].findAll { it }.join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

}
