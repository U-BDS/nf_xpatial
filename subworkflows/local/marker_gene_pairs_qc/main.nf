#!/usr/bin/env nextflow

include { GENERATE_GENE_PAIR_STATS   } from '../../../modules/local/generate_gene_pair_stats'
include { DETERMINE_MUTEX_GENE_PAIRS } from '../../../modules/local/determine_mutex_gene_pairs'

workflow MARKER_GENE_PAIRS_QC {
    take:
        ch_xenium_data // channel: annotated xenium data with associated marker genes

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Generate Gini Score and Spearman Correlation for each gene pair
        //
        GENERATE_GENE_PAIR_STATS (
            ch_xenium_data
        )
        ch_versions = ch_versions.mix(GENERATE_GENE_PAIR_STATS.out.versions)

        //
        // MODULE: Determine mutually exclusive gene pairs
        //
        DETERMINE_MUTEX_GENE_PAIRS (
            GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
        )
        ch_versions = ch_versions.mix(DETERMINE_MUTEX_GENE_PAIRS.out.versions)

        //
        // MODULE: Generate Barnyard Plot
        //


        //
        // MODULE: Generate Heatmap Plot
        //

    emit:
        gene_pair_stats = GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
        exclusive_gene_pair_stats = DETERMINE_MUTEX_GENE_PAIRS.out.mutex_gene_pair_stats
        barnyard_plot = null
        heatmap_plot = null
        versions = ch_versions


}