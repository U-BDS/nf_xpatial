#!/usr/bin/env nextflow

include { GENERATE_GENE_PAIR_STATS   } from '../../../modules/local/generate_gene_pair_stats'
include { DETERMINE_MUTEX_GENE_PAIRS } from '../../../modules/local/determine_mutex_gene_pairs'
include { QC_BARNYARD_PLOT           } from '../../../modules/local/qc_barnyard_plot'

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

        // Cartesian product the xenium data with ALL gene pairings
        ch_xenium_data
            .map{
                meta, xenium_rds, gene_list ->
                    [meta, xenium_rds]
            }
            .combine (
                GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
                    .flatMap{
                        meta, gene_pair_stats ->
                            gene_pair_stats.splitCsv(header: true)
                    }
                    .map{ row -> [row['gene1'], row['gene2']] }
                    .unique()
        )
        .map {
            meta, xenium_rds, gene1, gene2 ->
                [meta, xenium_rds, ['gene1': gene1, 'gene2': gene2]]
        }
        .set { ch_xenium_gene_pairs }

        //
        // MODULE: Generate Barnyard Plot
        //
        QC_BARNYARD_PLOT (
            ch_xenium_gene_pairs
        )

        //
        // MODULE: Generate Heatmap Plot
        //

    emit:
        gene_pair_stats = GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
        exclusive_gene_pair_stats = DETERMINE_MUTEX_GENE_PAIRS.out.mutex_gene_pair_stats
        barnyard_plot = QC_BARNYARD_PLOT.out.barnyard_plot
        heatmap_plot = null
        versions = ch_versions


}