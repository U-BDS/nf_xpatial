#!/usr/bin/env nextflow

include { GET_VARIABLE_FEATURES     } from '../../../subworkflows/local/get_variable_features'

include { GENERATE_GENE_PAIR_STATS   } from '../../../modules/local/generate_gene_pair_stats'
include { FILTER_GENE_PAIRS          } from '../../../modules/local/filter_gene_pairs'
include { QC_BARNYARD_PLOT           } from '../../../modules/local/qc_barnyard_plot'
include { CONCAT_CSV                 } from '../../../modules/local/concat_csv'
include { QC_HEATMAP_PLOT            } from '../../../modules/local/qc_heatmap_plot'

workflow MARKER_GENE_PAIRS_QC {
    take:
        ch_xenium_data            // channel: annotated xenium data
        marker_gene_list          // file: marker gene list
        skip_gene_list_filtering  // bool: whether to filter gene list used by barnyard plots
        skip_marker_barnyard_plot // bool: whether to skip the barnyard qc plots

    main:
        ch_versions = Channel.empty()

        // If no marker gene list is provided, identify variable genes to use as the gene list
        ch_gene_list = Channel.empty()
        if (!marker_gene_list) {
            GET_VARIABLE_FEATURES (
                ch_xenium_data
            )

            ch_gene_list = GET_VARIABLE_FEATURES.out.vf_xenium_obj
        } else {
            ch_gene_list = Channel.from(marker_gene_list)
        }

        //
        // MODULE: Generate Gini Score and Spearman Correlation for each gene pair
        //
        GENERATE_GENE_PAIR_STATS (
            ch_xenium_data.combine ( ch_gene_list )
        )
        ch_versions = ch_versions.mix ( GENERATE_GENE_PAIR_STATS.out.versions )

        ch_gene_list = GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
        if (!skip_gene_list_filtering) {
            //
            // MODULE: Determine mutually exclusive gene pairs
            //
            FILTER_GENE_PAIRS (
                GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
            )
            ch_versions = ch_versions.mix ( FILTER_GENE_PAIRS.out.versions )

            ch_gene_list = FILTER_GENE_PAIRS.out.filtered_gene_pair_stats
        }

        ch_marker_barnyard = Channel.empty()
        if (!skip_marker_barnyard_plot) {
            //
            // MODULE: Generate Barnyard Plot
            //
            QC_BARNYARD_PLOT (
                ch_xenium_data.join ( ch_gene_list )
            )

            ch_marker_barnyard = QC_BARNYARD_PLOT.out.barnyard_plot
        }

        // //
        // // MODULE: Concatenate CSVs
        // //
        // CONCAT_CSV (
        //    FILTER_GENE_PAIRS.out.gene_pair_stats
        //        .map{
        //            meta, gene_stat_csv -> [gene_stat_csv]
        //        }
        //        .collect()
        //        .map{
        //            [ [ 'id': 'compiled'], it ]
        //        }
        // )

        // //
        // // MODULE: Generate Heatmap Plot
        // //
        // QC_HEATMAP_PLOT (
        //    CONCAT_CSV.out.concat_csv
        // )

    emit:
        gene_pair_stats           = GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
        gene_list                 = ch_gene_list
        marker_barnyard_plot      = ch_marker_barnyard
        versions                  = ch_versions
}
