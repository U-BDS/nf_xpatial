#!/usr/bin/env nextflow

include { IDENTIFY_VARIABLE_GENES    } from '../../../modules/local/identify_variable_genes'
include { COMPILE_LISTS              } from '../../../modules/local/compile_lists'
include { GENERATE_GENE_PAIR_STATS   } from '../../../modules/local/generate_gene_pair_stats'
include { DETERMINE_MUTEX_GENE_PAIRS } from '../../../modules/local/determine_mutex_gene_pairs'
include { QC_BARNYARD_PLOT           } from '../../../modules/local/qc_barnyard_plot'
include { CONCAT_CSV                 } from '../../../modules/local/concat_csv'
include { QC_HEATMAP_PLOT            } from '../../../modules/local/qc_heatmap_plot'

workflow MARKER_GENE_PAIRS_QC {
    take:
        ch_xenium_data   // channel: annotated xenium data
        marker_gene_list // file: marker gene list

    main:
        ch_versions = Channel.empty()

        // If no marker gene list is provided, identify variable genes to use as the gene list
        ch_gene_list = Channel.empty()
        if (!marker_gene_list) {
            IDENTIFY_VARIABLE_GENES (
                ch_xenium_data
            )

            COMPILE_LISTS (
                IDENTIFY_VARIABLE_GENES.out.variable_gene_list
                    .map { meta, gene_list -> [gene_list] }
                    .collect()
                    .map { 
                        gene_list -> 
                            def new_meta = ['id': 'genes']
                        [new_meta, gene_list]
                    }
            )

            ch_gene_list = COMPILE_LISTS.out.compiled_list.map { meta, gene_list -> [gene_list] }
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

        //
        // MODULE: Determine mutually exclusive gene pairs
        //
        DETERMINE_MUTEX_GENE_PAIRS (
           GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
        )
        ch_versions = ch_versions.mix ( DETERMINE_MUTEX_GENE_PAIRS.out.versions )

        //
        // MODULE: Generate Barnyard Plot
        //

        QC_BARNYARD_PLOT (
           ch_xenium_data.join ( GENERATE_GENE_PAIR_STATS.out.gene_pair_stats )
        )

        //
        // MODULE: Concatenate CSVs
        //
        CONCAT_CSV (
           GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
               .map{
                   meta, gene_stat_csv -> [gene_stat_csv]
               }
               .collect()
               .map{
                   [ [ 'id': 'compiled'], it ]
               }
        )

        //
        // MODULE: Generate Heatmap Plot
        //
        QC_HEATMAP_PLOT (
           CONCAT_CSV.out.concat_csv
        )

    emit:
        gene_pair_stats           = GENERATE_GENE_PAIR_STATS.out.gene_pair_stats
        exclusive_gene_pair_stats = DETERMINE_MUTEX_GENE_PAIRS.out.mutex_gene_pair_stats
        barnyard_plot             = QC_BARNYARD_PLOT.out.barnyard_plot
        heatmap_plot              = QC_HEATMAP_PLOT.out.heatmap_plot
        gene_list                 = ch_gene_list
        versions                  = ch_versions


}
