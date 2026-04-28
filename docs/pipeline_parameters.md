

# U-BDS/nf_xpatial pipeline parameters                                                                                                                                                        
                                                                                                                                                                                              
This pipeline is used to automate analysis of Xenium data                                                                                                                                     
                                                                                                                                                                                              
## Input/output Options                                                                                                                                                                       
                                                                                                                                                                                              
Define where the pipeline should find input data and save output data.                                                                                                                        
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `input` | Path to samplesheet.csv <details><summary>Help</summary><small>Path to `samplesheet.csv`</small></details>| `string` |  | True |  |                                               
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |                               
| `marker_gene_list` | The list of marker genes to evaluate and plot | `string` |  |  |  |                                                                                                    
| `gene_whitelist` | A list of features to always keep if you subset by variable features. | `string` |  |  |  |                                                                              
                                                                                                                                                                                              
## General Pipeline Options                                                                                                                                                                   
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `clustering_methods` | This lets you choose the clustering methods to use for your data <details><summary>Help</summary><small>Options are: Seurat,BANKSY,BANKSYSeurat</small></details>|   
`string` | Seurat,BANKSY | True |  |                                                                                                                                                          
| `skip_qc` | Skip all QC Plots | `boolean` |  |  |  |                                                                                                                                        
                                                                                                                                                                                              
## Manual Annotation Options                                                                                                                                                                  
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `skip_man_ann_dim_plot` | Skip generating the Manual Annotation Dim Plot | `boolean` |  |  |  |                                                                                             
                                                                                                                                                                                              
## Filtering Options                                                                                                                                                                          
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `skip_filtering` | Skip filtering alltogether | `boolean` |  |  |  |                                                                                                                        
| `skip_column_removal` | Skip unneeded column removal | `boolean` |  |  |  |                                                                                                                 
| `columns_to_remove` | A comma-delimited list of unnecessary columns to remove from the xenium object | `string` | nCount_BlankCodeword,<br>nFeature_BlankCodeword,<br>nCount_ControlCodeword,<br>nFeature_ControlCodeword,<br>nCount_ControlProbe,<br>nFeature_ControlProbe |  |  |                                                 
| `skip_percentile_filtering` | Skip the percentile based filtering where cells are filtered based on nFeature percentiles instead of absolute nFeature value | `boolean` | True |  |  |      
| `min_percentile` | The minimum percentile threshold for nFeature percentile based filtering | `number` |  |  |  |                                                                           
| `max_percentile` | The maximum percentile threshold for nFeature percentile based filtering | `number` |  |  |  |                                                                           
| `skip_nFeature_filtering` | Skip the nFeature filtering | `boolean` |  |  |  |                                                                                                              
| `min_nFeature` | The minimum nFeature threshold for nFeature filtering | `number` | 10 |  |  |                                                                                              
| `max_nFeature` | The maximum nFeature threshold for nFeature filtering | `number` |  |  |  |                                                                                                
| `skip_nCount_filtering` | Skip the nCount filtering | `boolean` |  |  |  |                                                                                                                  
| `min_nCount` | The minimum nCount threshold for nCount filtering | `number` | 10 |  |  |                                                                                                    
| `max_nCount` | The maximum nCount threshold for nCount filtering | `number` |  |  |  |                                                                                                      
| `skip_cell_area_filtering` | Skip the cell area filtering | `boolean` | True |  |  |                                                                                                        
| `min_cell_area` | The minimum cell area for cell area filtering | `integer` |  |  |  |                                                                                                      
| `max_cell_area` | The maximum cell area for cell area filtering | `integer` |  |  |  |                                                                                                      
| `filter_features` | A list of features to filter from samples | `string` |  |  |  |                                                                                                         
                                                                                                                                                                                              
## Normalization Options                                                                                                                                                                      
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `normalization_method` | Choose the normalization method to use (log, and/or area). Multiple options can be chosen in a comma delimited list | `string` | log |  |  |                       
| `skip_norm_ncount` | Skip generating the nCount plot for normalized data | `boolean` |  |  |  |                                                                                             
| `skip_norm_nfeature` | Skip generating the nFeature plot for normalized data | `boolean` |  |  |  |                                                                                         
                                                                                                                                                                                              
## Seurat Clustering Options                                                                                                                                                                  
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `skip_tsne_plot` | Use this to specify to skip producing a tsne plot | `boolean` | True |  |  |                                                                                             
| `dim_Seurat` | The comma-delimited list of dimensions to use for Seurat clustering | `string` |  |  |  |                                                                                    
| `res_Seurat` | The comma-delimited list of resolutions to use for Seurat clustering | `string` |  |  |  |                                                                                   
                                                                                                                                                                                              
## BANKSY Clustering Options                                                                                                                                                                  
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `skip_banksy_vf_filter` | This will skip subsetting the xenium object by highly variable genes before going into BANKSY | `boolean` |  |  |  |                                              
| `vf_nfeatures` | The number of variable features to keep if filtering by variable feature is enabled | `number` | 2000 |  |  |                                                              
| `lambda_BANKSY` | The comma-delimited list of lambda values to pass to BANKSY | `string` |  |  |  |                                                                                         
| `k_geom_BANKSY` | The comma-delimited list of k_grom values to pass to BANKSY | `string` |  |  |  |                                                                                         
| `nPCs_BANKSY` | The comma-delimited list of nPCs to pass to BANKSY | `string` |  |  |  |                                                                                                    
| `res_BANKSY` | The comma-delimited list of resolutions to pass to BANKSY | `string` |  |  |  |                                                                                              
| `use_agf_BANKSY` | Whether or not to use agf for BANKSY | `boolean` | True |  |  |                                                                                                          
| `max_iter_BANKSY` | The maximum number of iterations to use for BANKSY | `integer` | 60 |  |  |                                                                                             
                                                                                                                                                                                              
## General QC Options                                                                                                                                                                         
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `skip_pre_general_qc` | Skip the General QC steps for pre-filtered data | `boolean` |  |  |  |                                                                                              
| `skip_pre_general_dim_plot` | Skip generating dim plot for pre-filtered data | `boolean` |  |  |  |                                                                                         
| `skip_pre_general_vln_plot` | Skip generating vln plot for pre-filtered data | `boolean` |  |  |  |                                                                                         
| `skip_pre_general_scatter_plot` | Skip generating nFeature/nCount scatter plot for pre-filtered data | `boolean` |  |  |  |                                                                 
| `skip_pre_general_nfeature_plot` | Skip generating nFeature plot for pre-filtered data | `boolean` |  |  |  |                                                                               
| `skip_pre_general_ncount_plot` | Skip generating nCount plot for pre-filtered data | `boolean` |  |  |  |                                                                                   
| `skip_post_general_qc` | Skip the General QC steps for post-filtered data | `boolean` |  |  |  |                                                                                            
| `skip_post_general_dim_plot` | Skip generating dim plot for post-filtered data | `boolean` |  |  |  |                                                                                       
| `skip_post_general_vln_plot` | Skip generating vln plot for post-filtered data | `boolean` |  |  |  |                                                                                       
| `skip_post_general_scatter_plot` | Skip generating nFeature/nCount scatter plot for post-filtered data | `boolean` |  |  |  |                                                               
| `skip_post_general_nfeature_plot` | Skip generating nFeature plot for post-filtered data | `boolean` |  |  |  |                                                                             
| `skip_post_general_ncount_plot` | Skip generating nCount plot for post-filtered data | `boolean` |  |  |  |                                                                                 
                                                                                                                                                                                              
## Marker Gene QC Options                                                                                                                                                                     
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `max_spearman_p` | The maximum spearman p value to determine if two genes are co-expressed | `number` | 0.05 |  |  |                                                                        
| `max_spearman_r` | The maximum spearman r value to determine if two genes are co-expressed | `number` | -0.2 |  |  |                                                                        
| `min_gini` | The minimum gini score to determine if two genes are co-expressed | `number` | 0.7 |  |  |                                                                                     
| `skip_gene_list_filtering` | Determines whether to perform gene list filtering or not | `boolean` |  |  |  |                                                                                
| `skip_gene_pair_stats_filtering` | Skip filtering data by the spearman r, spearman p, and gini index values | `boolean` | True |  |  |                                                      
| `skip_marker_gene_qc` | Skip marker gene QC | `boolean` |  |  |  |                                                                                                                          
| `skip_marker_barnyard_plot` | Skip generating Barnyard Plots for gene pairs | `boolean` |  |  |  |                                                                                          
                                                                                                                                                                                              
## Cell Shape QC                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `skip_pre_cell_shape_qc` | Skip Cell Shape QC steps for pre-filtered datea | `boolean` |  |  |  |                                                                                           
| `skip_pre_cell_segm_prop_plot` | Skip generating the Cell Segmentation Proportion Plot for pre-filtered data | `boolean` |  |  |  |                                                         
| `skip_pre_cell_shape_prop_plot` | Skip generating the Cell Shape Proportion Plot for pre-filtered data | `boolean` |  |  |  |                                                               
| `skip_post_cell_shape_qc` | Skip Cell Shape QC steps for pre-filtered datea | `boolean` |  |  |  |                                                                                          
| `skip_post_cell_segm_prop_plot` | Skip generating the Cell Segmentation Proportion Plot for post-filtered data | `boolean` |  |  |  |                                                       
| `skip_post_cell_shape_prop_plot` | Skip generating the Cell Shape Proportion Plot for post-filtered data | `boolean` |  |  |  |                                                             
                                                                                                                                                                                              
## Cell Area QC Options                                                                                                                                                                       
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `skip_pre_cell_area_qc` | Skip the Cell Area QC steps for pre-filtered data | `boolean` |  |  |  |                                                                                          
| `skip_pre_area_histogram_plot` | Skip generating the Area Histogram Plot for pre-filtered data | `boolean` |  |  |  |                                                                       
| `skip_pre_area_box_plot` | Skip generating the Area Box Plot for pre-filtered data | `boolean` |  |  |  |                                                                                   
| `skip_pre_area_overlap_histogram_plot` | Skip generating the Area Overlapping Histogram for pre-filtered data | `boolean` |  |  |  |                                                        
| `skip_post_cell_area_qc` | Skip the Cell Area QC steps for post-filtered data | `boolean` |  |  |  |                                                                                        
| `skip_post_area_histogram_plot` | Skip generating the Area Histogram Plot for post-filtered data | `boolean` |  |  |  |                                                                     
| `skip_post_area_box_plot` | Skip generating the Area Box Plot for post-filtered data | `boolean` |  |  |  |                                                                                 
| `skip_post_area_overlap_histogram_plot` | Skip generating the Area Box Plot for post-filtered data | `boolean` |  |  |  |                                                                   
                                                                                                                                                                                              
## Cluster QC Options                                                                                                                                                                         
                                                                                                                                                                                              
                                                                                                                                                                                              
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `skip_cluster_umap_plot` | Skip generating the UMAP plot for clustered data | `boolean` |  |  |  |                                                                                          
| `skip_cluster_split_plot` | Skip generating the UMAP + Spatial Dim plot split out by sample | `boolean` |  |  |  |                                                                          
| `skip_cluster_vln_plot` | Skip generating the marker violin plot for clustered data | `boolean` | True |  |  |                                                                              
| `skip_cluster_dot_plot` | Skip generating the marker dot plot for clustered data | `boolean` |  |  |  |                                                                                     
                                                                                                                                                                                              
## Institutional Config Options                                                                                                                                                               
                                                                                                                                                                                              
Parameters used to describe centralised config profiles. These should not be edited.                                                                                                          
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |                                                                                          
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional     
config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with
this parameter.</small></details>| `string` |  |  | True |                                                                                                                                    
| `config_profile_name` | Institutional config name. | `string` |  |  | True |                                                                                                                
| `config_profile_description` | Institutional config description. | `string` |  |  | True |                                                                                                  
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |                                                                                              
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |                                                                                                             
                                                                                                                                                                                              
## Generic Options                                                                                                                                                                            
                                                                                                                                                                                              
Less common options for the pipeline, typically set in a config file.                                                                                                                         
                                                                                                                                                                                              
| Parameter | Description | Type | Default | Required | Hidden |                                                                                                                              
|-----------|-----------|-----------|-----------|-----------|-----------|                                                                                                                     
| `version` | Display version and exit. | `boolean` |  |  | True |                                                                                                                            
| `publish_dir_mode` | Method used to save pipeline results to output directory. (accepted: `symlink`\|`rellink`\|`link`\|`copy`\|`copyNoFollow`\|`move`)                                     
<details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method
should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>
