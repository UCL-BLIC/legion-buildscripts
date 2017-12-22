---
output:
  html_document: default
  pdf_document: default
  word_document: default
---


<center> 
# Cytofpipe v1.0 
_________________
</center>

This pipeline was developed to by Lucia Conde at the BLIC - UCL Cancer Institute, in collaboration with Jake Henry from the Immune Regulation and Tumour Immunotherapy Research Group, for the automatic analysis of flow and 
mass cytometry data in the UCL cluster __Legion__. Currently, Cytofpipe v1.0 can be used to run standard cytometry data analysis for subset identification, for comparison of groups of samples, and to construct scaffold maps 
for visualizing complex relationships between samples. The methods underneath Cytofpipe v1.0 are based on publicly available R packages for flow/cytof data analysis

- Cytofpipe **--clustering** is based mainly on cytofkit (https://bioconductor.org/packages/release/bioc/html/cytofkit.html), which is used for preprocessing/clustering, and openCyto (https://www.bioconductor.org/packages/release/bioc/html/openCyto.html), for basic automated gating. 

- Cytofpipe **--scaffold** is based on scaffold (https://github.com/nolanlab/scaffold) 

- Cytofpipe **--citrus** is based on citrus (https://github.com/nolanlab/citrus) 


<br />

<div align=center>
# How to run the cytofpipe pipeline
</div>

<br>

##  {.tabset}


### 1. Connect to legion, bring inputfiles

<br />

You will need to connect to legion (apply for an account here: https://wiki.rc.ucl.ac.uk/wiki/Account_Services), and transfer there a folder with the input FCS files, a file with the list of markers (which will be 
the ones used for clustering), a file listing which samples belong to each condition/group (for citrus), and optionally a config file (at the moment only for cytofpipe --clustering).

To connect to legion, you can use Putty if you have Windows (check this UCL link: https://wiki.rc.ucl.ac.uk/wiki/Accessing_RC_Systems) or use SSH from your Mac terminal:

`$ ssh UCL_ID@legion.rc.ucl.ac.uk`

To transfer the files to legion, you can either use SCP from your laptop, for example:

`$ scp -r FILES UCL_ID@legion.rc.ucl.ac.uk:/home/user/Scratch/my_cytof_analysis/.`

or if you have a FTP transfer program (for example cyberduck: http://download.cnet.com/Cyberduck/3000-2160_4-10246246.html or WinSCP: https://winscp.net/eng/download.php) you can also transfer the files from/to legion simply 
by dragging them from one window to 
another.

<br />

### 2. Load modules

<br />

Once you are in legion, you will need to load the modules necessary for the pipeline to work.

`$ module load blic-modules`

`$ module load cytofpipe/v1.0`

<br />

### 3. Run pipeline

<br />

Let’s say you have a folder called my_cytof_analyses in your home in legion that contains a directory with the FCS files, a file that contains the markers that you want to use for clustering, a file listing which samples belong to 
each condition, and perhaps a config file, for example:

```
/home/user/Scratch/my_cytof_analyses/
/home/user/Scratch/my_cytof_analyses/inputfiles/
/home/user/Scratch/my_cytof_analyses/inputfiles/file1.fcs
/home/user/Scratch/my_cytof_analyses/inputfiles/file2.fcs
/home/user/Scratch/my_cytof_analyses/inputfiles/file3.fcs
/home/user/Scratch/my_cytof_analyses/inputfiles/gated/file1_cellType1.fcs
/home/user/Scratch/my_cytof_analyses/inputfiles/gated/file1_cellType2.fcs
/home/user/Scratch/my_cytof_analyses/inputfiles/gated/file1_cellType3.fcs
/home/user/Scratch/my_cytof_analyses/inputfiles/gated/file1_cellType4.fcs
/home/user/Scratch/my_cytof_analyses/markers.txt
/home/user/Scratch/my_cytof_analyses/conditions.txt
/home/user/Scratch/my_cytof_analyses/config.txt
```

To run the pipeline in "clustering" mode with default parameters, just go to the my_cytof_analyses folder and run:

`$ cytofpipe --clustering -i inputfiles -o results -m markers.txt`

That will crate a new folder called “results” that will contain all the results of the analysis.

Similarly, to run the pipeline in "scaffold" mode with default parameters:

`$ cytofpipe --scaffold -i inputfiles --ref file1.fcs -o results -m markers.txt`

Finally, to run the pipeline in "citrus" mode with default parameters:

`$ cytofpipe --citrus -i inputfiles --cond conditions.txt -o results -m markers.txt`

<br />

### 4. Errors before running

<br />

When you submit the job, before it actually runs, there is a script that checks that everything is in order. For example, that the inputfiles folder exists, that there is not a results folder already there (so that nothing is overwritten), that if there is a config.txt file, it has the appropriate format, etc... Only if everything looks fine, the job will be submitted. Otherwise, an error message will appear that will tell you that there is a problem. For example:

```
------------------------------------------------------------------
		 ** Cytofpipe v1.0 **
			--clustering
------------------------------------------------------------------
Usage: cytofpipe --clustering -i DIR -o DIR -m FILE [options]

Required:
	-i DIR		Input directory with the FCS files
	-o DIR		Output directory where results will be generated
	-m FILE		File with markers that will be selected for clustering

Options:
	--config FILE       Configuration file to customize the analysis (see below) 
	--flow | --cytof      Flow cytometry data (transformation = autoLgcl)
	                    or Cytof data (transformation = cytofAsinh)
	--all | --downsample NUM    Use all events in the analysis or downsample
	                    each FCS file to the specified number of events
	                    (with no replacement for sample with events < NUM)
	--displayAll       	Display all markers in output files

Unable to run job: Please check that you are providing a inputdir (-i), outputdir (-o) and markersfile (-m)
Exiting.

```

<br />

### 5. Check job is running

<br />

If there were no errors found, the job will be submitted to the queue through a qsub system. To check that the job is queued or running, use qstat:

`$ qstat`

```
job-ID  prior   name       user         state submit/start at     queue                     slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
2739095 3.50000 cytof-raw_ regmond      r     04/03/2017 10:52:31 Yorick@node-z00a-011          1        
2739177 0.00000 cytof-gate regmond      qw    04/03/2017 10:59:51                               1        
```

In the above example I have one job (with ID 2739095) that is already running (state = r), and a second job (with ID 2739177) that is in queue (state = q).

If you submit a job, and later on it does not show when you do qstat, that means that it finished.  You should then be able to see a new folder that has the results of the analysis.

<br />

##  {.tabset}

<div align=center>
# Cytofpipe v1.0 commands
</div>

<br>


### --clustering 

<br />

```
Usage: cytofpipe --clustering -i DIR -o DIR -m FILE [options]
```

<br />


Cytofpipe **--clustering** can be used to analyze data from multiple FCS files.
<br />

First, FCS files will be merged, expression data for selected markers will be transformed, and data will be downsampled according to the user's specifications. Then, clustering will be performed to detect cell types. Finally, the high dimensional flow/mass cytometry data will be visualized into a two-dimensional map with colors representing cell type, and heatmaps to visualize the median expression for each marker in each cell type will be generated.

- *Note 1*: The markers uploaded by the user should be the ones provided in the “Description” filed of the FCS. This will usually be a longer format (141Pr_CD38) in cytof data and a shorter format (CD38) in Flow data. However, the shorter version can be used when uploading cytod data.
- *Note 2*: Markers uploaded by the user are used for clustering and dimensional reduction, and by default they are the only ones displayed in the results (heatmaps, etc..). Using the **--displayAll** option will override this and all the markers (with exception of Time, Event, viability and FSC/SSC channels) will be included in the output plots and files. All markers are transformed to the user's specifications with exception of FSC/SSC that are linearly transformed. 
- *Note 3*: Cytofpipe v1.0 runs cytofkit_1_10_0
<br />

Cytofpipe assumes that the data has been properly preprocessed beforehand, i.e., that  normalisation, debarcoding and compensation (if flow) were done properly, and that all the debris, doublets, and live_neg events were removed before analysis. 
However, if manual gating has not been done, automatic gating of live events can be done by setting the "Gating" option in the config file to "YES". But importantly, the default Cytofpipe gating strategy (shown below) is based on the manual gating example provided by Jake. **The gating will not work properly on datasets coming from other sources where the staining panel, instrument, or other factors are different than the ones used by Jake at the Cancer Institute.**

<p align="center">
Root ---> BeadNeg ---> Cells ---> Live
</p>

Automated gating can be performed for more complex assays, however, a proper gating template needs to be created before using the pipeline for gating these dtasets. Any type of dataset can be used with the current pipeline if gating is not needed.

<br />


#### __Command arguments__

**Mandatory arguments**

<ul>
<li>**-i DIR**: A folder with FCS files</li>

<li>**-o DIR**: The name for the folder where you want to output the results. It can not be an existing folder.</li>

<li>**-m FILE**: A text file with the names of the markers, one per line. For example:

```
CD3
CD4
CD8
FOXP3
TCRgd
CD11b
CD56
HLA-DR
```
</li>
</ul>


**Optional arguments**


<ul>
<li> **--config FILE**: The config file is not mandatory. If is not provided, the pipeline will use a default config.txt file, which has GATING = no, TRANSFORM = cytofAsinh, MERGE = ceil (n = 10,000), PHENOGRAPH = yes (other clustering methods = no), DISAPLAY_ALL = no, TSNE parameters: perplexity = 30, theta = 0.5, max_iter = 1000. If provided, it has to have the following format: 

```
[ clustering ] 		        #-- MANDATORY FIELD, IT SHOULD BE THE FIRST LINE OF THE CONFIG FILE

GATING = yes|no 		    #-- MANDATORY FIELD:

TRANSFORM = autoLgcl, cytofAsinh, logicle, arcsinh or none  #-- MANDATORY FIELD
MERGE = ceil, all, min, or fixed			 	            #-- MANDATORY FIELD
DOWNSAMPLE = number between 500 and 100000 			        #-- MANDATORY FIELD if MERGE = fixed/ceil

#- DimRed method (tSNE) parameters:
PERPLEXITY = 30
THETA = 0.5
MAX_ITER = 1000

#- Clustering methods:
PHENOGRAPH = yes|no
CLUSTERX = yes|no
DENSVM = yes|no
FLOWSOM = yes|no
FLOWSOM_K = number between 2 and 50 				  #-- MANDATORY FIELD if FLOWSOM = YES:

#- Additional visualization methods:
PCA = yes|no
ISOMAP = yes|no

#- Other:
DISPLAY_ALL = yes|no

```

</li>

<li>**--flow**, **--cytof**: Shorcut to let cytofpipe know if the user is analyzing flow or cytof data, without having to provide a config file. If **--flow** is selected, the autoLgcl transformation will be used. If **--cytof** is selected, the cytofAsinh transformation will be used. **--flow** and **--cytof** cannot be used at the same time, and they will override the TRANSFORMATION option of the config file if a cofig file is supplied too. </li>

<li>**--all**, **--downsample NUM**: Shorcut to let cytofpipe know if we want to use all the events/downsample the data, without having to provide a config file. **--all** and **--downsample NUM** cannot be used at the same time, and they will override the DOWNSAMPLE option of the config file if a cofig file is supplied too. </li>

<li>**--displayAll**: Shorcut to let cytofpipe know if we want to display all the markers in the output files and plots, without having to provide a config file. **--displayAll** will override the DISPLAY_ALL option of the config file if a cofig file is supplied too. Please note that Time, Event, viability and FSC/SSC channels will not be displayed even if the **--displayAll** option is selected. Please contact me if you want to change this.</li>

</ul>


<br />

#### Outputfiles

<ul>
<li>
**Rphenograph**: Contains the data/plots from the clustering. i.e., the markers average values per cluster, cell percentages in each cluster per FCS file, heatmaps… If more than one clustering method was selected, then there will be several folders, one per clustering method (i.e., one Rphenograph folder, one clusterX folder, etc..)</li>

<li>**Gating**: If GATING was selected, there will be folder with plots that show the gating applied to each FCS file, a plot that shows the gating scheme, and a folder called “gating_fs_live” that will have the FCS files for the live population.</li>

<li>**cytofpipe_analyzedFCS**: The original FCS files (i.e., the original files if gating = YES, or the “gating_fs_live” files if gating = NO), with the clutering and tSNE added information.</li>

<li>**cytofpipe_tsne_level_plot.pdf**: Marker level plot (shows the expression level of markers on the tSNE data in one plot)</li>

<li>**cytofpipe_tsne_dimension_reduced_data.csv**: tSNE data (tSNE1 and tSNE2 values per event)</li>

<li>**cytofpipe_markerFiltered_transformed_merged_exprssion_data.csv**: Expression data (expression of each marker per event)</li>

<li>**summary_clustering.pdf**: This is a PDF with a summary of the analysis. It describes what files, markers and config options were used, and shows the main plots from the analysis (gates, cluster plots, marker level plots, heatmaps…).</li>

<li>**cytofpipe.RData**: The object resulting from the cytofkit analysis. i.e., a file that was saved and that can be used for loading to the cytofkit shiny APP to further exploration of the results.</li>

<li>**log_R.txt**: This is just the log file from the R script, just so that the user can see what R commands were run when doing the analysis. It will help me figure out what the problem is if the job finishes with an error. </li>
</ul>


<br />


### --scaffold

<br />
```
Usage: cytofpipe --scaffold -i DIR -o DIR -m FILE --ref FILE [options]
```
<br />

Cytofpipe **--scaffold** can be used to generate scaffold maps to compare cell population networks across different conditions.
<br />


Cytofpipe clusters each FCS sample file independently (currently set up to 200 clusters) using the clara function in R, as implemented in <a href="https://github.com/nolanlab/scaffold">scaffold</a>. A graph is then 
constructed connecting the nodes fom the manually gated populations and the clusters from the reference FCS file, with edge weights defined as the cosine similarity between the vectors of median marker values of each cluster. Edges of 
low weight are filtered out and the graph is then laid out (shaped) using a ForceAtlas2 algorithm implemented in <a href="https://github.com/nolanlab/scaffold">scaffold</a>. Graphs are generated for every FCS file, where the 
position of the landmark nodes stay constant, providing a visual reference that allows the comparison of the different datasets.

- *Note 1*: Because the clustering is very computationally intensive, by default cytofpipe downsamples the original FCS files to 10,000 events (with NO replacement if the total number of cell in the file is less than 10,000), and then **all the clustering and construction of maps are done on these downsampled files** to be able to run the jobs in a timely fashion. This can be changed with the **--all** and **--downsample NUM** arguments.
- *Note 2*: Cytofpipe v1.0 is running scaffold_0.1
 
Cytofpipe assumes that the FCS data has been properly preprocessed beforehand, i.e., that  normalisation, debarcoding and compensation (if flow) were done properly, and that all the debris, doublets, and live_negevents were removed before analysis. With regards to compensation, please note that the software will try to apply a compensation matrix if one is found in the FCS file. So if the data in the file is already compensated, it will be erroneously compensated again. If your data needs compensation make sure that the FCS file has a spillover matrix embedded and the data is uncompensated.

<br />


#### __Command arguments__

**Mandatory arguments**

<ul>
<li>**-i DIR**: Folder with FCS files. It should also contain a "gated" subfolder with the manually gated populations (see below).</li>

<li>**--ref FILE**: The FCS file that will be used for the construction of the reference map together with the landmark (manually gated) populations. It should be one of the FCS files inside INPUTDIR</li>

<li>**-o DIR**: Name for the folder where you want to output the results. It can not be an existing folder.</li>

<li>**-m FILE**: A text file with the names of the markers that will be used in the clustering step, one per line. For example:
```
CD3
CD4
CD8
FOXP3
TCRgd
CD11b
CD56
HLA-DR
```
</li>
</ul>

The landmark populations have to be provided as single FCS files (one for each population) that need to be located in a subdirectory called "gated" of the INPUTDIR directory. The program will split the name of the FCS file using "_" as separator and the last field will be used as the population name. For instance if you want an FCS file to define your "B cells" population you have to use the following naming scheme: XXXX_B cells.fcs



**Optional arguments**

<ul>

<li>**--flow**, **--cytof**: Shorcut to let cytofpipe know if the user is analyzing flow or cytof data, without having to provide a config file. If **--flow** is selected, arcsinh transformation will be used with asinh_cofactor = 150. If **--cytof** is selected, arcsinh transformation will be used with asinh_cofactor = 5. **--flow** and **--cytof** cannot be used at the same time </li>

<li>**--all**, **--downsample NUM**: Shorcut to let cytofpipe know if we want to use all the events/downsample the data, without having to provide a config file. **--all** and **--downsample NUM** cannot be used at the same time. </li>
</ul>


<br />


#### Outputfiles

<ul>

<li>**Clustering**: Contains the files created after the clustering step. For each FCS files two files will be created:
<ul><li>*.clustered.txt: this file contains the marker medians for each cluster</li>
<li>*.clustered.all_events.RData: this file is an RData object which contains all the events in the original FCS file but with an added column that specifies the cluster membership. The data in this file is arcsinh transformed</li></ul></li>

<li>**downsampled_X**: If downsampling, this folder will contain the FCS files created after downsampling the original FCS files to the  selected number of events. **If downsampling, these would the FCS files that are actually analised**.</li>

<li>**scaffold_map_XXX.pdf**: These are the PDFs with the scaffold maps, one for each input dataset. By default, landmark nodes are coloured in red and population clusters in blue. </li>

<li> **XXX.scaffold**: A .scaffold file with the same name of the dataset that you have used as reference. This is a single self-contained bundle that contains everything you need to browse the data. It can be loaded into the original scaffold software for further exploration of the results.</li>

<li>**summary_scaffold.pdf**: This is a PDF with a summary of the analysis. It describes what files, markers and config options were used, and shows the scaffold maps in reduced size.</li>

<li>**log_R.txt**: Log file from the R 
script, just so that the user can see what R commands were run when doing the analysis. It will help me figure out what the problem is if the job finishes with an error. </li>

</ul>

<br />


### --citrus

<br />
```
Usage: cytofpipe --citrus -i DIR -o DIR -m FILE --cond FILE [options]
```

<br />

Cytofpipe **--citrus** can be used to identify cell populations associated with an experimental or clinical endpoint in flow and cytof data. 
<br />

After indicating which samples belong to each condition, expression data are transformed (using the arcsin hyperbolic transform) and scaled, and cell populations are indentified in every sample using hierarchical clustering based on 
selected markers. Sample features are calculated (either abundance of cell populations in each sample [default], or median expression of specific markers in a cluster), and these are used 
to identify the features (cell populations or markers) that are likely to be predictive (pamr, glmnet) or correlated (sam) with the experimental/clinical 
endpoint. Briefly, predictive models identify the fewest number of markers needed to predict the experimental endpoint. Correlative models detect all the markers that are correlated with the experimental endpoint but that are not necessarily accurate predictors of an experimental outcome. Among the predictive models, 'pamr' can be used with 2 or more groups, whereas 'glmnet' can only deal with 2 groups or continuous endpoint measures. 


- *Note 1*: Citrus requires 8 or more samples in each experimental group. Running Citrus with fewer than 8 samples per group will likely produce spurious results.
- *Note 2*: To ensure that each sample is equally represented in the clustering, by default Citrus selects an equal number of events from each sample (10,000) that are combined and clustered together. This can be changed with the **--all** and **--downsample NUM** arguments.
- *Note 3*: By default, cytofpipe --citrus runs in "abundances" mode. You can change to "medians" mode using the **--medians FILE** parameter by supplying a list of amrkers used for median 
calculation. *Markers that were selected for clustering should not be selected again for statistics (for example, you might want to use surface markers for clustering and phospho-specific 
signal for analysis).*
- *Note 4*: By default, clusters (cell populations) of size lower than 5% of the total number of clustered events will be ignored. If you wish to change the minimum cluster size threshold 
please contact me. 
- *Note 5*: Parameters must be measured on the same channels in each file, the same parameters must be measured in all FCS files (no extras or missing parameters in any FCS file) and measured 
parameters and channels must appear in the same order in each FCS file.
- *Note 6*: Due to the stochastic nature of the citrus algorithm, is recommended to run (repeat) the analysis at least 3 times to make sure the results are not false positives.
- *Note 7*: Cytofpipe v1.0 runs citrus_0.08

Cytofpipe assumes that the data have been properly preprocessed beforehand, i.e., that  normalisation, debarcoding and compensation (if flow) were done properly, and that all the debris, doublets, and live_neg events were removed before analysis. For flow cytometry data, raw FCS file data must be compensated. Compensation matrices stored in FCS files will not be applied when running Cytofpipe in --citrus mode.

<br />


#### __Command arguments__

**Mandatory arguments**

<ul>
<li>**-i DIR**: Folder with FCS files.</li>

<li>**--cond FILE**: A text file that tells Citrus which samples belong to each condition/group, one sample per line. For example:

```
Patient1_FL2.fcs	Case
Patient2_FL2.fcs	Case
Patient3_FL2.fcs	Case
Patient4_FL2.fcs	Case
Patient5_Ref.fcs	Control
Patient6_Ref.fcs	Control
Patient7_Ref.fcs	Control
Patient8_Ref.fcs	Control
```

<li>**-o DIR**: Name for the folder where you want to output the results. It can not be an existing folder.</li>

<li>**-m FILE**: A text file with the names of the markers that will be used in the clustering step, one per line. For example:
```
CD3
CD4
CD8
FOXP3
TCRgd
CD11b
CD56
HLA-DR
```
</li>
</ul>


**Optional arguments**

<ul>

<li>**-medians FILE**: A text file with the names of the markers that will be used as features, one per line. For example:
```
p-NFκB
p-S6
p-PI3K
p-STAT5
```
</li>

<li>**--flow**, **--cytof**: Shorcut to let cytofpipe know if the user is analyzing flow or cytof data, without having to provide a config file. If **--flow** is selected, arcsin hyperbolic transformation will be used with asinh_cofactor = 150. If **--cytof** is selected, arcsin hyperbolic transformation will be used with asinh_cofactor = 5. **--flow** and **--cytof** cannot be used at the same time </li>

<li>**--all**, **--downsample NUM**: Shorcut to let cytofpipe know if we want to use all the events/downsample the data, without having to provide a config file. **--all** and **--downsample NUM** cannot be used at the same time. </li>
</ul>

<br />

#### Outputfiles


<ul>
<li>**MarkerPlots.pdf**: Plots of the clustering hierarchy, one marker per pdf. Helpful for determining the phenotype of identified clusters. </li>

<li>**MarkerPlotsAll.pdf**: Same content as MarkerPlots.pdf in a single PDF instead of many.</li>

<li>**sam_/pamr_/glmnet_results**: A directory containing results from model analysis. There will be one result for each model used to analyze the data. All models operate on the same clustering result. Within a mode result directory you will find:</li>

<ul>
<li>**clusters-{threshold}.pdf**: Shows phenotype of stratifying clusters at specified significance threshold.</li>
<li>**features-{threshold}.pdf**: Plot showing values stratifying features at specified significance threshold.</li>
<li>**features-{threshold}.csv**: Raw values of stratifying features at specified significance threshold.</li>
<li>**featurePlots-{thresholds}.pdf**: Shows the the location and relatedness of identified stratifying clusters in the clustering hierarchy. There will be a separate plot for each tested feature.</li>
<li>**ModelErrorRate.pdf**: Used to determine the accuracy of the constructed model. For predictive models only (PAMR and GLMNET). </li>
</ul>

<li>**exportedClusters**: Folder that contains data (FCS files) exported from each cluster, useful for further analysis on the clusters of interest.</li>

<li>**summary_citrus.pdf**: This is a PDF with a summary of the analysis. It describes what files, markers and config options were used, as well as the MarkerPlotsAll plot.</li>

<li>**cytofpipe_citrusClustering.rData**: Saved version of the clustered data that was computed during the citrus analysis, which can be loaded in R for further exploration of the results.</li>

<li>**log_R.txt**: This is just the log file from the R script, just so that the user can see what R commands were run when doing the analysis. It will help me figure out what the problem is if the job finishes with an error.</li>
</ul>

<br />


##  {.tabset}

<div align=center>
# Version changes
_________________
</div>

<br>


<ul>
<li><b>v1.0</b>
<ul>
<li>Changed command line usage</li>
<li>Added --flow, --cytof, --all, --downsample, --displayAll options</li>
<li>By default, outputfiles display only clustering Markers (use --displayAll to override)
<li>Marker level plots per sample are shown on the same scale on the x- and y-axis</li>
<li>Fixed bug with mergeMethod=fixed from previous cytofkit version (https://github.com/JinmiaoChenLab/cytofkit/issues/12) </li>
<li>More --citrus functionalities exposed ('medians' mode, clusters exported as new FCS files)  </li>
<li>Versions: cytofkit 1.10.0, scaffold 1.0, citrus 0.08</li>
</ul>
<li><b>v0.3</b>
<ul>
<li>Added --citrus mode</li>
<li>Versions: cytofkit 1.8.4, scaffold 1.0, citrus 0.08</li>
</ul>
<li><b>v0.2.1</b>
<ul><li>Output marker level plots per FCS fiyle</li></ul>
<li><b>v0.2</b>
<ul>
<li>Added --scaffold mode</li>
<li>Major changes to adapt cytofpipe to updated cytofkit 1.8.3</li>
</ul>
</ul>

<br />


##  {.tabset}

<div align=center>
# Questions?
_________________
</div>

<br>


Email me <a href="mailto:l.conde@ucl.ac.uk?">here</a>.
<br>Last modified Dec 2017.

<br />
