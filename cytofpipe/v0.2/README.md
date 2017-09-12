Cytofpipe v0.2
==============
_________________


This pipeline was developed to by Lucia Conde at the BLIC, in collaboration with Jake Henry from the Immune Regulation and Tumour Immunotherapy Research Group, for the automatic analysis of flow and mass cytometry data in the UCL cluster __Legion__. The methods underneath Cytofpipe v0.2 are based on available R packages. Currently, Cytofpipe v0.2 can be used to run standard cytometry data analysis for subset identification, and to construct scaffold maps for visualizing complex relationships between samples. 

- Cytofpipe **--clustering** is based mainly on cytofkit (https://bioconductor.org/packages/release/bioc/html/cytofkit.html), which is used for preprocessing/clustering, and openCyto (https://www.bioconductor.org/packages/release/bioc/html/openCyto.html), for basic automated gating. 

- Cytofpipe **--scaffold** is based on scaffold (https://github.com/nolanlab/scaffold) 

<br />

### 1. How to run the cytofpipe pipeline

_________________


#### 1.1. Connect to legion and bring the inputfiles

You will need to connect to legion(apply for an account here: https://wiki.rc.ucl.ac.uk/wiki/Account_Services), and transfer there a folder with the input FCS files, a file with the list of markers (which will be 
the ones used for clustering), and optionally a config file (at the moment only for cytofpipe --clustering).

To connect to legion, you can use Putty if you have Windows (check this UCL link: https://wiki.rc.ucl.ac.uk/wiki/Accessing_RC_Systems) or use SSH from your Mac terminal:

`$ ssh UCL_ID@legion.rc.ucl.ac.uk`

To transfer the files to legion, you can either use SCP from your laptop, for example:

`$ scp -r FILES UCL_ID@legion.rc.ucl.ac.uk:/home/user/Scratch/my_cytof_analysis/.`

or if you have a FTP transfer program (for example cyberduck: http://download.cnet.com/Cyberduck/3000-2160_4-10246246.html) you can also transfer the files from/to legion simply by dragging them from one window to another.

<br />

#### 1.2. Load the modules

Once you are in legion, you will need to load the modules necessary for the pipeline to work.

`$ module load blic-modules`

`$ module load cytofpipe/v0.2`

<br />

#### 1.3. Run the pipeline

Let’s say you have a folder called my_cytof_analyses in your home in legion that contains a directory with the FCS files, the parameters file and perhaps a config file, for example:

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
/home/user/Scratch/my_cytof_analyses/parameters.txt
/home/user/Scratch/my_cytof_analyses/config.txt
```


To run the pipeline in "clustering" mode, just go to the my_cytof_analyses folder and run:

`$ cytofpipe --clustering inputfiles results parameters.txt config.txt`

That will crate a new folder called “results” that will contain all the results of the analysis.

Similarly, to run the pipeline in "scaffold" mode (using file1.fcs as reference and asinh_cofactor = 5):

`$ cytofpipe --scaffold inputfiles file1.fcs results parameters.txt 5`

<br />

#### 1.4. Errors before running

When you submit the job, before it actually runs, there is a script that checks that everything is in order. For example, that the inputfiles folder exists, that there is not a results folder already there (so that nothing is overwritten), that the config.txt file has the appropriate format, etc... Only if everything looks fine, the job will be submitted. Otherwise, an error message will appear that will tell you that there is a problem. For example:

```
------------------------------------------------------------------------------------
USAGE: cytofpipe --clustering <INPUTDIR> <OUTPUTDIR> <MARKERSFILE> [<CONFIGFILE>]
------------------------------------------------------------------------------------
Where CONFIGFILE has the following format:

	[ cytofpipe ] 		        #-- MANDATORY FIELD, IT SHOULD BE THE FIRST LINE OF THE CONFIG FILE

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


Unable to run job: The outputdir <my_output> already exists, please choose a different outputdir.
Exiting.
```

<br />

#### 1.5. Check that the job is running

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

### 2. Running cytofpipe in --clustering mode

_________________


Cytofpipe **--clustering** can be used to analyze data from multiple FCS files.
<br />

First, FCS files will be merged, expression data for selected markers will be transformed, and data will be downsampled according to the user's specifications. Then, clustering will be performed to detect cell types. Finally, the high dimensional flow/mass cytometry data will be visualized into a two-dimensional map with colors representing cell type, and heatmaps to visualize the median expression for each marker in each cell type will be generated.
<br />

Cytofpipe assumes that the data has been properly preprocessed beforehand, i.e., that  normalisation, debarcoding and compensation (if flow) were done properly, and that all the debris, doublets, and live_neg events were removed before analysis. 
However, if manual gating has not been done, automatic gating of live events can be done by setting the "Gating" option in the config file to "YES". But importantly, the default Cytofpipe gating strategy (shown below) is based on the manual gating example provided by Jake. **The gating will not work properly on datasets coming from other sources where the staining panel, instrument, or other factors are different than the ones used by Jake at the Cancer Institute.**

<p align="center">
![](/shared/ucl/depts/cancer/apps/cytofpipe/v0.2/doc/gating.png)
</p>

Automated gating can be performed for more complex assays, however, a proper gating template needs to be created before using the pipeline for gating these dtasets. Any type of dataset can be used with the current pipeline if gating is not needed.

<br />

#### 2.1. Inputfiles

<br />
```
USAGE: cytofpipe --clustering <INPUTDIR> <OUTPUTDIR> <MARKERSFILE> [<CONFIGFILE>]
```
<br />

- **INPUTDIR**: Just a folder with FCS files

- **OUTPUTDIR**: This is just the name for the folder where you want to output the results. It can not be an existing folder.

- **MARKERSFILE**: A text file with the names of the markers, one per line. For example:

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

- **CONFIGFILE**: The config file has to have the following format: 

```
[ cytofpipe ] 		        #-- MANDATORY FIELD, IT SHOULD BE THE FIRST LINE OF THE CONFIG FILE

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

```

The config file is not mandatory. If is not provided, the pipeline will use a default config.txt file, which has GATING = no, TRANSFORM = cytofAsinh, MERGE = fixed (n = 10,000), PHENOGRAPH = yes (other clustering methods = no), TSNE parameters: perplexity =30, theta = 0.5, max_iter = 1000.


<br />

#### 2.2. Outputfiles


<p align="center">
![](/shared/ucl/depts/cancer/apps/cytofpipe/v0.2/doc/clustering.png)
</p>


- **Rphenograph**: Contains the data/plots from the clustering. i.e., the markers average values per cluster, cell percentages in each cluster per FCS file, heatmaps… If more than one clustering method was selected, then there will be several folders, one per clustering method (i.e., one Rphenograph folder, one clusterX folder, etc..)

- **Gating**: If GATING was selected, there will be folder with plots that show the gating applied to each FCS file, a plot that shows the gating scheme, and a folder called “gating_fs_live” that will have the FCS files for the live population.


- **log_R.txt**: This is just the log file from the R script, just so that the user can see what R commands were run when doing the analysis. It will help me figure out what the problem is if the job finishes with an error. 

- **summary.pdf**: This is a PDF with a summary of the analysis. It describes what files, markers and config options were used, and shows the main plots from the analysis (gates, cluster plots, marker level plots, heatmaps…).

- **cytofpipe_analyzedFCS**: The original FCS files (i.e., the original files if gating = YES, or the “gating_fs_live” files if gating = NO), with the clutering and tSNE added information.

- **cytofpipe_tsne_level_plot.pdf**: Marker level plot (shows the expression level of markers on the tSNE data in one plot)

- **cytofpipe_tsne_dimension_reduced_data.csv**: tSNE data (tSNE1 and tSNE2 values per event)

- **cytofpipe_markerFiltered_transformed_merged_exprssion_data.csv**: Expression data (expression of each marker per event)

- **cytofpipe.RData**: The object resulting from the cytofkit analysis. i.e., a file that was saved and that can be used for loading to the cytofkit shiny APP to further exploration of the results.


<br />


### 2. Running cytofpipe in --scaffold mode

_________________


Cytofpipe **--scaffold** can be used to generate scaffold maps to compare cell population networks across different conditions.
<br />


Cytofpipe clusters each FCS sample file independently (currently set up to 200 clusters) using the clara function in R, as implemented in <a href="https://github.com/nolanlab/scaffold">scaffold</a>. A graph is then constructed connecting the nodes fom the manually gated populations and the clusters from the reference FCS file, with edge weights defined as the cosine similarity between the vectors of median marker values of each cluster. Edges of low weight are filtered out and the graph is then laid out (shaped) using a ForceAtlas2 algorithm implemented in <a href="https://github.com/nolanlab/scaffold">scaffold</a>. Graphs are generated for every FCS file, where the position of the landmark nodes stay constant, providing a visual reference that allows the comparison of the different datasets.

Cytofpipe assumes that the FCS data has been properly preprocessed beforehand, i.e., that  normalisation, debarcoding and compensation (if flow) were done properly, and that all the debris, doublets, and live_neg events were removed before analysis. With regards to compensation, please note that the software will try to apply a compensation matrix if one is found in the FCS file. So if the data in the file is already compensated, it will be erroneously compensated again. If your data needs compensation make sure that the FCS file has a spillover matrix embedded and the data is uncompensated.

<br />

#### 2.1. Inputfiles

<br />
```
USAGE: cytofpipe --scaffold <INPUTDIR> <REF_FCS> <OUTPUTDIR> <MARKERSFILE> [<ASINH_COFACTOR>]
```
<br />

- **INPUTDIR**: Folder with FCS files. It should also contain a "gated" subfolder with the manually gated populations (see below).

- **REF_FCS**: The FCS file that will be used for the construction of the reference map together with the landmark (manually gated) populations. It should be one of the FCS files inside INPUTDIR

- **OUTPUTDIR**: Name for the folder where you want to output the results. It can not be an existing folder.

- **MARKERSFILE**: A text file with the names of the markers that will be used in the clustering step, one per line. For example:
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

- **ASINH_COFACTOR**: the cofactor used for the asinh transformation (optional, default = 5). Recommended values are 5 for mass cytometry data, and 150 for fluorescence-based flow cytometry. 

The landmark populations have to be provided as single FCS files (one for each population) that need to be located in a subdirectory called "gated" of the INPUTDIR directory. The program will split the name of the FCS file using "_" as separator and the last field will be used as the population name. For instance if you want an FCS file to define your "B cells" population you have to use the following naming scheme: XXXX_B cells.fcs


<br />

#### 2.2. Outputfiles


<p align="center">
![](/shared/ucl/depts/cancer/apps/cytofpipe/v0.2/doc/scaffold.png)
</p>


- **Clustering**: Contains the files created after the clustering step. For each FCS files two files will be created:
*.clustered.txt: this file contains the marker medians for each cluster
*.clustered.all_events.RData: this file is an RData object which contains all the events in the original FCS file but with an added column that specifies the cluster membership. The data in this file is arcsinh transformed

- **log_R.txt**: This is just the log file from the R script, just so that the user can see what R commands were run when doing the analysis. It will help me figure out what the problem is if the job finishes with an error. 

- **scaffold_map_XXX.pdf**: These are the PDFs with the scaffold maps, one for each input dataset. By default, landmark nodes are coloured in red and population clusters in blue. 

- **XXX.scaffold**: A .scaffold file with the same name of the dataset that you have used as reference. This is a single self-contained bundle that contains everything you need to browse the data. It can be loaded into the original scaffold software for further exploration of the results.


----------------

<p align="right">
Questions? Email me <a href="mailto:l.conde@ucl.ac.uk?">here</a>.
<br>Last modified Aug 2017.
</p>

<br />
