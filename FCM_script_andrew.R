# Flow Cytometry (FCM) Code
# Andrew Black
# Last updated 08-18-2021


# Watch this YouTube video to see what we are about to do:
	# https://www.youtube.com/watch?v=2INqQNMNaV0


# Run the code below to install packages
	#citation https://github.com/hally166/Cytometry-R-scripts/blob/master/Startingscript.R

	# Installing the "basic" flow cytometry packages. They all live on a repository called Bioconductor.
	install.packages("BiocManager")

	#You can either do library(BiocManager) and install("flowCore) or what I have done below.
	BiocManager::install("flowCore") #interpret .fcs files
	BiocManager::install("flowViz") #basic visualization
	BiocManager::install("ggcyto") #advanced visualization using the ggPlot nomenclature
	BiocManager::install("openCyto") #Used to link various analysis methodologies
	BiocManager::install("flowWorkspace") #used to build analysis templates
	BiocManager::install("CytoML") #imports FlowJo and DiVA workspaces

	#These packages largely require each other to work (except flowCore which is the "base package) 
	#so will often load each other without my help.  For simplicity I have loaded them all.

	#You will need to "clean" your data.  flowAI and flowCut are my recommendations.  
	#flowClean is the original, but is succeeded by flowCut
	BiocManager::install("flowAI")

	if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
	BiocManager::install("flowAI") # install failed, so doing it manually
	
	BiocManager::install("flowClean")

	#flowCut is not available on bioconductor and needs to be loaded straight from GitHub. To do this you need the package devtools.
	install.packages("devtools")
	devtools::install_github("jmeskas/flowCut")

	#An interesting project is CytoExploreR that tries to blend the power of R with the ease of use of a mouse.
	devtools::install_github("DillonHammill/CytoExploreRData")
	devtools::install_github("DillonHammill/CytoExploreR",build_vignettes = TRUE)


# Load data
	# Load a single fcs file
		myfile <- "C:/Users/Andrew/Desktop/Jade_IDR_Biofilter/A01 20210815_01_SYBR.FCS"
		fcsfile <- flowCore::read.FCS(myfile)
		library(flowCore)
		fcsfile1 <- read.FCS(myfile)
		fcsfile

	# Load many fcs files into a flow set
		setwd("C:/Users/Andrew/Desktop/Jade_IDR_Biofilter")
		myfiles <- list.files(pattern = "\\.fcs$") # load all files in this folder that are .fcs file type
		fs <- flowCore::read.flowSet(myfiles, path="C:/Users/Andrew/Desktop/Jade_IDR_Biofilter/")
		fs
		fs[[1]]

# Clean and visualize data
	# Citation for basic code to follow, most notes added are my own: https://github.com/hally166/Cytometry-R-scripts/blob/master/compensate_transform_clean.R
	# Citation for video tutorial: https://www.youtube.com/watch?v=WWa7dwwiLvI
	
	#Load the packages
		library(BiocManager)
		library(flowViz)
		library(openCyto)
		library(devtools)
		library(flowWorkspace)
		library(CytoML)
		library(flowCore)
		library(flowAI)
		library(ggcyto)

		#How to get help
		??flowCore

		#Load a single file
		myfile <- "C:/Users/Andrew/Desktop/Jade_IDR_Biofilter/A01 20210815_01_SYBR.FCS"
		fcsfile <- read.FCS(myfile)
		fcsfile
		names(fcsfile) # list of parameters
		exprs(fcsfile) # Object of class matrix containing the measured intensities. Rows correspond to cells, columns to the different measurement channels
		each_col(fcsfile, median) # show the median of each column
		keyword(fcsfile) # show the metadata for the file

		#Compensation
		spillover(fcsfile) # use this to see which keyword to use
		fcsfile_comp <-compensate(fcsfile, spillover(fcsfile)$'$SPILLOVER') # this was changed, from SPILL to '$SPILLOVER' keyword. . . this was tricky, watched the video many times to figure this change out
		fcsfile_comp

		#Cleaning
		fcsfile_comp_clean <- flow_auto_qc(fcsfile_comp) # will remove "bad events"
		fcsfile_comp_clean
		keyword(fcsfile_comp_clean) <- keyword(fcsfile)
		fcsfile_comp_clean
		??flowAI

		#Transformation
		trans <- estimateLogicle(fcsfile_comp_clean, colnames(fcsfile_comp_clean[,3:10]))
		fcsfile_comp_clean_trans <- transform(fcsfile_comp_clean, trans)

		#Visualize the results
		??ggcyto
		autoplot(fcsfile_comp_clean)
		autoplot(fcsfile_comp_clean_trans)
		autoplot(fcsfile_comp_clean_trans, x="FITC-A", y="PE-A")
		

# Stopped here, below is unfinished

		#In a flowSet
		myfiles <- list.files(path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files", pattern=".FCS$")
		fs <- flowCore::read.flowSet(myfiles, path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZU_files/")
		fs
		fs[[1]]
		spillover(fs[[1]])
		fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)
		fs_comp_clean <- flow_auto_qc(fs_comp)
		trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[[1]][,3:10]))
		fs_comp_clean_trans <- transform(fs_comp_clean, trans)

		#fsApply
		??fsApply
		fsApply(fs,each_col,median)