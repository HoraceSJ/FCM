# Flow Cytometry (FCM) Code
# Andrew Black
# Last updated 08-18-2021

# Adapted from
	# http://rprops.github.io/PhenoFlow/

# Step 1: Install R tools, here is tutorial:
	# https://cran.rstudio.com/bin/windows/Rtools/

	# After installation is complete, you need to perform one more step to be able to compile R packages: you need to put the location of the Rtools make utilities (bash, make, etc) on the PATH.
	writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
	# Now restart R, and verify that make can be found, which should show the path to your Rtools installation.
	Sys.which("make")
	## should show the following: "C:\\rtools40\\usr\\bin\\make.exe"

# Step 2: Install the required packages, first set:
	install.packages("mclust")
	install.packages("vegan")
	install.packages("MESS")
	install.packages("multcomp")
	install.packages("KernSmooth")
	install.packages("mvtnorm")
	install.packages("lattice")
	install.packages("survival")
	install.packages("TH.data")
	source("https://bioconductor.org/biocLite.R")
	biocLite("flowCore")
	source("https://bioconductor.org/biocLite.R")
	biocLite("flowViz")

# Step 3: Install second set of required packages:
	 packages = c("Biobase",
	 			"BiocGenerics",
	 			"flowViz",
	 			"flowFP",
	 			"graphics",
	 			"grDevices",
	 			"methods",
	 			"stats",
	 			"stats4",
	 			"MASS",
	 			"multcomp",
	 			"devtools")

	 ## Now load or install&load all
		package.check <- lapply(
		  packages,
		  FUN = function(x) {
		    if (!require(x, character.only = TRUE)) {
		      install.packages(x, dependencies = TRUE)
		      library(x, character.only = TRUE)
		    }
		  }
		)

# Step 4: Install flowCore package
	#source of code below: (http://bioconductor.org/packages/release/bioc/html/flowCore.html)
		if (!requireNamespace("BiocManager", quietly = TRUE))
		    install.packages("BiocManager")
		BiocManager::install("flowCore")

# Step 5: Install flowFDAE
	# Then install flowFDAE using the install_github function in the devtools package. (With build_vignettes=TRUE, the vignettes will be built and installed.) You first need to install the flowFDADataExample package for this purpose

		library(devtools)
		install_github("lievenclement/flowFDAExampleData")
		install_github("lievenclement/flowFDA", build_vignettes=TRUE)

# For help, see document:
		vignette("flowFDAExampleData")

		# The example dataset can be loaded using the data function. The fset data object is a flowCore flowSet object with the raw flow cytometric data of the experiment (Ellis et al., 2013).
			library(flowFDAExampleData)
			data(fset)
			fset

		# The following commands were used to create the fset data included with this package.
			library(flowFDA)
			fset<-read.flowSet(path="~/Dropbox/LabMet/flowcytometry/stress_test_2/",
			transformation=FALSE)
			fset
			##subset feet to reduce memory footprint
			param=c("SS Log","FL 1 Log","FL 3 Log")
			fset=fset[,param]
			fset


		# We will use the channels SS Log, FL 1 Log and FL 3 Log, which correspond to the side scatter, SYBR green and Propidium Iodide staining bandpass filters. The data have been transformed to fall within a range of 0 and 1 and extremely low intensities are removed using a rectangular gate so as to avoid artefacts. The flowcytometer used in this study returned log transformed intensities and had a maximum log transformed intensity of 2^16.
.
			mytrans<-function(x) x/2^16
			fset<-transform("FL 1 Log"=mytrans,"FL 3 Log"=mytrans,"SS Log"=mytrans)%on%fset
			rg <- rectangleGate(filterId="myRectGate", list("SS Log"=c(1/2^17, Inf), "FL 1 Log"=c(1/2^17, Inf),"FL 3 Log"=c(1/2^17,Inf)))
			fset<-Subset(fset,rg)


		# A good choice of the filename can enable an automated construction of the grouping variable
			#construct experiment factor
				files<-list.files(path="~/Dropbox/LabMet/flowcytometry/stress_test_2/",pattern=".fcs")
				expHlp<-unlist(strsplit(files,split="_replicate"))
				dim(expHlp)<-c(2,length(fset))
				group<-as.factor(expHlp[1,])
				nGroup<-nlevels(group)



# Begin here at normal session start, once entire installation above is done:

	packages_to_load <- c("mclust",
						"vegan",
						"MESS",
						"multcomp",
						"KernSmooth",
						"mvtnorm",
						"lattice",
						"survival",
						"TH.data",
						"flowCore",
						"flowViz",
						"Biobase",
						"BiocGenerics",
						"flowFP",
	 					"graphics",
	 					"grDevices",
	 					"methods",
	 					"stats",
	 					"stats4",
	 					"MASS",
	 					"devtools"
						)
	## Now load or install&load all
		package.check <- lapply(
		  packages_to_load,
		  FUN = function(x) {
		    if (!require(x, character.only = TRUE)) {
		      install.packages(x, dependencies = TRUE)
		      library(x, character.only = TRUE)
		    }
		  }
		)

# Test out with "IDR_Biofilter_backwash_try1.csv"