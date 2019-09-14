# R code for building this website

# TO RUN

#  R CMD BATCH --vanilla --quiet build_site_script.R 


###########

# add new plain markdown pages that aren't Rmd files here
	# so they can be rendered to html here

mdFiles <- c(
	"software/install_extra.md",
	"software/install_main.md",
	"pages/symbols.md",
	"pages/code_conduct.md",
	"pages/readings.md"
	)

	
#############################################	
	
# first, clean site
pkgdown::clean_site()

# make out files

# confusingly, render changes its internal working directory to match where the md is
	# so need to cut away directory locations
#mdOut <- sapply(mdFiles,function(x) rev(unlist(strsplit(x, split="/")))[[1]])

# OH but they need to be in the docs folder
wDir <- getwd()
mdOut <- paste0(wDir, "/docs/", mdFiles)
# replace md with html
mdOut <- gsub(".md$", ".html", mdOut)

# remake directories in docs
dir.create(path="docs/software")
dir.create(path="docs/pages")

# render them to html
for(i in 1:length(mdFiles)){
	rmarkdown::render(
		input = mdFiles[i],
		output_file = mdOut[i],
		quiet = TRUE, clean = TRUE
		)
	print(paste0("Created ", mdOut[i]))
	}

pkgdown::build_site(preview = TRUE)
