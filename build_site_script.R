# R code for building this website

# add new plain markdown pages that aren't Rmd files here
	# so they can be rendered to html here

mdFiles <- c(
	"./software/install_extra.md",
	"./software/install_main.md",
	"./pages/symbols.md",
	"./pages/code_conduct.md",
	"./pages/readings.md"
	)

	
#############################################	
	
# first, clean site
pkgdown::clean_site()

# make out files

# confusingly, render changes its internal working directory to match where the md is
	# so need to cut away directory locations
mdOut <- sapply(mdFiles,function(x) rev(unlist(strsplit(x, split="/")))[[1]])
# replace md with html
mdOut <- gsub(".md$", ".html", mdOut)
	
# render them to html
for(i in 1:length(mdFiles)){
	rmarkdown::render(
		input = mdFiles[i],
		output_file = mdOut[i]
		)
	}

pkgdown::build_site()
