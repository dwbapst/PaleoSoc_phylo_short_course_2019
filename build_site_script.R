# R code for building this website

# add new plain markdown pages that aren't Rmd files here
	# so they can be rendered to html here

mdFiles <- c(
	"software/install_main.md"
	)


	
	
#############################################	
	
# first, clean site
pkgdown::clean_site()

# make out files
mdOut <- lapply(mdFiles,
	
# render them to html
rmarkdown::render(
	input=mdFiles[i],
	output_file="install_main.html"
	)

