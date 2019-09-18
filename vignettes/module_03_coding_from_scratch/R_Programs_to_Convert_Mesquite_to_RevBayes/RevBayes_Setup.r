ZERO <- 1e-323;
MINEXPN <- 10^-10;
MINNO <- 5e-324;
MAXNO <- 1.797693e+308;
UNKNOWN <- -11;
INAP <- -22;
hell_no <- F;
taxonomic_rank <- c("subspecies","species","subgenus","genus","tribe","subfamily","family");
uncertains <- c("cf.","aff.","ex_gr.");

	##### ROUTINES TO READ CHARACTER MATRIX INFORMATION IN NEXUS FILE #######
# clean nexus file of characters that R hates....
accio_clade_reunion_given_paleodb_data <- function(clade_members)	{
paleodb_taxonomic_data <- accio_taxonomic_data_for_list_of_taxa(taxon_list=clade_members);
basic_taxonomic_data <- evanesco_na_from_matrix(data=paleodb_taxonomic_data$entered_taxa,replacement="");
paleodb_stored_ranks <- c("phylum","class","order","family","genus");
paleodb_not_entered_default <- c("NO_PHYLUM_SPECIFIED","NO_CLASS_SPECIFIED","NO_ORDER_SPECIFIED ","NO_FAMILY_SPECIFIED");
relv_columns <- (1:ncol(basic_taxonomic_data))[colnames(basic_taxonomic_data) %in% paleodb_stored_ranks];
basic_taxonomic_data[,relv_columns];
}

scourgify_nexus_text <- function(nexus_line)	{
nexus_line <- gsub("\xd4","",nexus_line);
nexus_line <- gsub("\xd5","",nexus_line);
nexus_line <- gsub("\x87","a",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\x93\xa0\xbc","ae",nexus_line);
nexus_line <- gsub("\xfc\xbe\x99\x83\xa0\xbc","c",nexus_line);
nexus_line <- gsub("\x8e","e",nexus_line);
nexus_line <- gsub("\x8f","e",nexus_line);
nexus_line <- gsub("\x92","i",nexus_line);
nexus_line <- gsub("\xbf","o",nexus_line);
nexus_line <- gsub("\x9a","o",nexus_line);
nexus_line <- gsub("\x97","o",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8e\x93\xa4\xbc","o",nexus_line);
nexus_line <- gsub("\x9f","ue",nexus_line);
nexus_line <- gsub("\xd0","-",nexus_line);
nexus_line <- gsub("\xd2","\"",nexus_line);
nexus_line <- gsub("\xd3","\"",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x86\x90\xbc","\'",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x86\x8c\xbc","ƒ",nexus_line);
nexus_line <- gsub("\xa7","ß",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8c\xa6\x88\xbc","≤",nexus_line);
nexus_line <- gsub("\xb3","≥",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x96\x8c\xbc","≈",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\xa6\x98\xbc","˚",nexus_line);
nexus_line <- gsub("\xb6","∂",nexus_line);
nexus_line <- gsub("\xc6","∆",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\xb6\x88\xbc","∑",nexus_line);
nexus_line <- gsub("\xfc\xbe\x99\x86\x88\xbc","Ω",nexus_line);
nexus_line <- gsub("\xa5"," ",nexus_line);
nexus_line <- gsub("Á","A",nexus_line);
nexus_line <- gsub("Ä","A",nexus_line);
nexus_line <- gsub("ä","a",nexus_line);
nexus_line <- gsub("á","a",nexus_line);
nexus_line <- gsub("å","a",nexus_line);
nexus_line <- gsub("Ç","C",nexus_line);
nexus_line <- gsub("ç","c",nexus_line);
nexus_line <- gsub("č","c",nexus_line);
nexus_line <- gsub("é","e",nexus_line);
nexus_line <- gsub("è","e",nexus_line);
nexus_line <- gsub("ê","e",nexus_line);
nexus_line <- gsub("ė","e",nexus_line);
nexus_line <- gsub("î","i",nexus_line);
nexus_line <- gsub("Î","I",nexus_line);
nexus_line <- gsub("ñ","n",nexus_line);
nexus_line <- gsub("Ö","O",nexus_line);
nexus_line <- gsub("Ø","O",nexus_line);
nexus_line <- gsub("ø","o",nexus_line);
nexus_line <- gsub("ó","o",nexus_line);
nexus_line <- gsub("ö","o",nexus_line);
nexus_line <- gsub("õ","o",nexus_line);
nexus_line <- gsub("Š","S",nexus_line);
nexus_line <- gsub("š","s",nexus_line);
nexus_line <- gsub("ů","u",nexus_line);
nexus_line <- gsub("ü","u",nexus_line);
nexus_line <- gsub("’","’",nexus_line);
nexus_line <- gsub("\x88","a",nexus_line);
nexus_line <- gsub("Ã„","A",nexus_line);
nexus_line <- gsub("Á","A",nexus_line);
nexus_line <- gsub("Ã¡","a",nexus_line);
nexus_line <- gsub("Ã¤","a",nexus_line);
nexus_line <- gsub("Ã¥","a",nexus_line);
nexus_line <- gsub("Ã§","c",nexus_line);
nexus_line <- gsub("Ã©","e",nexus_line);
nexus_line <- gsub("Ã¨","e",nexus_line);
nexus_line <- gsub("Ã±","n",nexus_line);
nexus_line <- gsub("Ã–","O",nexus_line);
nexus_line <- gsub("Ã¸","o",nexus_line);
nexus_line <- gsub("Ã¶","o",nexus_line);
nexus_line <- gsub("Ãµ","o",nexus_line);
nexus_line <- gsub("Ã´","o",nexus_line);
nexus_line <- gsub("Ã¼","u",nexus_line);
nexus_line <- gsub("√Æ","i",nexus_line);
nexus_line <- gsub("≈†","S",nexus_line);
nexus_line <- gsub("≈°","s",nexus_line);
nexus_line <- gsub("√•","a",nexus_line);
nexus_line <- gsub("&#367;","u",nexus_line);
nexus_line <- gsub("&#945;","α",nexus_line);
return(nexus_line);
}

# routine to read nexus file of Mesquite or Maclade format & return important infromation
accio_data_from_nexus_file <- function(nexus_file_name, polymorphs=T, UNKNOWN=-11, INAP=-22, rate_partitions="", trend_partitions="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
nexus <- scan(file=nexus_file_name,what=character(),sep="\n");
ml <- 0;
#i <- 1

for (i in 1:length(nexus))  {
	nexus[i] <- scourgify_nexus_text(nexus_line = nexus[i])
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]]
	if (length(j)>ml) ml <- length(j)
	}
ml <- ml+1;	# LENGTH OF LONGEST LINE
	
# file is now a vector of characters.  Turn it into a matrix with one char per cell
nexusfile <- matrix("\n",length(nexus),ml)
for (i in 1:length(nexus))  {
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]]
	for (k in 1:length(j))      nexusfile[i,k] <- j[k]
	if ((length(j)+2)<ml)
		for (k in (length(j)+2):ml) nexusfile[i,k] <- ""
	}

top <- 0;
ln <- 1;		# this is the row with the word "Matrix": character data starts next.
while (top==0)	{
	em_nexus <- gsub("\t","",nexus[ln]);
	nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
	if (!is.na(match("matrix",tolower(nexus_words))))	{
		top <- ln;
		}
	else	ln <- ln+1;
	}
top <- top+1;	# this will give the first row of data
# skip the comment text denoting character numbers (if present)
while(nexusfile[top,1]=="[" || nexusfile[top,1]==" ") top <- top+1

missing <- "?";
gap <- "-";
notu <- nchars <- strat <- range <- geog <- 0
for (i in 2:top)  {
	while ((nexusfile[i,1]=="[" || nexusfile[i,1]=="\n") && i<top)	i <- i+1;
	em_nexus <- gsub("\t","",nexus[i]);
	em_nexus <- gsub("="," = ",em_nexus);
	em_nexus <- gsub(";"," ; ",em_nexus);
	nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
	nexus_words <- nexus_words[nexus_words!=""];
	if (!is.na(match("ntax",tolower(nexus_words))) || !is.na(match("ntaxa",tolower(nexus_words))))	{
		j <- 1+match("ntax",tolower(nexus_words));
		if (is.na(j))	j <- 1+match("ntaxa",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		notu <- as.numeric(nexus_words[j]);
		}
	if (!is.na(match("nchar",tolower(nexus_words))) || !is.na(match("nchars",tolower(nexus_words))))	{
		j <- 1+match("nchar",tolower(nexus_words));
		if (is.na(j))	j <- 1+match("nchars",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		nchars <- as.numeric(nexus_words[j]);
		}
	if (!is.na(match("gap",tolower(nexus_words))))	{
		j <- 1+match("gap",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		gap <- nexus_words[j];
		}
	if (!is.na(match("missing",tolower(nexus_words))))	{
		j <- 1+match("missing",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		missing <- nexus_words[j];
		}
	if (!is.na(match("fa",tolower(nexus_words))))	{
		strat <- nexus_words[match("fa",tolower(nexus_words))-1];
		}
	if (!is.na(match("la",tolower(nexus_words))))	{
		range <- nexus_words[match("la",tolower(nexus_words))-1];
		}
	if (!is.na(match("geog",tolower(nexus_words))) || !is.na(match("geography",tolower(nexus_words))))	{
		geog <- nexus_words[match("la",tolower(nexus_words))-1];
		}
	}

extra <- 0;
if (strat>0)	{
	if (range>0)	{
		nchars <- nchars-2
		extra <- 2
		} else {
		nchars <- nchars-1
		extra <- 1
		} 
	strat_ranges <- matrix(0,notu,2)	
	}
if (geog>0)	{
	nchars <- nchars-1
	geography <- vector(length=notu)
	extra <- extra+1
	}
	
taxa <- vector(length=notu)
nstates <- array(0,dim=nchars)
chmatrix <- matrix(0,notu,nchars)
tx <- 1

# look for outgroup designation
exclude <- outgroup <- -1;
if (!is.na(match("BEGIN SETS;",nexus)))	{
	tx_pt <- match("BEGIN SETS;",nexus);	# look at taxon partitions
	look_for_outgroup <- TRUE;
	while (look_for_outgroup)	{
		tx_pt <- 1+tx_pt;
		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("-"," - ",yyy);
		yyy <- gsub("- "," - ",yyy);
		yyy <- gsub("  -  "," - ",yyy);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (!is.na(match("outgroup",tolower(xxx))))	{
			ttl_ln <- length(xxx);
			jj <- 1+match("outgroup",tolower(xxx));
			while (xxx[jj]==":" || xxx[jj]=="=")	jj <- jj+1;
			outgroup <- c();
			while (xxx[jj]!="," && jj<=ttl_ln)	{
				if (xxx[jj]=="-")	{
					jj <- jj+1;
					outgroup <- c(outgroup,((as.numeric(outgroup[length(outgroup)])+1):as.numeric(xxx[jj])));
					} else	{
					outgroup <- c(outgroup,xxx[jj]);
					}
				jj <- jj+1;
				}
			look_for_outgroup <- FALSE;
			} else	{
			if (tolower(nexus[tx_pt])=="end;" || tolower(nexus[tx_pt])=="\tend;")
				look_for_outgroup <- FALSE;
			}
		}

	# look for characters to exclude
	tx_pt <- match("BEGIN SETS;",nexus);
	xxx <- strsplit(paste(nexusfile[tx_pt-1,],collapse = "")," ");
	while(tolower(xxx[1])!="end")	{
		tx_pt <- tx_pt+1;
		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("- "," - ",yyy);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (length(xxx)==0 || is.na(xxx))
			xxx <- "";
#		if (!is.na(xxx) && !is.null(xxx) && xxx!="")	{
		if (xxx[1]=="charpartition")	{
			if (xxx[1]=="charpartition" && !is.na(match("exclude",tolower(xxx))))	{
				ttl_ln <- length(xxx);
				jj <- 1+match("exclude",tolower(xxx));
				while (xxx[jj]==":")	jj <- jj+1;
				exclude <- c();
				while (xxx[jj]!="," && jj<ttl_ln)	{
					if (xxx[jj]=="-")	{
						jj <- jj+1;
						exclude <- c(exclude,((as.numeric(exclude[length(exclude)])+1):as.numeric(xxx[jj])));
						} else	{
						exclude <- c(exclude,as.numeric(xxx[jj]));
						}
					jj <- jj+1;
					}
				}
			}
#		xxx[1];
#		tx_pt;
		}
	}

if (rate_partitions!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(rate_partitions,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(rate_partitions,breakup_this_line):length(breakup_this_line)];
			kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
			partition_names <- breakup_this_line[kk-1];
			kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
			character_rate_partitions <- rep("",nchars);
			for (pn in 1:length(partition_names))	{
				ll <- kk[pn]+1;
				this_part <- as.numeric(breakup_this_line[ll]);
				ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
				if (pn < length(partition_names))	{
					break_cell <- kk[pn+1]-1;
					} else	{
					break_cell <- kk[pn+1];
					}
				while (ll<break_cell)	{
					if (breakup_this_line[ll]=="-")	{
						ll <- ll+1;
						this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
						} else	{
						this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
						}
					ll <- ll+1;
					}
				character_rate_partitions[this_part] <- partition_names[pn];
				}
			got_splits<- T;
			}
		}
	} else	character_rate_partitions <- rep("imagine",nchars);

if (trend_partitions!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(trend_partitions,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(trend_partitions,breakup_this_line):length(breakup_this_line)];
			kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
			partition_names <- breakup_this_line[kk-1];
			kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
			character_trend_partitions <- rep("",nchars);
			for (pn in 1:length(partition_names))	{
				ll <- kk[pn]+1;
				this_part <- as.numeric(breakup_this_line[ll]);
				ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
				if (pn < length(partition_names))	{
					break_cell <- kk[pn+1]-1;
					} else	{
					break_cell <- kk[pn+1];
					}
				while (ll<break_cell)	{
					if (breakup_this_line[ll]=="-")	{
						ll <- ll+1;
						this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
						} else	{
						this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
						}
					ll <- ll+1;
					}
				character_trend_partitions[this_part] <- partition_names[pn];
				}
			got_splits<- T;
			}
		}
	} else	character_trend_partitions <- rep("square",nchars);

state_orders <- rep("unordered",nchars);

if (!is.na(match("BEGIN ASSUMPTIONS;",nexus)))	{
	tx_pt <- 1+match("BEGIN ASSUMPTIONS;",nexus);	# look at taxon partitions
	while (tolower(nexus[tx_pt])!="end;")	{
#		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("- "," - ",nexus[tx_pt]);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (!is.na(match("ord:",tolower(xxx))) && !is.na(match("revbayes",tolower(xxx))))	{
			ttl_ln <- length(xxx);
			jj <- 1+match("ord:",xxx);
			while (xxx[jj]==":")	jj <- jj+1;
			ordered <- c();
			while (xxx[jj]!="," && jj<=ttl_ln)	{
				if (xxx[jj]=="-")	{
					jj <- jj+1;
					ordered <- c(ordered,((as.numeric(ordered[length(ordered)])+1):as.numeric(xxx[jj])));
					} else	{
					ordered <- c(ordered,as.numeric(xxx[jj]));
					}
				jj <- jj+1;
				}
			state_orders[ordered] <- "ordered";
			}
		tx_pt <- 1+tx_pt;
		}
	}
mxln <- length(nexusfile[top,]);
s <- top;
# te all of the taxon names
for (tx in 1:notu)	{
	# first, read taxon name
	#### look for quotations###
	s <- top+tx-1
	if (nexusfile[s,1]=="'" || nexusfile[s,2]=="'")	{
		jj <- ((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"]);
		i <- max((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"])
		taxa[tx] <- pracma::strcat(nexusfile[s,(jj[1]+1):(jj[2]-1)])
		i <- i+1
		while (nexusfile[s,i]==" " && i<ncol(nexusfile))	i <- i+1
		}	else	{
		i <- 1
		if (nexusfile[s,1]!="\"")  {
			while (nexusfile[s,i]=="\t")	i <- i+1
			taxa[tx] <- nexusfile[s,i]
			i <- i+1
			while (nexusfile[s,i]!=" " && nexusfile[s,i]!='\t' && i<ncol(nexusfile))	{
				if (nexusfile[s,i]!="_")	{
					taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
					} else {
					taxa[tx] <- paste0(taxa[tx]," ")
					}
				i <- i+1
				}
			}	else {
			taxa[tx] <- nexusfile[s,2]
			i <- 3
			while (nexusfile[s,i]!=" " && nexusfile[s,i+1]!=" " && i<ncol(nexusfile))	{
				taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
				i <- i+1
				}
			}
		# now, get to characters
		while ((nexusfile[s,i]==" " || nexusfile[s,i]=="\t") && i<ncol(nexusfile))
			i <- i+1
		}
	k <- i;
	endline <- match("\n",nexusfile[s,])
	if (is.na(endline))	endline <- length(nexusfile[s,])
	if ((endline-k)==(nchars+extra))	{
		# true if there are no polymorphic characters for the taxon
		dummy <- nexusfile[s,k:(endline-1)]
		dummy[dummy==missing] <- UNKNOWN
		dummy[dummy==gap] <- INAP
		letterstate <- dummy
		dummy <- sapply(letterstate,switch_letter_state_to_numeric)
		chmatrix[tx,] <- as.numeric(dummy[1:nchars]);
		if (strat>0)	{
			strat_ranges[tx,1] <- strat_ranges[tx,2] <- as.numeric(dummy[strat])
			if (range>0)	strat_ranges[tx,2] <- as.numeric(dummy[range])
			}
		if (geog>0)	geography[tx]=as.numeric(nexusfile[geog,i])
		for (c in 1:nchars)	{
			if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1
			}
		} else	{
#		for (c in 1:(nchars+extra))	{
		c <- 0;
		while (c < (nchars+extra))	{
			c <- c+1;
			if (c<=nchars)	{
				if (nexusfile[s,i]=="(" || nexusfile[s,i]=="{")	{
					if (polymorphs==TRUE || polymorphs==1)	{
						i <- i+1
						w <- as.numeric(nexusfile[s,i])
						chmatrix[tx,c] <- -1*as.numeric(nexusfile[s,i])
						if ((1+w)>nstates[c])  nstates[c] <- 1+w
						i <- i+1
						j <- 1
						while (nexusfile[s,i]!=")" && nexusfile[s,i]!="}" && i<ncol(nexusfile))	{
							if (nexusfile[s,i]!="," && nexusfile[s,i]!=" ")	{
								w <- as.numeric(nexusfile[s,i])
								if ((w+1)>nstates[c])	nstates[c] <- w+1
								chmatrix[tx,c] <- chmatrix[tx,c]-((10^j)*w)
								i <- i+1
								j <- j+1
								} else {
								i <- i+1
								}
							}
						}	else {
						chmatrix[tx,c] <- UNKNOWN;
						while (nexusfile[s,i]!=')' && nexusfile[s,i]!="}")	i <- i+1;
						}
					} else if (nexusfile[s,i]==missing)	{
					chmatrix[tx,c] <- UNKNOWN;
					}	else if (nexusfile[s,i]==gap)	{
					chmatrix[tx,c] <- INAP;
					} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")  {
					chmatrix[tx,c] <- switch_letter_state_to_numeric(nexusfile[s,i]);
					}	else if (nexusfile[s,i]>="0" && nexusfile[s,i]<="9") {
					chmatrix[tx,c] <- as.numeric(nexusfile[s,i]);
					}
				if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1
				i <- i+1
				}  else {
				if (c==strat)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,1]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,1]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					if (range==0)	strat_ranges[tx,2] <- strat_ranges[tx,1]
					i <- i+1
					} else if (c==range)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,2]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,2]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					i <- i+1
					} else if (c==geog)	{
						if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
							geography[tx]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
							geography[tx]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					}
				} # end non-morphological data
#			print(nexusfile[s,k:83]);
#			print(chmatrix[tx,])
			if (nexusfile[s,i+1]=="\n" || i==(mxln-1)) c <- nchars+extra;
			}
		}
#	chmatrix[tx,];
#	tx <- tx+1;
#	s <- s+1
	}
	#x <- list(taxa,chmatrix,strat_ranges,geography)
	#return (list(taxa,chmatrix,strat_ranges,geography))

chmatrix <- scourgify_character_matrix(chmatrix,minst=0,UNKNOWN,INAP);	# clean up coding
nstates <- count_states(chmatrix,UNKNOWN,INAP);

tree_found <- 0;
while (s<length(nexus) && tree_found==0)	{
	while (nexus[s]!= "BEGIN TREES; " && s<length(nexus))
		s <- s+1;
	if (s<length(nexus))	{
		while (tree_found==0 && s<length(nexus))	{
			s <- s+1
			jj <- strsplit(nexus[s],split=c("\t"," "),fixed=TRUE)[[1]];
			jj <- paste(jj,collapse="")
			jj <- strsplit(jj,split=" ",fixed=TRUE)[[1]];
			if (sum(jj=="TREE")>0 || sum(jj=="tree")>0)	tree_found <- 1;
			}
#		s <- s+notu;
#		while (jj[i]=="")	jj[i] <- NULL;
#		while (j[1]=="\t")	j <- j[2:length(j)];
#		if (j[1]=="T" && j[2]=="R" && j[3]=="E")	{
#			while (j!="(")	j <- j[2:length(j)];
#			}
		newick_string <- jj[length(jj)];
		tree <- read_newick_string(newick_string);
		tree_found <- 1
		s <- length(nexus);
		}
	}

row.names(chmatrix) <- taxa

unscored_taxa <- c();
for (n in 1:notu)	{
	if (sum(chmatrix[n,]==UNKNOWN)==nchars)
		unscored_taxa <- c(unscored_taxa,n);
	}

if (nchars<10)	{
	colnames(chmatrix) <- 1:nchars;
	} else if (nchars<100)	{
	colnames(chmatrix) <- c(paste(0,(1:9),sep=""),10:nchars);
	} else if (nchars<1000)	{
	colnames(chmatrix) <- c(paste(00,(1:9),sep=""),paste(0,(10:99),sep=""),100:nchars);
	}
if (exclude[1]!=-1)	{
	keepers <- (1:nchars)[!(1:nchars) %in% exclude];
	chmatrix <- chmatrix[,keepers];
	nstates <- nstates[keepers];
	state_orders <- state_orders[keepers];
	character_rate_partitions <- character_rate_partitions[keepers];
	}

if (strat!=0 && geog!=0 && tree_found==1)  {
	output <- list(taxa,chmatrix,nstates,state_orders,strat_ranges,geography,tree,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Order","Stratigraphic_Ranges","Geography","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	} else if (strat!=0)  {
	if (geog!=0)  {
		output <- list(taxa,chmatrix,nstates,state_orders,strat_ranges,geography,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Stratigraphic_Ranges","Geography","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else if (tree_found!=0)	{
		output <- list(taxa,chmatrix,nstates,state_orders,strat_ranges,tree,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Stratigraphic_Ranges","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else	{
		output <- list(taxa,chmatrix,nstates,state_orders,strat_ranges,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Stratigraphic_Ranges","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		}
	} else if (geog!=0)  {
	if (tree_found!=0)	{
		output <- list(taxa,chmatrix,nstates,state_orders,geography,tree,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Geography","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else	{
		output <- list(taxa,chmatrix,nstates,state_orders,geography,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Geography","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		}
	} else if (tree_found!=0) {
	output <- list(taxa,chmatrix,nstates,state_orders,tree,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Order","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	} else	{
	output <- list(taxa,chmatrix,nstates,state_orders,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Order","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	}

return(output)
}

# routine to prompt user for a nexus file of Mesquite or Maclade format & return important infromation
accio_data_from_chosen_nexus_file <- function(polymorphs=T, UNKNOWN=-11, INAP=-22, rate_partition="", trend_partition="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
#print("Choose the nexus file you with to analyze: ");
nexus_file_name <- file.choose();
nexus <- scan(file=nexus_file_name,what=character(),sep="\n");
ml <- 0;
#i <- 1

for (i in 1:length(nexus))  {
	nexus[i] <- scourgify_nexus_text(nexus_line = nexus[i])
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]]
	if (length(j)>ml) ml <- length(j)
	}
ml <- ml+1;	# LENGTH OF LONGEST LINE
	
# file is now a vector of characters.  Turn it into a matrix with one char per cell
nexusfile <- matrix("\n",length(nexus),ml)
for (i in 1:length(nexus))  {
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]]
	for (k in 1:length(j))      nexusfile[i,k] <- j[k]
	if ((length(j)+2)<ml)
		for (k in (length(j)+2):ml) nexusfile[i,k] <- ""
	}

top <- 0;
ln <- 1;		# this is the row with the word "Matrix": character data starts next.
while (top==0)	{
	em_nexus <- gsub("\t","",nexus[ln]);
	nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
	if (!is.na(match("matrix",tolower(nexus_words))))	{
		top <- ln;
		}
	else	ln <- ln+1;
	}
top <- top+1;	# this will give the first row of data
# skip the comment text denoting character numbers (if present)
while(nexusfile[top,1]=="[" || nexusfile[top,1]==" ") top <- top+1

missing <- "?";
gap <- "-";
notu <- nchars <- strat <- range <- geog <- 0
for (i in 2:top)  {
	while ((nexusfile[i,1]=="[" || nexusfile[i,1]=="\n") && i<top)	i <- i+1;
	em_nexus <- gsub("\t","",nexus[i]);
	em_nexus <- gsub("="," = ",em_nexus);
	em_nexus <- gsub(";"," ; ",em_nexus);
	em_nexus <- gsub(",","",em_nexus);
	nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
	nexus_words <- nexus_words[nexus_words!=""];
	if (!is.na(match("ntax",tolower(nexus_words))) || !is.na(match("ntaxa",tolower(nexus_words))))	{
		j <- 1+match("ntax",tolower(nexus_words));
		if (is.na(j))	j <- 1+match("ntaxa",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		notu <- as.numeric(nexus_words[j]);
		}
	if (!is.na(match("nchar",tolower(nexus_words))) || !is.na(match("nchars",tolower(nexus_words))))	{
		j <- 1+match("nchar",tolower(nexus_words));
		if (is.na(j))	j <- 1+match("nchars",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		nchars <- as.numeric(nexus_words[j]);
		}
	if (!is.na(match("gap",tolower(nexus_words))))	{
		j <- 1+match("gap",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		gap <- nexus_words[j];
		}
	if (!is.na(match("missing",tolower(nexus_words))))	{
		j <- 1+match("missing",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		missing <- nexus_words[j];
		}
	if (!is.na(match("fa",tolower(nexus_words))))	{
		strat <- as.numeric(nexus_words[match("fa",tolower(nexus_words))-1]);
		}
	if (!is.na(match("la",tolower(nexus_words))))	{
		range <- as.numeric(nexus_words[match("la",tolower(nexus_words))-1]);
		}
	if (!is.na(match("geog",tolower(nexus_words))) || !is.na(match("geography",tolower(nexus_words))))	{
		geog <- c(nexus_words[match("geog",tolower(nexus_words))-1],nexus_words[match("geography",tolower(nexus_words))-1]);
		geog <- as.numeric(geog[!is.na(geog)]);
		}
	}

extra <- 0;
if (strat>0)	{
	if (range>0)	{
		nchars <- nchars-2
		extra <- 2
		} else {
		nchars <- nchars-1
		extra <- 1
		} 
	strat_ranges <- data.frame(FA=as.numeric(rep(0,notu)),LA=as.numeric(rep(0,notu)));
	}
if (geog>0)	{
	nchars <- nchars-1
	geography <- vector(length=notu)
	extra <- extra+1
	}
	
taxa <- vector(length=notu);
nstates <- array(0,dim=nchars);
chmatrix <- matrix(0,notu,nchars);
tx <- 1;

# look for outgroup designation
exclude <- outgroup <- -1;
if (!is.na(match("BEGIN SETS;",nexus)))	{
	tx_pt <- match("BEGIN SETS;",nexus);	# look at taxon partitions
	look_for_outgroup <- TRUE;
	while (look_for_outgroup)	{
		tx_pt <- 1+tx_pt;
		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("-"," - ",yyy);
		yyy <- gsub("- "," - ",yyy);
		yyy <- gsub("  -  "," - ",yyy);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (!is.na(match("outgroup",tolower(xxx))))	{
			ttl_ln <- length(xxx);
			jj <- 1+match("outgroup",tolower(xxx));
			while (xxx[jj]==":" || xxx[jj]=="=")	jj <- jj+1;
			outgroup <- c();
			while (xxx[jj]!="," && jj<=ttl_ln)	{
				if (xxx[jj]=="-")	{
					jj <- jj+1;
					outgroup <- c(outgroup,((as.numeric(outgroup[length(outgroup)])+1):as.numeric(xxx[jj])));
					} else	{
					outgroup <- c(outgroup,xxx[jj]);
					}
				jj <- jj+1;
				}
			look_for_outgroup <- FALSE;
			} else	{
			if (tolower(nexus[tx_pt])=="end;" || tolower(nexus[tx_pt])=="\tend;")
				look_for_outgroup <- FALSE;
			}
		}

	# look for characters to exclude
	tx_pt <- match("BEGIN SETS;",nexus);
	xxx <- strsplit(paste(nexusfile[tx_pt-1,],collapse = "")," ");
	while(tolower(xxx[1])!="end")	{
		tx_pt <- tx_pt+1;
		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("- "," - ",yyy);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (length(xxx)==0 || is.na(xxx))
			xxx <- "";
#		if (!is.na(xxx) && !is.null(xxx) && xxx!="")	{
		if (xxx[1]=="charpartition")	{
			if (xxx[1]=="charpartition" && !is.na(match("exclude",tolower(xxx))))	{
				ttl_ln <- length(xxx);
				jj <- 1+match("exclude",tolower(xxx));
				while (xxx[jj]==":")	jj <- jj+1;
				exclude <- c();
				while (xxx[jj]!="," && jj<ttl_ln)	{
					if (xxx[jj]=="-")	{
						jj <- jj+1;
						exclude <- c(exclude,((as.numeric(exclude[length(exclude)])+1):as.numeric(xxx[jj])));
						} else	{
						exclude <- c(exclude,as.numeric(xxx[jj]));
						}
					jj <- jj+1;
					}
				}
			}
#		xxx[1];
#		tx_pt;
		}
	}

# look for rate variation partitions
if (rate_partition!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(rate_partition,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			nexus[ln] <- gsub(":"," : ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  :"," :",nexus[ln]);
			nexus[ln] <- gsub(":  ",": ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(rate_partition,breakup_this_line):length(breakup_this_line)];
			kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
			partition_names <- breakup_this_line[kk-1];
			kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
			character_rate_partitions <- rep("",nchars);
			for (pn in 1:length(partition_names))	{
				ll <- kk[pn]+1;
				this_part <- as.numeric(breakup_this_line[ll]);
				ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
				if (pn < length(partition_names))	{
					break_cell <- kk[pn+1]-1;
					} else	{
					break_cell <- kk[pn+1];
					}
				while (ll<break_cell)	{
					if (breakup_this_line[ll]=="-")	{
						ll <- ll+1;
						this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
						} else	{
						this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
						}
					ll <- ll+1;
					}
				character_rate_partitions[this_part] <- partition_names[pn];
				}
			got_splits<- T;
			}
		}
	} else	character_rate_partitions <- rep("imagine",nchars);

if (trend_partition!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(trend_partition,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(trend_partition,breakup_this_line):length(breakup_this_line)];
			kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
			partition_names <- breakup_this_line[kk-1];
			kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
			character_trend_partitions <- rep("",nchars);
			for (pn in 1:length(partition_names))	{
				ll <- kk[pn]+1;
				this_part <- as.numeric(breakup_this_line[ll]);
				ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
				if (pn < length(partition_names))	{
					break_cell <- kk[pn+1]-1;
					} else	{
					break_cell <- kk[pn+1];
					}
				while (ll<break_cell)	{
					if (breakup_this_line[ll]=="-")	{
						ll <- ll+1;
						this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
						} else	{
						this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
						}
					ll <- ll+1;
					}
				character_trend_partitions[this_part] <- partition_names[pn];
				}
			got_splits<- T;
			}
		}
	} else	character_trend_partitions <- rep("square",nchars);

state_orders <- rep("unordered",nchars);

if (!is.na(match("BEGIN ASSUMPTIONS;",nexus)))	{
	tx_pt <- 1+match("BEGIN ASSUMPTIONS;",nexus);	# look at taxon partitions
	while (tolower(nexus[tx_pt])!="end;")	{
#		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("- "," - ",nexus[tx_pt]);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (!is.na(match("ord:",tolower(xxx))) && !is.na(match("revbayes",tolower(xxx))))	{
			ttl_ln <- length(xxx);
			jj <- 1+match("ord:",xxx);
			while (xxx[jj]==":")	jj <- jj+1;
			ordered <- c();
			while (xxx[jj]!="," && jj<=ttl_ln)	{
				if (xxx[jj]=="-")	{
					jj <- jj+1;
					ordered <- c(ordered,((as.numeric(ordered[length(ordered)])+1):as.numeric(xxx[jj])));
					} else	{
					ordered <- c(ordered,as.numeric(xxx[jj]));
					}
				jj <- jj+1;
				}
			state_orders[ordered] <- "ordered";
			}
		tx_pt <- 1+tx_pt;
		}
	}
mxln <- length(nexusfile[top,]);
s <- top;
# te all of the taxon names
for (tx in 1:notu)	{
	# first, read taxon name
	#### look for quotations###
	s <- top+tx-1;
	if (nexusfile[s,1]=="'" || nexusfile[s,2]=="'")	{
		jj <- ((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"]);
		i <- max((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"])
		taxa[tx] <- pracma::strcat(nexusfile[s,(jj[1]+1):(jj[2]-1)])
		i <- i+1
		while (nexusfile[s,i]==" " && i<ncol(nexusfile))	i <- i+1
		}	else	{
		i <- 1
		if (nexusfile[s,1]!="\"")  {
			while (nexusfile[s,i]=="\t")	i <- i+1
			taxa[tx] <- nexusfile[s,i]
			i <- i+1
			while (nexusfile[s,i]!=" " && nexusfile[s,i]!='\t' && i<ncol(nexusfile))	{
				if (nexusfile[s,i]!="_")	{
					taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
					} else {
					taxa[tx] <- paste0(taxa[tx]," ")
					}
				i <- i+1
				}
			}	else {
			taxa[tx] <- nexusfile[s,2]
			i <- 3
			while (nexusfile[s,i]!=" " && nexusfile[s,i+1]!=" " && i<ncol(nexusfile))	{
				taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
				i <- i+1
				}
			}
		# now, get to characters
		while ((nexusfile[s,i]==" " || nexusfile[s,i]=="\t") && i<ncol(nexusfile))
			i <- i+1
		}
	k <- i;
	endline <- match("\n",nexusfile[s,])
	if (is.na(endline))	endline <- length(nexusfile[s,])
	if ((endline-k)==(nchars+extra))	{
		# true if there are no polymorphic characters for the taxon
		dummy <- nexusfile[s,k:(endline-1)]
		dummy[dummy==missing] <- UNKNOWN
		dummy[dummy==gap] <- INAP
		letterstate <- dummy
		dummy <- sapply(letterstate,switch_letter_state_to_numeric)
		chmatrix[tx,] <- as.numeric(dummy[1:nchars]);
		if (strat>0)	{
			strat_ranges$FA[tx] <- strat_ranges$LA[tx] <- as.numeric(dummy[strat])
			if (range>0)	strat_ranges$LA[tx] <- as.numeric(dummy[range])
			}
		if (geog>0)	geography[tx]=as.numeric(nexusfile[geog,i])
		for (c in 1:nchars)	{
			if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1
			}
		} else	{
#		for (c in 1:(nchars+extra))	{
		c <- 0;
		while (c < (nchars+extra))	{
			c <- c+1;
			if (c<=nchars)	{
				if (nexusfile[s,i]=="(" || nexusfile[s,i]=="{")	{
					if (polymorphs==TRUE || polymorphs==1)	{
						i <- i+1
						w <- as.numeric(nexusfile[s,i])
						chmatrix[tx,c] <- -1*as.numeric(nexusfile[s,i])
						if ((1+w)>nstates[c])  nstates[c] <- 1+w
						i <- i+1
						j <- 1
						while (nexusfile[s,i]!=")" && nexusfile[s,i]!="}" && i<ncol(nexusfile))	{
							if (nexusfile[s,i]!="," && nexusfile[s,i]!=" ")	{
								w <- as.numeric(nexusfile[s,i])
								if ((w+1)>nstates[c])	nstates[c] <- w+1
								chmatrix[tx,c] <- chmatrix[tx,c]-((10^j)*w)
								i <- i+1
								j <- j+1
								} else {
								i <- i+1
								}
							}
						}	else {
						chmatrix[tx,c] <- UNKNOWN;
						while (nexusfile[s,i]!=')' && nexusfile[s,i]!="}")	i <- i+1;
						}
					} else if (nexusfile[s,i]==missing)	{
					chmatrix[tx,c] <- UNKNOWN;
					}	else if (nexusfile[s,i]==gap)	{
					chmatrix[tx,c] <- INAP;
					} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")  {
					chmatrix[tx,c] <- switch_letter_state_to_numeric(nexusfile[s,i]);
					}	else if (nexusfile[s,i]>="0" && nexusfile[s,i]<="9") {
					chmatrix[tx,c] <- as.numeric(nexusfile[s,i]);
					}
				if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1
				i <- i+1
				}  else {
				if (c==strat)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,1]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,1]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					if (range==0)	strat_ranges[tx,2] <- strat_ranges[tx,1]
					i <- i+1
					} else if (c==range)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,2]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,2]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					i <- i+1
					} else if (c==geog)	{
						if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
							geography[tx]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
							geography[tx]=switch_letter_state_to_numeric(nexusfile[s,i])
						}
					}
				} # end non-morphological data
#			print(nexusfile[s,k:83]);
#			print(chmatrix[tx,])
			if (nexusfile[s,i+1]=="\n" || i==(mxln-1)) c <- nchars+extra;
			}
		}
#	chmatrix[tx,];
#	tx <- tx+1;
#	s <- s+1
	}

if (trend_partition!="")	print("Choose a file giving the order in which taxa appear (with low numbers = early): ");
chmatrix <- scourgify_character_matrix(chmatrix,minst=0,UNKNOWN,INAP);	# clean up coding
if (trend_partition!="")	{
	appearance_order_file <- file.choose();
	appearance_order <- read.csv(appearance_order_file,header=T);
	otu_fas <- match(appearance_order$appearance,sort(unique(appearance_order$appearance)));
	for (nch in 1:ncol(chmatrix))	{
		otu_states <- chmatrix[,nch];
		chmatrix[,nch] <- rescore_states_by_first_appearances(otu_states,otu_fas,UNKNOWN,INAP);
		}
	}
nstates <- count_states(chmatrix,UNKNOWN,INAP);

tree_found <- 0;
while (s<length(nexus) && tree_found==0)	{
	while (nexus[s]!= "BEGIN TREES; " && s<length(nexus))
		s <- s+1;
	if (s<length(nexus))	{
		while (tree_found==0 && s<length(nexus))	{
			s <- s+1
			jj <- strsplit(nexus[s],split=c("\t"," "),fixed=TRUE)[[1]];
			jj <- paste(jj,collapse="")
			jj <- strsplit(jj,split=" ",fixed=TRUE)[[1]];
			if (sum(jj=="TREE")>0 || sum(jj=="tree")>0)	tree_found <- 1;
			}
		newick_string <- jj[length(jj)];
		tree <- read_newick_string(newick_string);
		tree_found <- 1
		s <- length(nexus);
		}
	}

row.names(chmatrix) <- taxa

unscored_taxa <- c();
for (n in 1:notu)	{
	if (sum(chmatrix[n,]==UNKNOWN)==nchars)
		unscored_taxa <- c(unscored_taxa,n);
	}

if (nchars<10)	{
	colnames(chmatrix) <- 1:nchars;
	} else if (nchars<100)	{
	colnames(chmatrix) <- c(paste(0,(1:9),sep=""),10:nchars);
	} else if (nchars<1000)	{
	colnames(chmatrix) <- c(paste(00,(1:9),sep=""),paste(0,(10:99),sep=""),100:nchars);
	}
if (exclude[1]!=-1)	{
	keepers <- (1:nchars)[!(1:nchars) %in% exclude];
	chmatrix <- chmatrix[,keepers];
	nstates <- nstates[keepers];
	state_orders <- state_orders[keepers];
	character_rate_partitions <- character_rate_partitions[keepers];
	character_trend_partitions <- character_trend_partitions[keepers];
	}

if (strat!=0 && geog!=0 && tree_found==1)  {
	output <- list(taxa,chmatrix,nstates,state_orders,strat_ranges,geography,tree,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Order","Stratigraphic_Ranges","Geography","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	} else if (strat!=0)  {
	if (geog!=0)  {
		output <- list(taxa,chmatrix,nstates,state_orders,strat_ranges,geography,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Stratigraphic_Ranges","Geography","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else if (tree_found!=0)	{
		output <- list(taxa,chmatrix,nstates,state_orders,strat_ranges,tree,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Stratigraphic_Ranges","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else	{
		output <- list(taxa,chmatrix,nstates,state_orders,strat_ranges,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Stratigraphic_Ranges","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		}
	} else if (geog!=0)  {
	if (tree_found!=0)	{
		output <- list(taxa,chmatrix,nstates,state_orders,geography,tree,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Geography","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else	{
		output <- list(taxa,chmatrix,nstates,state_orders,geography,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Order","Geography","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		}
	} else if (tree_found!=0) {
	output <- list(taxa,chmatrix,nstates,state_orders,tree,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Order","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	} else	{
	output <- list(taxa,chmatrix,nstates,state_orders,outgroup,unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Order","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	}

return(output)
}

# routine to count varying states for each character
count_states <- function(chmatrix,UNKNOWN=-11,INAP=-22,include_polymorphs=T)	{
if (!is.matrix(chmatrix))	{
	chmatrix <- data.frame(ch=chmatrix);
	}
nchars <- ncol(chmatrix);
nstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	unique_states <- char_states[char_states>=0]
	if (sum(char_states<0)>0 && include_polymorphs)	{
		polys <- char_states[char_states<0];
		polystates <- c();
		for (pp in polys)	{
			polystates <- c(polystates,unravel_polymorph(poly=pp));
			}
		unique_states <- sort(unique(c(unique_states,polystates)));
		}
	nstates <- c(nstates,length(unique_states));
	}	# pick up here!!!
return(nstates);
}

# list characters that have autapomorphic states;
list_autapomorphic_characters <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(chmatrix);
#scored <- vector(length=nch)
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- c()
for (ch in 1:nchars)	{
	if(sum(notus_per_chstate[ch,1:chstates[ch]]==1)>0)	{
		autaps <- c(autaps,ch);
		}
	}
return(autaps);
}

# count missing and/or inapplicable per otu
count_scored_otus_per_character_state <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
nch <- ncol(chmatrix);
#notu <- nrow(chmatrix);
#scored <- vector(length=nch)
scored <- array(0,dim=c(nch,max(chstates)))
for (c in 1:nch)	{
	for (st in 1:chstates[c])	{
		stt <- st-1
		scored[c,st] <- sum(chmatrix[,c]==stt)
		}
	}
return(scored);
}

# convert polytomies to individual states
unravel_polymorph <- function(poly)	{
combo <- -1*poly
sts <- 1+floor(log10(abs(combo)))
polymorphics <- vector(length=sts)

base <- 10^(sts-1)
for (s in 1:sts)	{
	polymorphics[s] <- floor(abs(combo)/base)
	combo <- combo%%base
	base <- base/10
	}
return (polymorphics)
}

count_characters_in_string <- function(string_to_count)	{
j <- strsplit(string_to_count,split="",fixed=TRUE)[[1]];
return(length(j));
}

scourgify_character_matrix <- function(chmatrix,minst=0,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(chmatrix)	# replaces spc to standardize coding.
nchars <- ncol(chmatrix)
for (ch in 1:nchars)	{
#ch <- 0
#while (ch<=nchars)	{
#	ch <- ch+1;
	rem <- c((1:notu)[chmatrix[,ch]==UNKNOWN],(1:notu)[chmatrix[,ch]==INAP]);
	relv <- (1:notu)[!(1:notu) %in% rem];
	relv_np <- relv[chmatrix[relv,ch]>=0];
	relv_pl <- relv[chmatrix[relv,ch]<0];
	if (length(relv)>0)	{
		test <- chmatrix[relv_np,ch];
		coded <- sort(unique(test[test>=0]));
		if (length(relv_pl)>0)	{
			polystates <- c();
			for (rp in relv_pl)	{
				polystates <- c(polystates,unravel_polymorph(poly=chmatrix[rp,ch]));
				}
			polystates <- sort(unique(polystates));
			coded <- sort(unique(c(coded,polystates)));
			}
		new_codes <- (1:length(coded))+(minst-1);
		chmatrix[relv_np,ch] <- new_codes[match(chmatrix[relv_np,ch],coded)];
		if (length(relv_pl)>0)	{
			for (rp in relv_pl)	{
				polycode <- unravel_polymorph(poly=chmatrix[rp,ch]);
				newpolycode <- new_codes[match(polycode,coded)];
				chmatrix[rp,ch] <- ravel_polymorph(newpolycode);
				}
			}
		}
	}
return(chmatrix)
}

#### convert (1,(2,3)) to vector_tree = 4 5 5 -1 4
#### 	where number is the htu number of the clade to which a species or htu belong
read_newick_string <- function(newick_string) {
nodes <- 0;
if (length(newick_string)==1)	newick_string <- strsplit(newick_string,split="",fixed=TRUE)[[1]];
for (i in 1:length(newick_string))		if (newick_string[i]=="(")	nodes <- nodes+1;
# get clades
clades <- vector(length=nodes);
for (c in 1:nodes)	clades[c] <- c;
# get taxa
notu <- p <- 0
for (i in 1:length(newick_string))	{
	if (newick_string[i]>="0" && newick_string[i]<="9")	{
		otu <- as.numeric(newick_string[i])+(otu * (10^p))
		p <- p+1
		if (otu>notu)	notu <- otu
		} else {
		p <- otu <- 0
		}
	}
vector_tree <- vector(length=notu+max(clades))
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(newick_string))	{
	if (newick_string[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (newick_string[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (newick_string[i]==",")	{
		sp <- p <- 0
		} else if (newick_string[i]>="0" && newick_string[i]<="9")	{
		sp <- as.numeric(newick_string[i])+(sp*(10^p))
		p <- p+1
		if (newick_string[i+1]<"0" || newick_string[i]>"9")	vector_tree[sp] <- notu+clades[c]
		}
	}

return(vector_tree)
}

	##### ROUTINES TO WRITE REVBAYES SCRIPTS #######
diffindo_character_matrix_by_state_numbers_and_ordering <- function(nexus_file_name,analysis_name="",write_data_directory="",read_data_directory,data_file_lead="data/",polymorphs=TRUE, UNKNOWN=-11, INAP=-22)	{
orig_nexus_file_name <- nexus_file_name;
if (analysis_name=="")	{
	j <- strsplit(orig_nexus_file_name,split="",fixed=TRUE)[[1]];
	j <- j[1:(max((1:length(j))[j %in% "."])-1)];
	analysis_name <- paste(j,collapse="");
	}

initial_data <- accio_data_from_nexus_file(nexus_file_name=paste(read_data_directory,nexus_file_name,sep=""),polymorphs,UNKNOWN,INAP);

n_states <- initial_data$States;
n_states[n_states<2] <- 2;
taxon_names <- initial_data$OTUs;
chmatrix <- initial_data$Matrix;
rownames(chmatrix) <- initial_data$OTUs;
n_chars <- length(n_states);
state_orders <- initial_data$State_Order;
notu <- nrow(chmatrix);
#if (write_rev_bayes_source)	{
nstates <- sort(unique(n_states));
state_numbers <- state_ordering <- matrix_file_names <- c();
#	}
#for (st in 1:length(state_numbers))	{
st <- 0;
while (st < length(nstates))	{
	st <- st+1;
	relv_chars <- (1:n_chars)[n_states==nstates[st]];
	orderings <- unique(state_orders[relv_chars]);
	for (sto in 1:length(orderings))	{
		relv_chars <- (1:n_chars)[n_states==nstates[st]];
		relv_chars <- relv_chars[state_orders[relv_chars]==orderings[sto]];
		ch_matrix <- initial_data$Matrix[,relv_chars];
		new_file_name <- paste(analysis_name,"_Matrix_",nstates[st],"_States_",orderings[sto],".nex",sep="");
		if (write_data_directory!="")
			new_file_name <- paste(write_data_directory,new_file_name,sep="");
		if (!is.matrix(ch_matrix))	{
			ch_matrix <- matrix(ch_matrix);
			rownames(ch_matrix) <- taxon_names;
			}
		scribio_rev_bayes_nexus_file_from_character_matrix(ch_matrix,new_file_name,no_states=nstates[st],UNKNOWN,INAP);
#		if (write_rev_bayes_source)	{
			matrix_file_names <- c(matrix_file_names,paste(data_file_lead,analysis_name,"_Matrix_",nstates[st],"_States_",orderings[sto],".nex",sep=""));
#			matrix_file_names <- c(matrix_file_names,new_file_name);
			state_numbers <- c(state_numbers,nstates[st]);
			state_ordering <- c(state_ordering,orderings[sto]);
#			}
		}
	}
output <- list(initial_data,matrix_file_names,state_numbers,state_ordering);
names(output) <- c("initial_data","matrix_file_names","state_numbers","state_ordering");
return(output);
}

diffindo_character_matrix_by_state_numbers_and_other_partitions <- function(analysis_name="",write_data_directory="",rate_partition="",trend_partition="",taxa_subset="",data_file_lead="data/",polymorphs=TRUE, UNKNOWN=-11, INAP=-22)	{
# rate_partition: name of acharacter group that should have separate rates (e.g., CHARPARTITION Functional_Partitions  =  Nonfeeding :  1- 45 58- 60, Feeding :  46- 57; in a Nexus file)
# rattrend_partition: name of character group that for which one group ("Driven") has biased change (e.g., Driven_Trend  =  Unbiased :  1- 32 34- 40 43- 50 52- 60, Driven :  33 41- 42 51; in a Nexus file)
# taxa_subset: vector listing the subset of species to include in the final matrices;
# write_data_directory: the directory to which the partitioned matrices should go
# data_file_lead: the directory addendum to paste in front of data files so that RevBayes can find them
if (analysis_name=="")	{
	analysis_name <- "Wombat_Rock";
	}
#initial_data <- accio_data_from_nexus_file(nexus_file_name=paste(read_data_directory,nexus_file_name,sep=""),polymorphs,UNKNOWN,INAP,character_partitions=rate_partition);
#if (rate_partition!="")	{
print("Choose a nexus file to convert to RevBayes preferred format:")
initial_data <- accio_data_from_chosen_nexus_file(polymorphs,UNKNOWN,INAP,rate_partition,trend_partition);
#	} else if (trend_partition!="")	{
#	initial_data <- accio_data_from_chosen_nexus_file(polymorphs,UNKNOWN,INAP,character_partitions=trend_partition);
#	} else	{
#	initial_data <- accio_data_from_chosen_nexus_file(polymorphs,UNKNOWN,INAP,character_partitions="");
#	}

n_states <- initial_data$States;
n_states[n_states<2] <- 2;
taxon_names <- initial_data$OTUs;
taxon_names <- gsub("\\(\\?\\)","",taxon_names);
taxon_names <- gsub("\\?","",taxon_names);
taxon_names <- gsub("n\\. sp\\. ","",taxon_names);
taxon_names <- gsub("n\\. gen\\. ","",taxon_names);
taxon_names <- gsub("  "," ",taxon_names);
initial_data$OTUs <- taxon_names;

chmatrix <- initial_data$Matrix;
rownames(chmatrix) <- initial_data$OTUs;
n_chars <- length(n_states);
state_orders <- initial_data$State_Order;
notu <- nrow(chmatrix);
#if (write_rev_bayes_source)	{
nstates <- sort(unique(n_states));
rate_partitions <- initial_data$Rate_Partitions;
trend_partitions <- initial_data$Trend_Partitions;

# separate out all unique combinations that we might use to partition the characters
if (length(unique(rate_partitions))>1 && length(unique(trend_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),ordering=as.character(state_orders),rate_partitions=as.character(rate_partitions),trend_partitions=as.character(trend_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$ordering),];
	} else if (length(unique(rate_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),ordering=as.character(state_orders),rate_partitions=as.character(rate_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$ordering),];
	} else if (length(unique(trend_partitions))>1)	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),ordering=as.character(state_orders),trend_partitions=as.character(trend_partitions),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$ordering),];
	} else	{
	partition_combos <- data.frame(n_states=as.numeric(n_states),ordering=as.character(state_orders),stringsAsFactors=hell_no);
	partition_combos <- unique(partition_combos);
	partition_combos <- partition_combos[order(partition_combos$n_states,partition_combos$ordering),];
	}

if (taxa_subset!="" && length(taxa_subset)<notu)	{
	chmatrix <- chmatrix[match(taxa_subset,rownames(chmatrix)),];
	notu <- nrow(chmatrix);
	}

pc <- 0;
character_rate_partitions <- character_trend_partitions <- state_numbers <- state_ordering <- matrix_file_names <- coding_bias <- c();
while (pc < nrow(partition_combos))	{
	# check coding_bias here!!!!
	pc <- pc+1;
	# find characters with the appropriate numbers of states
	relv_characters <- (1:n_chars)[n_states %in% partition_combos$n_states[pc]];
	if (length(unique(state_orders))>1)	{
		# find characters with this ordering
		right_ordering <- (1:n_chars)[state_orders %in% partition_combos$ordering[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_ordering];
		}
	if (length(unique(rate_partitions))>1)	{
		# find characters in this rate partitioning
		right_partition <- (1:n_chars)[rate_partitions %in% partition_combos$rate_partitions[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_partition];
		}
	if (length(unique(trend_partitions))>1)	{
		# find characters in this trend partitioning
		right_partition <- (1:n_chars)[trend_partitions %in% partition_combos$trend_partitions[pc]];
		relv_characters <- relv_characters[relv_characters %in% right_partition];
		}
	# look for invariant & autapomorphic characters
	chmatrix_red <- chmatrix[,relv_characters];
	if (!is.matrix(chmatrix_red))
		chmatrix_red <- as.matrix(chmatrix_red);
	## RevBayes ignores polymorphs; so, if the only cases of derived states are in polymorphic taxa, recode them to derived states;
	observed_states_sans_poly <- count_states(chmatrix_red,UNKNOWN,INAP,include_polymorphs = F);
	poss_polymorph_fixs <- (1:ncol(chmatrix_red))[observed_states_sans_poly<partition_combos$n_states[pc]];
	ppf <- 0;
	while (ppf < length(poss_polymorph_fixs))	{
		ppf <- ppf+1;
		nch <- poss_polymorph_fixs[ppf];
		poly_otus <- (1:notu)[!chmatrix_red[,nch] %in% c(INAP,UNKNOWN,0:partition_combos$n_states[pc])];
		orig_state <- obs_states <- sort(unique(chmatrix_red[chmatrix_red[,nch]>=0,nch]));
		potu <- 0;
		while (potu < length(poly_otus))	{
			potu <- potu+1;
			polystates <- unravel_polymorph(chmatrix_red[poly_otus[potu],nch]);
			if (length(polystates[!polystates %in% obs_states])>0)	{
				chmatrix_red[poly_otus[potu],nch] <- max(polystates[!polystates %in% obs_states]);
				} else if (length(polystates[!polystates %in% orig_state])>0)	{
				chmatrix_red[poly_otus[potu],nch] <- min(polystates[!polystates %in% orig_state]);
				}
			}
		}
	observed_states <- count_states(chmatrix_red,UNKNOWN,INAP);
	
	if (sum(observed_states==1)>0)	{
		coding_bias <- c(coding_bias,"all");
		} else	{
		chstates <- initial_data$States;
		autaps <- list_autapomorphic_characters(chmatrix,chstates,UNKNOWN,INAP);
		if (length(autaps)>0 || partition_combos$n_states[pc]>2)	{
			coding_bias <- c(coding_bias,"variable");
			} else	{
			coding_bias <- c(coding_bias,"informative");
			}
		}
	new_file_name <- paste(analysis_name,"_Matrix_",partition_combos$n_states[pc],"_States",sep="");
	if (length(unique(state_orders))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$ordering[pc],sep="");
	if (length(unique(rate_partitions))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$rate_partitions[pc],sep="");
	if (length(unique(trend_partitions))>1)
		new_file_name <- paste(new_file_name,"_",partition_combos$trend_partitions[pc],sep="");
	orig_file_name <- new_file_name <- paste(new_file_name,".nex",sep="");
	if (write_data_directory!="")
		new_file_name <- paste(write_data_directory,new_file_name,sep="");
#	if (!is.matrix(ch_matrix))	{
#		ch_matrix <- matrix(ch_matrix);
#		rownames(ch_matrix) <- taxon_names;
#		}
	scribio_rev_bayes_nexus_file_from_character_matrix(ch_matrix=chmatrix_red,new_file_name,no_states=partition_combos$n_states[pc],UNKNOWN,INAP);
	matrix_file_names <- c(matrix_file_names,paste(data_file_lead,orig_file_name,sep=""));
#			matrix_file_names <- c(matrix_file_names,new_file_name);
	}

state_numbers <- partition_combos$n_states;
state_ordering <- partition_combos$ordering;
if (!is.null(partition_combos$rate_partitions))	{
	character_rate_partitions <- partition_combos$rate_partitions;
	} else	{
	character_rate_partitions <- rep("imagine",nrow(partition_combos));
	}
if (!is.null(partition_combos$trend_partitions))	{
	character_trend_partitions <- partition_combos$trend_partitions;
	} else	{
	character_trend_partitions <- rep("square",nrow(partition_combos));
	}

output <- list(initial_data,matrix_file_names,state_numbers,state_ordering,character_rate_partitions,character_trend_partitions,coding_bias);
names(output) <- c("initial_data","matrix_file_names","state_numbers","state_ordering","rate_partitions","trend_partitions","coding_bias");
return(output);
}

#scribio_rev_bayes_script(analysis_name,taxon_names=otu_names,matrix_file_names,state_numbers,state_ordering,outgroup_taxa,unscored_taxa,fbd_parameterization_script,no_runs=4,write_scripts_directory=write_scripts_directory,set_wdir)
#scribio_Rev_Bayes_script_for_partitioned_character_data(            analysis_name,   initial_data,matrix_file_names,state_numbers,state_ordering,coding,           write_scripts_directory,fbd_parameterization_script,character_rate_partitions,fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);
scribio_Rev_Bayes_script_for_partitioned_character_data <- function(analysis_name="",initial_data,matrix_file_names,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script="",character_rate_partitions="",character_trend_partitions="",fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",data_file_lead="data/",no_runs=4)	{
# analysis_name: name that specificies this particular analysis.  (The clade name is a good choice)
# initial_data: data from nexus file.  
# matrix_file_names: a vector giving the list of all character matrices to be used
# state_numbers: vector giving the number of states for each matrix in matrix_file_names
# state_ordering: vector designating unordered or ordered state evolution
# write_scripts_directory: tell R where to send script file
# FBD: if "true", then add script to initiate FBD analyses from a separate file
# set_wdir: tell RevBayes where to set working directory
# output_file_lead: tell RevBayes script where to send/find output (default: "output/")
# script_file_lead: tell RevBayes script where to send/find scripts (default: "script/")
# no_runs: tell RevBayes how many runs to execut (default=4)
#filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Partitioned_Analysis.Rev",sep="");
filename <- paste(write_scripts_directory,analysis_name,sep="");
if (length(unique(character_rate_partitions))>1)
	filename <- paste(filename,"_Rate_Partitioned",sep="");
if (length(unique(character_trend_partitions))>1)
	filename <- paste(filename,"_Driven_Trend_Test",sep="");
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename <- paste(filename,"_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename <- paste(filename,"_FBD_Analysis.Rev",sep="");
	}
revbayes_source <- "clear();"
if (set_wdir!="")	{
	revbayes_source <- c(revbayes_source,paste("setwd(\"",set_wdir,"\");\t#CHANGE THIS TO THE FOLDER IN WHICH YOU HAVE REVBAYES SCRIPTS & DATA!!!",sep=""));
	revbayes_source <- c(revbayes_source,"### This director needs three subdirectories (folders):");
	revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/scripts (additional RevBayes routines that will be used)");
	revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/data (holds data matrices & taxonomic information)");
	revbayes_source <- c(revbayes_source,"#     RevBayes_Projects/output (where trees & logs will be sent)");
	}
revbayes_source <- c(revbayes_source,paste("source(\"",paste(script_file_lead,"Imperio_Default_Settings.Rev",sep=""),"\");",sep=""));
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"###############################################################################");
revbayes_source <- c(revbayes_source,"# This is (these are) the nexus file(s) that you are using for this analysis  #");
revbayes_source <- c(revbayes_source,"#     Make sure that filenames & directories are correct!!!");                #
revbayes_source <- c(revbayes_source,"###############################################################################");
file_name_string <- paste(matrix_file_names, collapse = "\", \"");
file_name_string <- paste("filenames <- v(\"",file_name_string,"\");",sep="");
file_name_string <-gsub("~/","",file_name_string);
revbayes_source <- c(revbayes_source,file_name_string);
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,paste("partition_states <- v(",paste(state_numbers, collapse = ","),");",sep=""));
revbayes_source <- c(revbayes_source,paste("partition_ordering <- v(\"",paste(state_ordering, collapse = "\",\""),"\");",sep=""));
revbayes_source <- c(revbayes_source,paste("coding_bias <- v(\"",paste(coding_bias, collapse = "\",\""),"\");\t## prepare for ascertainment bias in binary characters; 'all': invariant & autapomorphies present; 'variable': all vary & autapomorphies present; 'informative': all vary & no autapomorphies.",sep=""));
if (length(unique(character_rate_partitions))>1)	{
	revbayes_source <- c(revbayes_source,paste("rate_partitions <- v(\"",paste(character_rate_partitions, collapse = "\",\""),"\");\t# rate partition for each character paritition",sep=""));
	partition_labels <- unique(character_rate_partitions)
	revbayes_source <- c(revbayes_source,paste("rate_partition_labels <- v(\"",paste(partition_labels, collapse = "\",\""),"\");\t# names of rate partitions",sep=""));
	revbayes_source <- c(revbayes_source,paste("ttl_rate_partitions <- ",length(partition_labels),";\t\t\t\t\t\t# number of rate partitions among character types",sep=""));
	} else	{
#	revbayes_source <- c(revbayes_source,paste("rate_partitions <- v(\"",paste(rep("none",length(state_ordering)), collapse = "\",\""),"\");",sep=""));
#	revbayes_source <- c(revbayes_source,"rate_partition_labels[1] <- \"none\";\t\t# we are not partitioning characters for different general rates");
#	revbayes_source <- c(revbayes_source,paste("ttl_rate_partitions <- 1;\t\t# no rate rate partitions among characters",sep=""));
	}
if (length(unique(character_trend_partitions))>1)	{
	revbayes_source <- c(revbayes_source,paste("driven_trend_partitions <- v(\"",paste(tolower(character_trend_partitions), collapse = "\",\""),"\");\t# use 'driven' for characters with biased change",sep=""));
	partition_labels <- unique(tolower(character_trend_partitions));
	revbayes_source <- c(revbayes_source,paste("trend_partition_labels <- v(\"",paste(partition_labels, collapse = "\",\""),"\");\t# names of rate partitions",sep=""));
	revbayes_source <- c(revbayes_source,paste("ttl_trend_partitions <- ",length(partition_labels),";\t\t\t\t\t\t# number of rate partitions among character types",sep=""));
	} else	{
#	revbayes_source <- c(revbayes_source,paste("driven_trend_partitions <- v(\"",paste(rep("none",length(state_ordering)), collapse = "\",\""),"\");",sep=""));
#	revbayes_source <- c(revbayes_source,"trend_partition_labels[1] <- \"none\";\t\t# we are not partitioning characters for different general rates");
#	revbayes_source <- c(revbayes_source,paste("ttl_trend_partitions <- 1;\t\t# no rate rate partitions among characters",sep=""));
	}
revbayes_source <- c(revbayes_source,paste("max_age <- ",max_age,";\t\t\t\t\t\t# used if big_bang==TRUE;",sep=""));
revbayes_source <- c(revbayes_source,"");
taxon_names <- initial_data$OTUs;
notu <- length(initial_data$OTUs);
if(as.numeric(initial_data$Outgroup[1])!=-1)	{
	if (length(outgroup_taxa)>1)	{
		outies <- paste(outgroup_taxa,collapse="\",\"");
#			outies <- paste("v(\"",outies,sep="");
		outies <- paste("v(\"",outies,"\")",sep="");
		outies <- gsub(" ","_",outies);
		revbayes_source <- c(revbayes_source,paste("outgroup = clade(",outies,");",sep=""));
		} else	{
		outies <- taxon_names[as.numeric(initial_data$Outgroup)];
		outies <- gsub(" ","_",outies);
		revbayes_source <- c(revbayes_source,paste("outgroup = clade(\"",outies,"\");",sep=""));
		}
	hip_crowd <- "\"";
	hip_crowd <- paste(hip_crowd,paste(ingroup_taxa,collapse="\",\""),sep="");
	hip_crowd <- paste(hip_crowd,"\"",sep="");
	hip_crowd <- gsub(" ","_",hip_crowd);
	revbayes_source <- c(revbayes_source,paste("ingroup = clade(",hip_crowd,");",sep=""));
	} else	{
	outies <- "ENTER_AN_OUTGROUP_HERE!"
	revbayes_source <- c(revbayes_source,paste("outgroup = clade(\"",outies,"\");",sep=""));
	hipsters <- "ENTER_THE_INGROUP_HERE!"
	revbayes_source <- c(revbayes_source,paste("ingroup = clade(\"",hipsters,"\");",sep=""));
	}
if (length(initial_data$Unscored_Taxa)>0)
	revbayes_source <- c(revbayes_source,paste("unscored_taxa <- v(",paste(initial_data$Unscored_Taxa,collapse=","),");",sep=""));

revbayes_source <- c(revbayes_source,"among_char_var <- \"ENTER_THE_AMONG-CHARACTER_RATE_DISTRIBUTION_YOU_WISH_TO_USE_HERE\";\t# enter \"gamma\" or \"lognormal\"");
#revbayes_source <- c(revbayes_source,"clock_model <- \"ENTER_THE_CLOCK_MODEL_YOU_WISH_TO_USE_HERE\";\t# enter \"strict\" for strict clock, or \"lognormal\" for relaxed clock with lognormal variation; we'll add \"dirichlet\ eventually;");

revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"#                  Get basic information about the clade                   #");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"n_data_subsets <- filenames.size();");
if (fbd_parameterization_script!="")	{
#	revbayes_source <- c(revbayes_source,paste("intervals = readDataDelimitedFile(file=\"data/",fossil_interval_file,"\",header=true);",sep=""));
	revbayes_source <- c(revbayes_source,paste("taxa <- readTaxonData(file=\"",data_file_lead,fossil_interval_file,"\");",sep=""));
	revbayes_source <- c(revbayes_source,paste("n_taxa <- taxa.size();"));
	} else {
	revbayes_source <- c(revbayes_source,"dummy <- readDiscreteCharacterData(filenames[1]);");
	revbayes_source <- c(revbayes_source,"taxa <- dummy.taxa();");
	revbayes_source <- c(revbayes_source,"n_taxa <- dummy.ntaxa();");
	}
revbayes_source <- c(revbayes_source,"n_branches <- (2 * n_taxa) - 2;");
revbayes_source <- c(revbayes_source,"moves = VectorMoves();");

revbayes_source <- c(revbayes_source,"");
if (fbd_parameterization_script!="")	{
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"# Set up appropriate parameters for speciation, extinction & sampling.     #");
	revbayes_source <- c(revbayes_source,"#      We also set up the tree search here.                                #");
	revbayes_source <- c(revbayes_source,"#                                                                          #");
	revbayes_source <- c(revbayes_source,"# NOTE: This will sometimes freeze; if it does, then edit the script so    #");
	revbayes_source <- c(revbayes_source,"#      origination & extinction are set to 1.0. This usually works!        #");
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,paste("source(\"",paste(fbd_parameterization_script,sep=""),"\");",sep=""));
#	revbayes_source <- c(revbayes_source,paste("source(\"",paste(script_file_lead,fbd_parameterization_script,sep=""),"\");",sep=""));
	} else	{
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"# Set up tree-search moves");
	revbayes_source <- c(revbayes_source,"############################################################################");
	revbayes_source <- c(revbayes_source,"topology ~ dnUniformTopology(taxa,outgroup);");
	revbayes_source <- c(revbayes_source,"moves.append = mvNNI(topology, weight=1.0);   # nearest neighbor interchange");
	revbayes_source <- c(revbayes_source,"moves.append = mvSPR(topology, weight=1.0);   # subtree pruning");
	revbayes_source <- c(revbayes_source,"for (b in 1:n_branches) {");
	revbayes_source <- c(revbayes_source,"    bl[b] ~ dnExponential(10.0);");
	revbayes_source <- c(revbayes_source,"    moves.append = mvScale(bl[b]);");
 	revbayes_source <- c(revbayes_source,"   }");
	revbayes_source <- c(revbayes_source,"tau := treeAssembly(topology, bl);  # assign branch lengths to trees");
	revbayes_source <- c(revbayes_source,"");
	}

revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Set up appropriate Q-matrices for the partitions");
revbayes_source <- c(revbayes_source,"#   as well as the among-character and among-branch");
revbayes_source <- c(revbayes_source,"#   rate variation models");
revbayes_source <- c(revbayes_source,"#  (Again, make sure that the directory is OK)");
revbayes_source <- c(revbayes_source,"############################################################################");
#if (length(unique(character_rate_partitions))>1)	{
	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accio_Parameters_for_Analysis_Partitioned_by_States_and_Ordering_and_Class.Rev\");",sep=""));
#	} else	{
#	revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accio_Parameters_for_Analysis_Partitioned_by_States_and_Ordering.Rev\");",sep=""));
#	}
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Wrap it all into your model");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"mymodel = model(tau,Q);\t\t# tau should have FBD & character evolution models attached to it");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Add monitors & commence MCMC'ing");
revbayes_source <- c(revbayes_source,"#  (Again, make sure that the source directory is OK)");
revbayes_source <- c(revbayes_source,"# NOTE: the program saves trees once every printgen generations; so, the");
revbayes_source <- c(revbayes_source,"#   lower the number, the more trees you save.");
revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"monitors = VectorMonitors();");
tree_file_name <- analysis_name;
if (length(unique(character_rate_partitions))>1)
	tree_file_name <- paste(tree_file_name,paste(partition_labels,collapse="_vs_"),sep="_");
mcmc_output <- paste("output/",tolower(tree_file_name),sep="");
revbayes_source <- c(revbayes_source,paste("## alter the file names to better reflect the analysis (e.g., ",tolower(analysis_name),"_strict_clock.log and ",tolower(analysis_name),"_strict_clock.trees)",sep=""));
revbayes_source <- c(revbayes_source,paste("monitors.append(mnModel(filename=\"",mcmc_output,".log\", printgen=10));",sep=""));
revbayes_source <- c(revbayes_source,paste("monitors.append(mnFile(tau, filename=\"",mcmc_output,".trees\",printgen=10,separator=TAB,tau));",sep=""));
revbayes_source <- c(revbayes_source,"monitors.append(mnScreen(printgen=500,mean_rt,alpha,speciation_rate,extinction_rate,psi,summed_gaps,num_samp_anc,origin_time));");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"    ################################################################################");
revbayes_source <- c(revbayes_source,"    # Here are some starting parameters for your MCMC analysis: but use your own!");
revbayes_source <- c(revbayes_source,"    # NOTE: as the number of moves increases, the greater the number of generations");
revbayes_source <- c(revbayes_source,"    #     we need to make a thorough search of parameter space.  So, as taxa and ");
revbayes_source <- c(revbayes_source,"    #     and complexity of character evolution models increases, the greater the ");
revbayes_source <- c(revbayes_source,"    #     number of generations you should use. ");
revbayes_source <- c(revbayes_source,"    ################################################################################");
revbayes_source <- c(revbayes_source,paste("no_runs=",no_runs,";\t\t# Number of independent MCMC analyses. (Even MCMC can get stuck in local optima!)",sep=""));
revbayes_source <- c(revbayes_source,"burnin_gens=10000;\t# Number of generations for the burnin pre-analysis (to tune parameters).");
revbayes_source <- c(revbayes_source,"tuning_int=200;\t\t# Frequency at which burnin analysis will tune parameters (in generations).");
revbayes_source <- c(revbayes_source,"running_gens=1000000;\t# Number of generations for the real analysis; the bigger the analysis, the more you usually need.");
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,"# Now, go read Anna Karenina.....");
revbayes_source <- c(revbayes_source,paste("source(\"scripts/Expecto_MCMC_with_Partitioned_Characters.Rev\");",sep=""));
revbayes_source <- c(revbayes_source,"# .......");
revbayes_source <- c(revbayes_source,"# Sigh: nobody remembers elementary train safety anymore.  Oh, your trees are done.");
revbayes_source <- c(revbayes_source,"");

revbayes_source <- c(revbayes_source,"############################################################################");
revbayes_source <- c(revbayes_source,"# Prepare MCMC output to get consensus tree(s) and the most probable tree(s)");
revbayes_source <- c(revbayes_source,"#    As always, double check the directories....");
revbayes_source <- c(revbayes_source,"############################################################################");

analyses_outputs <- paste(output_file_lead,tolower(tree_file_name),sep="");
if (no_runs>1)	{
	trees <- maj_rule <- max_prob_tree <- c();
	for (nr in 1:no_runs)	{
		trees <- c(trees,paste("\"",analyses_outputs,"_run_",nr,".trees\"",sep=""));
		maj_rule <- c(maj_rule,paste("\"",analyses_outputs,"_run_",nr,"_maj_rule.tre\"",sep=""));
		max_prob_tree <- c(max_prob_tree,paste("\"",analyses_outputs,"_run_",nr,"_simple_map.tre\"",sep=""));
		}
	trees <- paste(trees,collapse=",");
	maj_rule <- paste(maj_rule,collapse=",");
	max_prob_tree <- paste(max_prob_tree,collapse=",");
	revbayes_source <- c(revbayes_source,paste("tree_files <- v(",trees,");",sep=""));
	revbayes_source <- c(revbayes_source,paste("maj_rule_files <- v(",maj_rule,");",sep=""));
	revbayes_source <- c(revbayes_source,paste("most_probable_files <- v(",max_prob_tree,");",sep=""));
	} else	{
	trees <- paste("\"",analyses_outputs,".trees\"",sep="");
	maj_rule <- paste("\"",analyses_outputs,"_maj_rule.tre\"",sep="");
	max_prob_tree <- paste("\"",analyses_outputs,"_simple_map.tre\"",sep="");
	revbayes_source <- c(revbayes_source,paste("tree_files <- ",trees,";",sep=""));
	revbayes_source <- c(revbayes_source,paste("maj_rule_files <- ",maj_rule,";",sep=""));
	revbayes_source <- c(revbayes_source,paste("most_probable_files <- ",max_prob_tree,";",sep=""));
	}
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accio_Consensus_Tree.Rev\");",sep=""));
revbayes_source <- c(revbayes_source,"");
revbayes_source <- c(revbayes_source,paste("source(\"",script_file_lead,"Accio_Most_Probable_Tree.Rev\");",sep=""));
write(revbayes_source,file=filename);
return(revbayes_source);
}

# write out script to setup FBD analysis using origination, extinction, preservation and divergence bound dates calculated by some other routine.
scribio_fbd_portion_of_Rev_Bayes_script <- function(analysis_name,write_scripts_directory,origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_taxa,otu_names,uncoded_taxa="",script_file_lead="scripts/")	{
fbd_script <- c();
fbd_script <- c(fbd_script,"########################################################################");
fbd_script <- c(fbd_script,"# Set up appropriate parameters for speciation, extinction & sampling  #");
fbd_script <- c(fbd_script,"#   \"Seed\" numbers based on analyses of Paleobiology Database data.    #");
fbd_script <- c(fbd_script,"########################################################################");
fbd_script <- c(fbd_script,paste("# Diversification Rates based on ",control_taxon,sep=""));
diversification <- 1/((origination+extinction)/2);
fbd_script <- c(fbd_script,paste("speciation_rate ~ dnExponential(",round(diversification,3),");",sep=""));
fbd_script <- c(fbd_script,paste("extinction_rate ~ dnExponential(",round(diversification,3),");",sep=""));
fbd_script <- c(fbd_script,"# NOTE: If it gets stuck in this script, then set origination & extinction to 1.0");
fbd_script <- c(fbd_script,"moves.append(mvScale(speciation_rate, lambda=0.01, weight=5));");
fbd_script <- c(fbd_script,"moves.append(mvScale(speciation_rate, lambda=0.10, weight=3));");
fbd_script <- c(fbd_script,"moves.append(mvScale(speciation_rate, lambda=1.00, weight=1));");
fbd_script <- c(fbd_script,"moves.append(mvScale(extinction_rate, lambda=0.01, weight=5));");
fbd_script <- c(fbd_script,"moves.append(mvScale(extinction_rate, lambda=0.10, weight=3));");
fbd_script <- c(fbd_script,"moves.append(mvScale(extinction_rate, lambda=1.00, weight=1));");
fbd_script <- c(fbd_script,"diversification := speciation_rate - extinction_rate;");
fbd_script <- c(fbd_script,"turnover := extinction_rate/speciation_rate;");
fbd_script <- c(fbd_script,"");
if (sampling_unit=="rock")
	sampling_unit <- "rock units"
fbd_script <- c(fbd_script,paste("# Fossil Sampling Rates based on ",sampling_unit," occupied by ",control_taxon,sep=""));
fbd_script <- c(fbd_script,paste("psi ~ dnExponential(",round(1/psi,3),");",sep=""));
fbd_script <- c(fbd_script,"completeness := psi/(extinction_rate+psi);");
fbd_script <- c(fbd_script,"moves.append(mvScale(psi, lambda=0.01, weight=5));");
fbd_script <- c(fbd_script,"moves.append(mvScale(psi, lambda=0.10, weight=3));");
fbd_script <- c(fbd_script,"moves.append(mvScale(psi, lambda=1.00, weight=1));");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"# Proportional Taxon Sampling of Youngest Time Slice");
fbd_script <- c(fbd_script,paste("rho <- ",round(rho,3),";	# 'extant' sampling.",sep=""));
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"# Establish Basal Divergence Time");
fbd_script <- c(fbd_script,paste("origin_time ~ dnUnif(",divergence_bounds[1],", ",divergence_bounds[2],");",sep=""));
fbd_script <- c(fbd_script,"moves.append(mvSlide(origin_time, delta=0.01, weight=5));");
fbd_script <- c(fbd_script,"moves.append(mvSlide(origin_time, delta=0.10, weight=3));");
fbd_script <- c(fbd_script,"moves.append(mvSlide(origin_time, delta=1.00, weight=1));");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"fbd_dist = dnFBDRP(origin=origin_time, lambda=speciation_rate, mu=extinction_rate, psi=psi, rho=rho, taxa=taxa);");

fbd_script <- c(fbd_script,"############################################################################");
fbd_script <- c(fbd_script,"#                               Set up tree                                #");
fbd_script <- c(fbd_script,"############################################################################");
fbd_script <- c(fbd_script,"# create the vector of clade constraints");
fbd_script <- c(fbd_script,"constraints = v(ingroup);");
fbd_script <- c(fbd_script,"tau ~ dnConstrainedTopology(fbd_dist,constraints=constraints);");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"moves.append(mvFNPR(tau, weight=n_branches/2));                              # time-tree pruning & grafting");
fbd_script <- c(fbd_script,"moves.append(mvNNI(tau, weight=n_branches/2));                               # nearest-neighbor interchanges");
fbd_script <- c(fbd_script,"moves.append(mvCollapseExpandFossilBranch(tau,origin_time,weight=n_taxa/4)); # consider ancestor-descendant rather than sister species");
fbd_script <- c(fbd_script,"moves.append(mvNodeTimeSlideUniform(tau, weight=n_branches/2));              # adjust divergence times");
fbd_script <- c(fbd_script,"moves.append(mvRootTimeSlideUniform(tau, origin_time, weight=5));            # adjust basaldivergence time.");
fbd_script <- c(fbd_script,"");
fbd_script <- c(fbd_script,"num_samp_anc := tau.numSampledAncestors();");
fbd_script <- c(fbd_script,"for (bn in 1:n_branches)\t{");
fbd_script <- c(fbd_script,"\tdivergence_dates[bn]:=tau.nodeAge(bn)                   # this is when a hypothesized ancestor diverges or an OTU is first seen;");
fbd_script <- c(fbd_script,"\tbranch_lengths[bn]:=tau.branchLength(bn);               # this is branch *duration* not expected change!");
fbd_script <- c(fbd_script,"\torigin_dates[bn]:=tau.branchLength(bn)+tau.nodeAge(bn); # this is when a lineage diverged from its ancestor");
fbd_script <- c(fbd_script,"\t}");
fbd_script <- c(fbd_script,"summed_gaps := sum(branch_lengths);");
fbd_script <- c(fbd_script,"");

zombies <- "\"";
zombies <- paste(zombies,paste(extant_taxa,collapse="\",\""),sep="");
zombies <- paste(zombies,"\"",sep="");
zombies <- gsub(" ","_",zombies);
fbd_script <- c(fbd_script,paste("clade_extant = clade(",zombies,");",sep=""));
fbd_script <- c(fbd_script,"age_extant := tmrca(tau, clade_extant);\t# There is no particularly good reason to keep this!");
fbd_script <- c(fbd_script,"");
if (length(uncoded_taxa)>0)	{
	fbd_script <- c(fbd_script,paste("pruned_tau := fnPruneTree(tau, prune=v(\"",paste(gsub(" ","_",otu_names[uncoded_taxa]),collapse="\",\""),"\"))",sep=""));
#	} else	{
#	fbd_script <- c(fbd_script,paste(fbd_script,"pruned_tau := fnPruneTree(tau, prune=v(\",\"))\t\t#All taxa coded!",sep=""));
	}
output_file <- paste(write_scripts_directory,"Accio_",analysis_name,"_Range_Based_FBD_Parameterization.Rev",sep="");
write(fbd_script,output_file);
output <- list(fbd_script,paste("Accio_",analysis_name,"_Range_Based_FBD_Parameterization.Rev",sep=""));
names(output) <- c("script","filename");
#return(paste(study,"_Range_Based_FBD_Parameterization.Rev",sep=""));
return(output);
}

# script that lists "recent" (= latest contemporaneous) taxa
list_faux_extant_taxa <- function(analysis_name,write_scripts_directory,fossil_intervals)	{
extant_intervals <- subset(fossil_intervals,fossil_intervals$min==min(abs(fossil_intervals$min)));
taxon_names_for_file <- gsub(" ","_",extant_intervals$taxon);
clade_extant <- c();
clade_extant <- c(clade_extant,"###################################################################");
clade_extant <- c(clade_extant,"#    Read in the \"Recent\" taxa (i.e., latest co-extant taxa)    #");
clade_extant <- c(clade_extant,"###################################################################");
clade_extant <- c(clade_extant,"");
clade_extant_info <- "clade_extant = clade(";
for (tt in 1:nrow(extant_intervals))	{
	clade_extant_info <- paste(clade_extant_info,"\"",taxon_names_for_file[tt],sep="");
	if (tt < nrow(extant_intervals))
		clade_extant_info <- paste(clade_extant_info,"\",",sep="");
	}
clade_extant_info <- paste(clade_extant_info,"\");",sep="");
clade_extant <- c(clade_extant,clade_extant_info);
output_file <- paste(write_scripts_directory,analysis_name,"_Read_Faux_Extant.Rev",sep="");
write(clade_extant,file=output_file);
return(paste(analysis_name,"_Read_Faux_Extant.Rev",sep=""));
}

#fossil_intervals <- read.table(file=fossil_interval_file,header=T);
#fossil_intervals <- read.table(file=file.choose(),header=T);
faux_extant_taxa <- function(fossil_intervals)	{
return(gsub(" ","_",as.character(fossil_intervals$taxon[fossil_intervals$min==min(abs(fossil_intervals$min))])));
#extant_taxa <- gsub(" ","_",extant_intervals$taxon);
#return(extant_taxa);
}

#scribio_RevBayes_scripts_from_nexus_file_and_PaleoDB_download <- function(analysis_name,nexus_file_name,taxon_subset_file="",rate_partition="",trend_partition="",taxon_level,lump_subgenera=F,species_only=T,bogarted="",rock_unit_databases="",chronostratigraphic_databases="",paleodb_fixes="",control_taxon="",zone_taxa="",onset="Cambrian",end="Holocene",end_FBD="",exclude_uncertain_taxa=T,basic_environments=c("terr","marine","unknown"),sampling_unit="collections",time_scale_stratigraphic_scale="International",temporal_precision=0.1,read_data_directory="",write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",UNKNOWN=-11,INAP=-22)	{
scribio_RevBayes_scripts_from_nexus_file_and_PaleoDB_download <- function(analysis_name,taxon_subset_file=F,rate_partition="",trend_partition="",taxon_level,lump_subgenera=F,species_only=T,bogarted="",rock_unit_databases="",chronostratigraphic_databases="",paleodb_fixes="",control_taxon="",zone_taxa="",onset="Cambrian",end="Holocene",end_FBD="",exclude_uncertain_taxa=T,basic_environments=c("terr","marine","unknown"),sampling_unit="collections",time_scale_stratigraphic_scale="International",temporal_precision=0.1,write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",UNKNOWN=-11,INAP=-22)	{
#### PART 0: Commence ####
print("This program will read a Nexus file and then create scripts that RevBayes can use to conduct phylogenetic analyses,");
print("   using data downloaded from the Paleobiology Database (https://www.paleobiodb.org/) for stratigraphic data and then");
print("   refining/cleaning/updating those data with updated time scales and biozonation information.");
print("");
print("NOTE: The Paleobiology Database should always be considered a STARTING point for these data. This program will also");
print("   output summaries of the data that you should check. We encourage you to contribute improvements to these data");
print("   (including collections, occurrences and taxonomy) to the Paleobiology Database.  Improvements to the stratigraphic");
print("   database used to refine PaleoDB data should be sent to pjwagner@gmail.com");
print("");
print("The program will partition the matrix based on numbers of states and whether states are ordered or unordered. The output");
print("   RevBayes scripts create appropriate Q-matrices for each partition. The program also outputs stratigraphic information");
print("   that RevBayes uses in Fossilized Birth-Death analyses (a type of birth-death-sampling analyses). It provides initial");
print("   estimates of sampling and diversification rates that are used to seed the FBD analyses, although these are varied by");
print("   RevBayes over the MCMC searches. It also provides a \"recent\" sampling rate that reflects the most-likely per-taxon");
print("   sampling rate for the latest interval relevant to the study.  (This pretends that the latest taxa are 'recent'). ");
print("");
if (taxon_level=="genus" && !lump_subgenera)
	taxon_level <- "subgenus";
#### PART 1: GET CHARACTER DATA & BASIC TAXON INFORMATION ####
#if (taxon_subset_file && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print("Choose file giving subset of taxa that you wish to be analyzed");
	for (i in 1:100)	j <- 1;
	} else	{   
	print("Choose the nexus file that you wish to analyze: ");
	taxa_subset <- "";
	}
#if (taxon_subset_file!="" && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print(".....");
	taxon_subset_file_name <- file.choose();
	taxa_subset <- read.table(taxon_subset_file_name,header = F,stringsAsFactors=hell_no)[,1];
	verboten <- c("taxon","taxa","species","genus","otu");
	taxa_subset <- taxa_subset[!taxa_subset %in% verboten];
	print("Choose the nexus file that you wish to analyze: ");
	}

basic_data <- diffindo_character_matrix_by_state_numbers_and_other_partitions(analysis_name,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead="data/",polymorphs=T,UNKNOWN,INAP);
initial_data <- basic_data$initial_data;
matrix_file_names <- basic_data$matrix_file_names;
state_numbers <- basic_data$state_numbers;
state_ordering <- basic_data$state_ordering;
character_rate_partitions <- basic_data$rate_partitions;
character_trend_partitions <- basic_data$trend_partitions;
otu_names <- otu_names_used <- initial_data$OTUs;
chmatrix <- initial_data$Matrix;
coding_bias <- basic_data$coding_bias; 
initial_data$Outgroup <- as.numeric(initial_data$Outgroup);
if (initial_data$Outgroup[1]!=-1)	{
	outgroup_taxa <- otu_names[initial_data$Outgroup];
	} else	{
	outgroup_taxa  <- "";
	}
if (taxa_subset[1]=="")	{
	ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% initial_data$Outgroup];
	} else	{
	ingroup_taxa <- taxa_subset[(1:length(taxa_subset))[!taxa_subset %in% outgroup_taxa]];
	}
# we now have all of the information that we need for the character-based part of FBD analyses.
# However, let's see if there are any taxa that belong to the ingroup-clade that are excluded!
if (species_only)	{
	taxon_names <- ingroup_taxa;
	clade_members <- unique(sapply(taxon_names,diffindo_genus_names_from_species_names));
	} else	{
	clade_members <- ingroup_taxa;
	}

#### PART 2: LOAD EXTERNAL DATA FOR CLEANING & REFINING PALEODB DATA  ####
fossilworks_collections <- paleodb_fixes$fossilworks_collections;
paleodb_rock_reidentifications <- paleodb_fixes$paleodb_rock_reidentifications;
paleodb_collection_edits <- paleodb_fixes$paleodb_collection_edits;
if (!is.null(paleodb_collection_edits$X))
	paleodb_collection_edits$X <- NULL;
time_scale <- chronostratigraphic_databases$time_scale;
zone_database <- chronostratigraphic_databases$zones;
if (is.list(rock_unit_databases))	{
	rock_database <- rock_unit_databases$rock_unit_database;
	rock_to_zone_database <- rock_unit_databases$rock_to_zone_database;
	rock_to_zone_database$ma_lb <- temporal_precision*round(rock_to_zone_database$ma_lb/temporal_precision,0);
	rock_to_zone_database$ma_ub <- temporal_precision*round(rock_to_zone_database$ma_ub/temporal_precision,0);
	}
time_scale$ma_lb <- temporal_precision*round(time_scale$ma_lb/temporal_precision,0);
time_scale$ma_ub <- temporal_precision*round(time_scale$ma_ub/temporal_precision,0);
zone_database$ma_lb <- temporal_precision*round(zone_database$ma_lb/temporal_precision,0);
zone_database$ma_ub <- temporal_precision*round(zone_database$ma_ub/temporal_precision,0);

zone_database <- subset(zone_database,zone_database$ma_lb<=time_scale$ma_lb[match(onset,time_scale$interval)]+5);
zone_database <- subset(zone_database,zone_database$ma_ub>=time_scale$ma_ub[match(end,time_scale$interval)]-5);

#### PART 3: GET INFORMATION NEEDED TO DOWNLOAD, 'CLEAN' AND ANALYZE STRATIGRAPHIC DATA  ####
compendium <- accio_updated_taxonomy_for_analyzed_taxa(otu_names=otu_names_used,local_directory=local_directory,study=analysis_name);
if (bogarted)	{
	print("Choose the file with your private stash information: ");
	for (i in 1:100)
		j <- i;
	}
if (bogarted)	{
	bogarted_info <- file.choose();
	print("Reading your private stash now....");
	bogarted_finds <- read.csv(bogarted_info,header = T,stringsAsFactors=hell_no);
	bogarted_finds <- evanesco_na_from_matrix(bogarted_finds,replacement="");
	bogarted_finds <- subset(bogarted_finds,bogarted_finds$identified_name!="");
	if (taxon_level=="genus" || taxon_level=="subgenus")	{
		taxon_name <- bogarted_finds$identified_name;
		bogarted_finds$genus <- as.character(sapply(taxon_name,diffindo_genus_names_from_species_names));
		if (taxon_level=="subgenus")	{
			genus_name <- bogarted_finds$genus;
			subgenus_results <- sapply(genus_name,diffindo_subgenus_names_from_genus_names);
#			bogarted_finds$genus <- subgenus_results[1,];
			bogarted_finds$subgenus <- subgenus_results[2,];
			bogarted_finds$subgenus[bogarted_finds$subgenus==""] <- bogarted_finds$genus[bogarted_finds$subgenus==""];
			}
	#	add occurrences
		unique_genera <- unique(bogarted_finds$genus);
		if (!is.null(bogarted_finds$subgenus))
			unique_genera <- unique(c(bogarted_finds$genus,bogarted_finds$subgenus));
		for (u_g in 1:length(unique_genera))	{
			if (!is.na(match(unique_genera[u_g],compendium$taxon_name)))	{
				compendium$n_occs[match(unique_genera[u_g],compendium$taxon_name)] <- length(unique(bogarted_finds$collection_no[unique(which(bogarted_finds==unique_genera[u_g],arr.ind = T)[,1])]));
				}
			}
		
		}
	if (!is.null(bogarted_finds$direct_ma))	{
		bogarted_finds$direct_ma <- as.numeric(bogarted_finds$direct_ma);
		bogarted_finds$direct_ma[is.na(bogarted_finds$direct_ma)] <- 0;
		}
	if (!is.null(bogarted_finds$direct_ma_error))	{
		bogarted_finds$direct_ma_error <- as.numeric(bogarted_finds$direct_ma_error);
		bogarted_finds$direct_ma_error[is.na(bogarted_finds$direct_ma_error)] <- 0;
		}
	# add something to get genus names if this is a genus-level study.
	bogarted_finds$max_ma <- time_scale$ma_lb[match(bogarted_finds$early_interval,time_scale$interval)];
	bogarted_finds$late_interval[bogarted_finds$late_interval==""] <- bogarted_finds$early_interval[bogarted_finds$late_interval==""];
	bogarted_finds$min_ma <- time_scale$ma_ub[match(bogarted_finds$late_interval,time_scale$interval)];
	bogarted_taxa <- unique(bogarted_finds$identified_name);
	notu <- nrow(compendium);
	compendia_that_are_new <- match(bogarted_finds$identified_name,compendium$taxon_name);
	compendia_that_are_new <- compendia_that_are_new[!is.na(compendia_that_are_new)];
	if (length(compendia_that_are_new)>0)
		compendium$n_occs <- compendium$n_occs+hist(compendia_that_are_new,breaks=(0:notu),plot=F)$counts;
	}

if (sum(compendium$n_occs==0)>0)	{
	missing_taxa <- subset(compendium,compendium$n_occs==0);
	missing_taxa <- subset(missing_taxa,missing_taxa$accepted_name=="?");
	taxon_name <- missing_taxa$taxon_name;
	missing_taxa_rows <- match(taxon_name,compendium$taxon_name);
	taxon_list <- genera <- unique(sapply(taxon_name,diffindo_genus_names_from_species_names));
	if (length(taxon_list)>0)	{
		taxonomyless_finds <- accio_occurrences_for_list_of_taxa(taxon_list);
		for (mt in 1:length(missing_taxa))
			compendium$n_occs[missing_taxa_rows[mt]] <- sum(taxonomyless_finds$occurrences_compendium$identified_name==taxon_name[mt]);
		}
	### insert command for getting occurrences & collections for lists of taxa here.
	}

if (sum(compendium$n_occs==0)>0)	{
#	print(paste("The following taxa have no occurrences:",paste(compendium$taxon_name[compendium$n_occs==0],collapse=", ")));
	print("The following taxa currently have no occurrences entered into the PaleoDB:");
	print(compendium$taxon_name[compendium$n_occs==0]);
#	print(paste("The following taxa are not entered into the PaleoDB:",paste(compendium$taxon_name[compendium$taxon_no==""],collapse=", ")));
	print("");
	if (sum(compendium$taxon_no=="")>0)	{
		print("The following taxa are not entered into the PaleoDB taxonomy tables:");
		print(compendium$taxon_name[compendium$taxon_no==""]);
		}
	print("Enter Data for these into the PaleoDB and try again tomorrow or setup a separate 'bogarted' file with occurrences for these taxa!");
	print("   Make sure the file as formation, member, stage, zonation, etc., information, too. (And consider entering it into the PaleoDB later.)");
	print("Also, make sure that there are no misspellings in your nexus matrix. (Computers do not autocorrect!)");
	return();
	}
otu_names[compendium$accepted_name!="?"] <- compendium$accepted_name[compendium$accepted_name!="?"];
## get paleodb data!!!!
paleodb_data <- accio_paleodb_data_for_Rev_Bayes(otu_names,analysis_name=analysis_name,local_directory,control_taxon,zone_taxa,exclude_uncertain_taxa,taxon_level,onset,end,basic_environments,time_scale,zone_database,fossilworks_collections,paleodb_rock_reidentifications,paleodb_collection_edits,lump_subgenera,species_only);
control_collections <- paleodb_data$control_collections;
control_occurrences <- paleodb_data$control_occurrences;
control_collections$collection_no <- as.numeric(control_collections$collection_no);
control_occurrences$collection_no <- as.numeric(control_occurrences$collection_no);
control_occurrences$occurrence_no <- as.numeric(control_occurrences$occurrence_no);
if (bogarted)	{
	print("Adding your private stash to the PaleoDB data....");
	bogarted_finds$collection_no <- as.numeric(bogarted_finds$collection_no);
	bogarted_finds$paleodb_collection_no[bogarted_finds$paleodb_collection_no==""] <- 0;
	bogarted_finds$paleodb_collection_no <- as.numeric(bogarted_finds$paleodb_collection_no);
	bogarted_finds$paleodb_collection_no[bogarted_finds$paleodb_collection_no==0] <- bogarted_finds$collection_no[bogarted_finds$paleodb_collection_no==0]+ceiling(max(control_collections$collection_no)/10^(floor(log10(max(control_collections$collection_no)))-1))*10^(floor(log10(max(control_collections$collection_no)))-1);

	colnames(bogarted_finds)[match("collection_no",colnames(bogarted_finds))] <- "my_collection_no"
	colnames(bogarted_finds)[match("paleodb_collection_no",colnames(bogarted_finds))] <- "collection_no"
	column_matches <- match(colnames(bogarted_finds),colnames(control_collections))
	bogarted_coll_info_in_paleodb <- (1:ncol(bogarted_finds))[!is.na(column_matches)];
	column_matches <- column_matches[!is.na(column_matches)];
	new_paleodb_coll <- control_collections[1:length(unique(bogarted_finds$collection_no)),];
	for (nc in 1:ncol(new_paleodb_coll))	{
		if (is.numeric(new_paleodb_coll[,nc]))	{
			new_paleodb_coll[,nc] <- 0;
			} else if (is.character(new_paleodb_coll[,nc]))	{
			new_paleodb_coll[,nc] <- "";
			}
		}
#	new_paleodb_coll <- control_collections[length(unique(bogarted_finds$collection_no)),];
	new_paleodb_coll[,column_matches] <- bogarted_finds[match(unique(bogarted_finds$collection_no),bogarted_finds$collection_no),bogarted_coll_info_in_paleodb];
	control_collections <- rbind(control_collections,new_paleodb_coll);
	
	# set up occcurrences in two steps;
	# first get collection part of occurrences
	column_matches <- match(colnames(bogarted_finds),colnames(control_occurrences))
	bogarted_occr_info_in_paleodb <- (1:ncol(bogarted_finds))[!is.na(column_matches)];
	column_matches <- column_matches[!is.na(column_matches)];
	new_paleodb_occr <- control_occurrences[1:nrow(bogarted_finds),];
	for (nc in 1:ncol(new_paleodb_occr))	{
		if (is.numeric(new_paleodb_occr[,nc]))	{
			new_paleodb_occr[,nc] <- as.numeric(0);
			} else if (is.character(new_paleodb_occr[,nc]))	{
			new_paleodb_occr[,nc] <- as.character("");
			new_paleodb_occr[,nc] <- as.character(new_paleodb_occr[,nc]);
			}
		}
	new_paleodb_occr[,column_matches] <- bogarted_finds[,bogarted_occr_info_in_paleodb];
	
	# now get the taxonomy part....
	taxon <- bogarted_taxa;
	for (tt in 1:length(bogarted_taxa))	{
		if (tt==1)	{
			bogarted_taxonomy <- revelio_taxonomy_for_one_taxon(taxon=bogarted_taxa[tt],settle=T);
			} else	{
			bogarted_taxonomy <- rbind(bogarted_taxonomy,revelio_taxonomy_for_one_taxon(taxon=bogarted_taxa[tt],settle=T));
			}
		}
	bogarted_taxonomy$accepted_name[bogarted_taxonomy$taxon_name!=bogarted_taxa] <- bogarted_taxonomy$taxon_name[bogarted_taxonomy$taxon_name!=bogarted_taxa] <- bogarted_taxa[bogarted_taxonomy$taxon_name!=bogarted_taxa];
	bogarted_taxonomy$accepted_rank[match(bogarted_taxonomy$accepted_rank,taxonomic_rank)>match(taxon_level,taxonomic_rank)] <- taxon_level;
	bogarted_taxonomy <- evanesco_na_from_matrix(bogarted_taxonomy,replacement = "");
	bogarted_taxonomy$record_type <- bogarted_taxonomy$flags <- bogarted_taxonomy$early_interval <- bogarted_taxonomy$late_interval <- NULL;
	bogarted_row_to_paledob <- match(bogarted_finds$identified_name,bogarted_taxonomy$taxon_name);
	paleodb_col_to_edit <- match(colnames(bogarted_taxonomy),colnames(control_occurrences));
	paleodb_col_to_edit <- paleodb_col_to_edit[!is.na(paleodb_col_to_edit)];
	bogarted_col_w_fix <- match(colnames(control_occurrences)[paleodb_col_to_edit],colnames(bogarted_taxonomy));
	new_paleodb_occr[,paleodb_col_to_edit] <- bogarted_taxonomy[bogarted_row_to_paledob,bogarted_col_w_fix];
	new_paleodb_occr$record_type <- control_occurrences$record_type[1];
	
	#(1:nrow(control_occurrences))[is.na(as.numeric(control_occurrences$occurrence_no))]
	new_paleodb_occr$occurrence_no <- (1:nrow(new_paleodb_occr))+10^ceiling(log10(max(control_occurrences$occurrence_no)))
	control_occurrences <- rbind(control_occurrences,new_paleodb_occr);
	}

if (taxon_level=="genus" || taxon_level=="subgenus")	
	control_occurrences <- add_subgenus_names_to_paleodb_finds(paleodb_finds = control_occurrences);

if (is.data.frame(paleodb_data$zone_occurrences))	{
	zone_occurrences <- paleodb_data$zone_occurrences;
	if (taxon_level=="genus" || taxon_level=="subgenus")
		zone_occurrences <- add_subgenus_names_to_paleodb_finds(paleodb_finds = zone_occurrences);
	}
hierarchical_chronostrat <- accio_hierarchical_timescale(chronostrat_units=unique(c(unique(as.character(control_collections$early_interval)),unique(as.character(control_collections$late_interval)))),time_scale,regional_scale=time_scale_stratigraphic_scale);
hierarchical_chronostrat$ma_lb <- temporal_precision*round(hierarchical_chronostrat$ma_lb/temporal_precision,0);
hierarchical_chronostrat$ma_ub <- temporal_precision*round(hierarchical_chronostrat$ma_ub/temporal_precision,0);
finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];
if (is.na(match(end_FBD,finest_chronostrat$interval)))
	end_FBD <- finest_chronostrat$interval[max(which(finest_chronostrat==end_FBD,arr.ind = T)[,1])];
downloaded_collections <- reset_paleodb_intervals_to_desired_time_scale(collections=control_collections,finest_chronostrat = finest_chronostrat,time_scale);
downloaded_collections$min_ma <- round(temporal_precision*round(downloaded_collections$min_ma/temporal_precision,0),floor(-log10(temporal_precision)));
downloaded_collections$max_ma <- round(temporal_precision*round(downloaded_collections$max_ma/temporal_precision,0),floor(-log10(temporal_precision)));
ncolls <- nrow(downloaded_collections);

#### PART 4: REFINE CHRONOSTRATIGRAPHY OF PALEODB DATA  ####
if (!is.null(zone_occurrences))	{
	paleodb_finds <- rbind(control_occurrences,zone_occurrences);
	paleodb_finds <- paleodb_finds[order(paleodb_finds$collection_no,paleodb_finds$occurrence_no),];
	paleodb_finds <- paleodb_finds[match(unique(paleodb_finds$occurrence_no),paleodb_finds$occurrence_no),];
	} else	{
	paleodb_finds <- control_occurrences;
	}
if (is.list(rock_unit_databases))	{
	print("Refining PaleoDB data with rock-unit and biozonation databases...");
	paleodb_data_refined <- refine_collection_dates_with_external_database(study=analysis_name,collections=downloaded_collections,rock_database,zone_database,rock_to_zone_database,time_scale,directory=local_directory);
	refined_collections <- paleodb_data_refined$Recalibrated_Collections;
	chronostrat_units <- unique(c(hierarchical_chronostrat$interval,refined_collections$interval_lb,refined_collections$interval_ub));
	hierarchical_chronostrat <- accio_hierarchical_timescale(chronostrat_units,time_scale,regional_scale=time_scale_stratigraphic_scale);
	finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
	} else {
	if (is.data.frame(zone_database))	{
		print("Refining PaleoDB data with biozonation databases...");
		refined_collections <- refine_paleodb_collection_dates_with_zone_data_only(paleodb_collections=downloaded_collections,paleodb_finds,zone_database,time_scale,hierarchical_chronostrat,finest_chronostrat,examine_finds=T,temporal_precision=0.05);
		} else	{
		refined_collections <- downloaded_collections;
		refined_collections$ma_lb <- refined_collections$max_ma;
		refined_collections$ma_ub <- refined_collections$min_ma;
		refined_collections$interval_lb <- refined_collections$early_interval;
		refined_collections$interval_ub <- refined_collections$late_interval;
		}
	}

age <- temporal_precision*round(refined_collections$ma_lb/temporal_precision,0);
refined_collections$interval_lb <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "onset",fine_time_scale = finest_chronostrat);
age <- temporal_precision*round(refined_collections$ma_ub/temporal_precision,0);
refined_collections$interval_ub <- sapply(age,rebin_collection_with_time_scale,onset_or_end = "end",fine_time_scale = finest_chronostrat);

# use quantitative biostratigraphy 101 to refine dates here.
print("Using basic biostratigraphy to minimize gaps for uncertainly aged collections...");
optimized_collections <- optimo_paleodb_collection_and_occurrence_stratigraphy(paleodb_finds,paleodb_collections=refined_collections,hierarchical_chronostrat,zone_database,update_search=T);
#ddd <- (1:ncolls)[optimized_collections$ma_lb<=optimized_collections$ma_ub]
# use radiometric dates here
if (!is.null(optimized_collections$direct_ma) && sum(optimized_collections$direct_ma>0)>0)	{
	print("Using radiometric data for final ages...");
	optimized_collections <- redate_collections_with_direct_dates(collections=optimized_collections,finest_chronostrat);
	}

# get rock unit numbers if we do not have a stratigraphic database
if (is.null(optimized_collections$rock_no))	{
	print("Putting numbers on unique rock units...");
	optimized_collections <- number_unique_rock_units(paleodb_collections = optimized_collections,zone_database,time_scale);
	}

# for unentered rock units that are unique to their time and location, create dummy numbers.
optimized_collections <- name_unnamed_rock_units(paleodb_collections=optimized_collections,finest_chronostrat);

paleodb_finds <- control_occurrences;
ncolls <- nrow(optimized_collections);
finest_chronostrat$ma_lb <- temporal_precision*round(finest_chronostrat$ma_lb/temporal_precision,0)
finest_chronostrat$ma_ub <- temporal_precision*round(finest_chronostrat$ma_ub/temporal_precision,0)
paleodb_collections <- completely_rebin_collections_with_uniform_time_scale(collections=optimized_collections,uniform_time_scale = finest_chronostrat);
print(paste("Saving",paste(analysis_name,"_Refined_Collections.csv",sep=""),"..."));
write.csv(paleodb_collections,paste(local_directory,analysis_name,"_Refined_Collections.csv",sep=""),row.names=F);
#print(paste("Saving",paste(analysis_name,"_Plus_Control_Finds.csv",sep=""),"..."));
#write.csv(paleodb_finds,paste(local_directory,analysis_name,"_Plus_Control_Finds.csv",sep=""),row.names=F);

#### PART 5: GET STRATIGRAPHIC DATA THAT REVBAYES CAN USE  ####
# at this point, it becomes a little easier if we have chronostratigraphic data directly tied to finds
#paleodb_collections <- completely_rebin_collections_with_uniform_time_scale(collections = paleodb_collections,uniform_time_scale = finest_chronostrat);
if (is.null(paleodb_finds$ma_lb))	{
	paleodb_finds$ma_lb <- paleodb_collections$ma_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	paleodb_finds$ma_ub <- paleodb_collections$ma_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	}
if (is.null(paleodb_finds$bin_lb))	{
	if (is.null(paleodb_collections$bin_lb))	{
		paleodb_collections$bin_lb <- match(paleodb_collections$interval_lb,finest_chronostrat$interval);
		paleodb_collections$bin_ub <- match(paleodb_collections$interval_ub,finest_chronostrat$interval);
		}
	paleodb_finds$bin_lb <- paleodb_collections$bin_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	paleodb_finds$bin_ub <- paleodb_collections$bin_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	}
print(paste("Saving",paste(analysis_name,"_Plus_Control_Finds_recalibrated.csv",sep=""),"..."));
write.csv(paleodb_finds,paste(local_directory,analysis_name,"_Plus_Control_Finds_recalibrated.csv",sep=""),row.names=F);

## output for only the ingroup
if (taxon_level=="species")	{
	otu_finds <- paleodb_finds[paleodb_finds$accepted_name %in% otu_names_used,];
	} else if (lump_subgenera==T)	{
	otu_finds <- paleodb_finds[paleodb_finds$genus %in% otu_names_used,];
	} else	{
	taxon_names <- otu_names_used;
	genus_names <- sapply(taxon_names,diffindo_genus_names_from_species_names);
	xxx <- sapply(genus_names,diffindo_subgenus_names_from_genus_names)
	subgenera_used <- otu_names_used;
	subgenera_used[xxx[2,]!=""] <- xxx[2,xxx[2,]!=""];
	otu_finds <- paleodb_finds[paleodb_finds$subgenus %in% subgenera_used,];
	find_order <- match(otu_finds$subgenus,subgenera_used)
	otu_finds <- otu_finds[order(find_order,-otu_finds$ma_lb),];
	}
print(paste("Saving",paste(analysis_name,"_Finds_recalibrated.csv",sep=""),"..."));
write.csv(otu_finds,paste(local_directory,analysis_name,"_Finds_recalibrated.csv",sep=""),row.names=F);

strat_for_Rev_Bayes <- accio_stratigraphic_information_for_Rev_Bayes(taxa=as.character(otu_names),paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank=taxon_level,sampling_unit,lump_cooccr=T,constrain=T,end_FBD = end_FBD,temporal_precision=temporal_precision);
end_FBD_b <- rebin_collection_with_time_scale(age=min(strat_for_Rev_Bayes$fossil_information_detailed$latest_poss_fa),onset_or_end = "end",fine_time_scale = finest_chronostrat);
if (end_FBD!="" && hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)] < match(end_FBD_b,finest_chronostrat$interval))	{
	end_FBD <- end_FBD_b;
	strat_for_Rev_Bayes <- accio_stratigraphic_information_for_Rev_Bayes(taxa=as.character(otu_names),paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank=taxon_level,sampling_unit,lump_cooccr=T,constrain=T,end_FBD,temporal_precision);
	}
per_bin_info_ingroup <- accio_per_stratigraphic_interval_sampling_information_for_Rev_Bayes(taxa=otu_names,paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank=taxon_level,sampling_unit,lump_cooccr=T,constrain=T,end_FBD = end_FBD,temporal_precision=temporal_precision);

# output data for RevBayes & for you to examine 
fossil_interval_file <- paste(tolower(analysis_name),"_fossil_intervals.tsv",sep="");
fossil_interval_file_FA <- paste(tolower(analysis_name),"_fossil_intervals_FA.tsv",sep="");
print(paste("Saving",fossil_interval_file,"for RevBayes to use and other files for you to examine...."));
fossil_intervals <- strat_for_Rev_Bayes$fossil_intervals;
fossil_intervals$taxon <- gsub(" ","_",otu_names_used);	# make sure that this file uses the same names as the original analysis!
if (taxon_subset_file)	{
	keeper_taxa <- gsub(" ","_",taxa_subset);	# make sure that this file uses the same names as the original analysis!
	keeper_rows <- match(taxa_subset,fossil_intervals$taxon);
	} else	{
	keeper_rows <- 1:nrow(fossil_intervals);
	}
write.table(fossil_intervals[keeper_rows,],file=paste(write_data_directory,fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals[keeper_rows,],file=paste(local_directory,fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);

fossil_intervals_FA <- fossil_intervals;
fossil_intervals_FA$min <- fossil_intervals_FA$max-min(fossil_intervals_FA$max);
fossil_intervals_FA$max <- fossil_intervals_FA$max-min(fossil_intervals_FA$max);
print(paste("Saving",fossil_interval_file_FA,"for RevBayes to use and other files for you to examine...."));
write.table(fossil_intervals_FA[keeper_rows,],file=paste(write_data_directory,fossil_interval_file_FA,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_FA[keeper_rows,],file=paste(local_directory,fossil_interval_file_FA,sep=""),row.names = F,sep="\t",quote = F);

fossil_intervals_fuzzy <- strat_for_Rev_Bayes$fossil_intervals_fuzzy;
fossil_intervals_fuzzy$taxon <- gsub(" ","_",otu_names_used);	# make sure that this file uses the same names as the original analysis!
fuzzy_fossil_interval_file <- paste(tolower(analysis_name),"_fossil_intervals_fuzzy.tsv",sep="");
print(paste("Saving",fuzzy_fossil_interval_file,"for RevBayes to use and other files for you to examine...."));
write.table(fossil_intervals_fuzzy[keeper_rows,],file=paste(write_data_directory,fuzzy_fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_fuzzy[keeper_rows,],file=paste(local_directory,fuzzy_fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);

detailed_information_file <- paste(tolower(analysis_name),"_detailed_fossil_information.csv",sep="");
detailed_information <- strat_for_Rev_Bayes$fossil_information_detailed;
write.csv(detailed_information[keeper_rows,],file=paste(local_directory,detailed_information_file,sep=""),row.names = F);

write.csv(per_bin_info_ingroup$finds_per_bin[keeper_rows,],file=paste(local_directory,analysis_name,"_",lettercase::str_ucfirst(sampling_unit),"_Finds_per_Bin.csv",sep=""),row.names = T);
write.csv(per_bin_info_ingroup$definite_finds_per_bin[keeper_rows,],file=paste(local_directory,analysis_name,"_",lettercase::str_ucfirst(sampling_unit),"_Definite_Finds_per_Bin.csv",sep=""),row.names = T);
write.csv(per_bin_info_ingroup$sampled_in_bin[keeper_rows,],file=paste(local_directory,analysis_name,"_",lettercase::str_ucfirst(sampling_unit),"_Sampled_in_Bin.csv",sep=""),row.names = T);

# we use the "accepted_name" field for analysis; so, if this is not a species-level analysis, rewrite the field with the appropriate genus or subgenus
if (taxon_level!="species")
	paleodb_finds$accepted_name <- as.character(paleodb_finds[,match(taxon_level,colnames(paleodb_finds))]);

#### PART 6: GET INITIAL SAMPLING ESTIMATES ####
print("Now getting initial estimates of sampling rates per million years.....");
if (taxon_level=="species") {
	taxon_list2 <- sort(unique(paleodb_finds$accepted_name));
	} else	{
	taxon_list2 <- sort(unique(paleodb_finds[,match(taxon_level,colnames(paleodb_finds))]));
	}
end_FBD_z <- finest_chronostrat$interval[nrow(finest_chronostrat)];
per_bin_info_total <- accio_per_stratigraphic_interval_sampling_information_for_Rev_Bayes(taxa=taxon_list2,paleodb_finds=paleodb_finds,paleodb_collections=paleodb_collections,hierarchical_chronostrat=hierarchical_chronostrat,taxon_rank=taxon_level,sampling_unit=sampling_unit,lump_cooccr=T,constrain=T,end_FBD=end_FBD_z,temporal_precision=temporal_precision);
sampled_in_bin_richness <- colSums(per_bin_info_total$finds_per_bin>=0.5);

taxon_ranges <- sepkoskify_paleodb_data(pbdb_finds=paleodb_finds,taxon_names=taxon_list2,interval_names=finest_chronostrat$interval);
taxon_lives <- 1+taxon_ranges$bin_ub-taxon_ranges$bin_lb;
names(taxon_lives) <- rownames(taxon_ranges);
interval_richness <- accio_synoptic_richness(taxon_ranges=cbind(taxon_ranges$bin_lb,taxon_ranges$bin_ub),interval_names=finest_chronostrat$interval);
if (sampling_unit=="collection" || sampling_unit=="collections")	{
	sample_units_per_bin <- ceiling(tally_collections_occupied_by_subinterval(taxon_collections=paleodb_collections,hierarchical_chronostrat=hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision));
	} else	{
	sample_units_per_bin <- ceiling(tally_rock_units_occupied_by_subinterval(taxon_collections=paleodb_collections,hierarchical_chronostrat=hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision));
	}
if (length(sample_units_per_bin) != length(interval_richness))	{
	sample_units_per_bin <- sample_units_per_bin[names(sample_units_per_bin) %in% names(interval_richness)];
	interval_richness <- interval_richness[names(interval_richness) %in% names(sample_units_per_bin)];
	}
# get the most likely uniform, exponential, beta & lognormal distributions for per-collection or per-rock sampling
print("Estimating best uniform, exponential, beta & lognormal sampling distributions for each interval.....");
fpb <- per_bin_info_total$finds_per_bin;
fpb <- fpb[,(1:ncol(fpb))[colnames(fpb) %in% names(interval_richness)]];
if (ncol(fpb)!=length(sample_units_per_bin))	{
	# put in fix here!
	}

sampling_distributions <- accio_sampling_distributions_for_RevBayes(finds_per_bin=fpb,sample_units_per_bin,end_FBD="");
# get the expected finds per bin given the distributions found above.
bin_spans <- (finest_chronostrat$ma_lb-finest_chronostrat$ma_ub);
names(bin_spans) <- finest_chronostrat$interval;
bin_spans <- bin_spans[names(bin_spans) %in% names(interval_richness)];
psi_bin <- per_interval_psis(sampling_distributions,sample_units_per_bin);
psi_bin_pma <- per_interval_per_ma_psis(sampling_distributions,sample_units_per_bin,bin_spans);
#psi_bin <- psi_bin*(finest_chronostrat$ma_lb-finest_chronostrat$ma_ub);

# USE THE MEDIAN OVERALL SAMPLING RATE AS THE STARTING SAMPLING RATE PER MILLION YEARS
#psi <- sum(psi_bin*bin_spans)/sum(bin_spans);

# THIS IS OUR SAMPLING PROBABILITY FOR THE LATEST SET ("RECENT") TAXA
if (end_FBD=="")	{
	youngest <- min(strat_for_Rev_Bayes$fossil_information_detailed$latest_poss_la);
	end_FBD <- rebin_collection_with_time_scale(age=youngest,onset_or_end = onset,fine_time_scale = finest_chronostrat);
	}
if (is.na(match(end_FBD,names(psi_bin))))	{
	end_FBD <- finest_chronostrat$interval[max(which(finest_chronostrat==end_FBD,arr.ind = T)[,1])];
	if (is.na(match(end_FBD,names(psi_bin))))	{
		end_FBD <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb >= min(strat_for_Rev_Bayes$fossil_information_detailed$latest_poss_fa)),];
		}
	}
last_bin <- hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)];
#psi <- median(psi_bin[1:last_bin]*bin_spans[1:last_bin]);
psi <- sum(psi_bin[1:match(end_FBD,names(psi_bin))])/sum(bin_spans[1:match(end_FBD,names(psi_bin))]);	# total median expected finds divided by total time
print(paste("The median per-ma sampling rate (psi) is: ",round(psi,4),".",sep=""));
faux_recent_bin <- hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)];
if (is.na(psi_bin[faux_recent_bin]))
	faux_recent_bin <- match(end_FBD,names(psi_bin));
rho <- Poisson_rate_to_probability((psi_bin[faux_recent_bin]*bin_spans[faux_recent_bin]));	# per-taxon sampling probability given parameters for the last interval
#fpb_frb <- sort(ceiling(fpb[fpb[,match(end_FBD,colnames(fpb))]>0,match(end_FBD,colnames(fpb))]),decreasing = T);
#rho2 <- chao2(abundance=fpb_frb);
print(paste("The ML per-taxon sampling rate for the final interval (rho) is: ",round(rho,4),".",sep=""));

#### PART 7: GET INITIAL ORIGINATION & EXTINCTION ESTIMATES ####
print("Estimating origination & extinction (given sampling) rates for each interval.....")
sampled_in_bin <- 1*per_bin_info_total$finds_per_bin>0.5;
for (sb in 1:ncol(sampled_in_bin))	{
	sampled_in_bin[sampled_in_bin[,sb]==T,sb] <- 1;
	sampled_in_bin[sampled_in_bin[,sb]==F,sb] <- 0;
	}
if (ncol(sampled_in_bin) != length(interval_richness))
	sampled_in_bin <- sampled_in_bin[,(1:ncol(sampled_in_bin))[colnames(sampled_in_bin) %in% names(interval_richness)]];
synoptic_richness <- interval_richness;
diversification <- accio_initial_diversification_rates(sampled_in_bin,synoptic_richness=interval_richness,psi_bin=psi_bin,chronostrat=finest_chronostrat);
origination <- diversification[1];
extinction <- diversification[2];
#if (extinction>origination)	extinction <- 0.99*origination;
print(paste("Median ML estimates for origination is: ",round(origination,4)," per myr and extinction is: ",round(extinction,4)," per myr.",sep=""));

# GET INITIAL BOUNDS FOR DIVERGENCE TIMES;
print("Get basic estimate of initial divergence time using Bapst's cal-3 method.....")
phi <- prob_sampling_clade_bapst(p=origination,q=extinction,r=psi);
bound_1 <- max(fossil_intervals_FA$max[keeper_rows]);
bound_2 <- sort(fossil_intervals_FA$max[keeper_rows],decreasing=T)[2];
#initial_divergence <- (simple_probability_divergence(bound_1,bound_2,phi,psi) + simple_likelihood_divergence(bound_1,bound_2,psi))/2;
initial_divergence <- simple_probability_divergence(bound_1,bound_2,phi,psi);
divergence_bounds <- c(bound_1,initial_divergence);
print(paste("Initial divergence bounds are ",round(divergence_bounds[1],2)," to ",round(divergence_bounds[2],2)," Ma before end of study.",sep=""));

#### PART 8: START WRITING SCRIPTS ####
extant_file <- list_faux_extant_taxa(analysis_name,write_scripts_directory,fossil_intervals=strat_for_Rev_Bayes$fossil_intervals[keeper_rows,]);
#extant_taxa <- faux_extant_taxa(fossil_intervals=strat_for_Rev_Bayes$fossil_intervals[keeper_rows,]);
extant_taxa <- faux_extant_taxa(fossil_intervals=fossil_intervals_FA[keeper_rows,]);
#																	 study,              write_scripts_directory,origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_file,otu_names,uncoded_taxa="",script_file_lead="script/"
fbd_parameterization <- scribio_fbd_portion_of_Rev_Bayes_script(analysis_name,write_scripts_directory,origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_taxa,otu_names,uncoded_taxa=initial_data$Unscored_Taxa,script_file_lead="scripts/");
fbd_parameterization_script <- paste("scripts/",fbd_parameterization$filename,sep="");

max_age <- max(fossil_intervals_FA$max[keeper_rows]);
#scribio_Rev_Bayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,write_scripts_directory,fbd_parameterization_script,extant_file,set_wdir,output_file_lead="output/",script_file_lead="script/",no_runs=4);
revbayes_babble <- scribio_Rev_Bayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script,character_rate_partitions,character_trend_partitions,fossil_interval_file=fossil_interval_file_FA,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=3);

write(fbd_parameterization$script,file=paste(paste(local_directory,"Accio_",analysis_name,sep=""),"_Range_Based_FBD_Parameterization.Rev",sep=""));
filename <- paste(local_directory,analysis_name,sep="");
if (length(unique(character_rate_partitions))>1)
	filename <- paste(filename,"_Rate_Partitioned",sep="");
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename <- paste(filename,"_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename <- paste(filename,"_FBD_Analysis.Rev",sep="");
	}
write(revbayes_babble,file=filename);
#revbayes_babble <- scribio_Rev_Bayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,write_scripts_directory=write_scripts_directory,fbd_parameterization_script,extant_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);
}

scribio_RevBayes_scripts_from_chosen_nexus_file_and_existing_FBD_script_and_data <- function(analysis_name,taxon_subset_file=F,rate_partition="",trend_partition="",write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",UNKNOWN=-11,INAP=-22)	{
print("This program will read a Nexus file and then create scripts that RevBayes can use to conduct phylogenetic analyses.");
print("\tThe program will prompt you for:");
print("\t\t1. The original nexus file;");
print("\t\t2. A .tsv file giving first and last appearance dates of each taxon (in Ma before the youngest taxa);");
print("\t\t3. A RevBayes script setting up the parameters for diversification & sampling;");
print("");
#if (taxon_subset_file && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print("Choose file giving subset of taxa that you wish to be analyzed");
	for (i in 1:100)	j <- 1;
	} else	{   
#	print("Choose the nexus file that you wish to analyze: ");
	print("Choose a nexus file to convert to RevBayes preferred format:")
	taxa_subset <- "";
	}
#if (taxon_subset_file!="" && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print(".....");
	taxon_subset_file_name <- file.choose();
	taxa_subset <- read.table(taxon_subset_file_name,header = T,stringsAsFactors=hell_no)[,1];
#	print("Choose the nexus file that you wish to analyze: ");
	print("Choose a nexus file to convert to RevBayes preferred format:")
	}

#basic_data <- diffindo_character_matrix_by_state_numbers_and_other_partitions(analysis_name,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead="data/",polymorphs=T,UNKNOWN,INAP);
basic_data <- diffindo_character_matrix_by_state_numbers_and_other_partitions(analysis_name,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead="data/",polymorphs=T,UNKNOWN,INAP);
initial_data <- basic_data$initial_data;
matrix_file_names <- basic_data$matrix_file_names;
state_numbers <- basic_data$state_numbers;
state_ordering <- basic_data$state_ordering;
character_rate_partitions <- basic_data$rate_partitions;
character_trend_partitions <- basic_data$trend_partitions;
otu_names <- otu_names_used <- initial_data$OTUs;
chmatrix <- initial_data$Matrix;
coding_bias <- basic_data$coding_bias; 
initial_data$Outgroup <- as.numeric(initial_data$Outgroup);
if (initial_data$Outgroup[1]!=-1)	{
	outgroup_taxa <- otu_names[initial_data$Outgroup];
	} else	{
	outgroup_taxa  <- "";
	}
if (taxa_subset[1]=="")	{
	ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% initial_data$Outgroup];
	} else	{
	ingroup_taxa <- taxa_subset[(1:length(taxa_subset))[!taxa_subset %in% outgroup_taxa]];
	}
# we now have all of the information that we need for the character-based part of FBD analyses.
# However, let's see if there are any taxa that belong to the ingroup-clade that are excluded!
print("Select .tsv file with first and last appearance dates for the FBD analysis:");
clade_members <- ingroup_taxa;

for (x in 1:10^6)	x <- x;
fossil_interval_file <- file.choose();
print("Choose the FBD parameterization script: ");
fossil_intervals <- read.table(fossil_interval_file,header=T);
max_age <- max(fossil_intervals$max);

for (x in 1:10^6)	x <- x;
fbd_parameterization_script <- file.choose();
revbayes_babble <- scribio_Rev_Bayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script,character_rate_partitions,character_trend_partitions,fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",data_file_lead="",no_runs=4);

filename <- paste(local_directory,analysis_name,sep="");
if (length(unique(character_rate_partitions))>1)
	filename <- paste(filename,"_Rate_Partitioned",sep="");
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename <- paste(filename,"_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename <- paste(filename,"_FBD_Analysis.Rev",sep="");
	}
write(revbayes_babble,file=filename);
#write(revbayes_babble,file=paste(paste(local_directory,analysis_name,sep=""),"_Partitioned_FBD_Analysis.Rev",sep=""));
#revbayes_babble <- scribio_Rev_Bayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,write_scripts_directory=write_scripts_directory,fbd_parameterization_script,extant_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);
}

scribio_RevBayes_scripts_from_nexus_file_and_existing_FBD_script_and_data <- function(analysis_name,fbd_parameterization_script,fossil_interval_file,taxon_subset_file=F,rate_partition="",trend_partition="",write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",UNKNOWN=-11,INAP=-22)	{
print("This program will read a Nexus file and then create scripts that RevBayes can use to conduct phylogenetic analyses.");
print("   If conducting FBD analyses, then it relies on the user to provide the name of an FBD parameterization script as.");
print("   well as a file giving fossil intervals. IF you do not have these yet, then you should use another routine:" );
print("\t\tscribio_RevBayes_scripts_from_nexus_file_and_PaleoDB_download()");
print("");
#if (taxon_subset_file && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print("Choose file giving subset of taxa that you wish to be analyzed");
	for (i in 1:100)	j <- 1;
	} else	{   
#	print("Choose the nexus file that you wish to analyze: ");
	print("Choose a nexus file to convert to RevBayes preferred format:")
	taxa_subset <- "";
	}
#if (taxon_subset_file!="" && tolower(taxon_subset_file)!="n")	{
if (taxon_subset_file)	{
	print(".....");
	taxon_subset_file_name <- file.choose();
	taxa_subset <- read.table(taxon_subset_file_name,header = T,stringsAsFactors=hell_no)[,1];
#	print("Choose the nexus file that you wish to analyze: ");
	print("Choose a nexus file to convert to RevBayes preferred format:")
	}

#basic_data <- diffindo_character_matrix_by_state_numbers_and_other_partitions(analysis_name,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead="data/",polymorphs=T,UNKNOWN,INAP);
basic_data <- diffindo_character_matrix_by_state_numbers_and_other_partitions(analysis_name,write_data_directory,rate_partition,trend_partition,taxa_subset,data_file_lead="data/",polymorphs=T,UNKNOWN,INAP);
print("Select the .tsv file with first and last appearance dates for the FBD analysis:");
initial_data <- basic_data$initial_data;
matrix_file_names <- basic_data$matrix_file_names;
state_numbers <- basic_data$state_numbers;
state_ordering <- basic_data$state_ordering;
character_rate_partitions <- basic_data$rate_partitions;
character_trend_partitions <- basic_data$trend_partitions;
otu_names <- otu_names_used <- initial_data$OTUs;
chmatrix <- initial_data$Matrix;
coding_bias <- basic_data$coding_bias; 
initial_data$Outgroup <- as.numeric(initial_data$Outgroup);
if (initial_data$Outgroup[1]!=-1)	{
	outgroup_taxa <- otu_names[initial_data$Outgroup];
	} else	{
	outgroup_taxa  <- "";
	}
if (taxa_subset[1]=="")	{
	ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% initial_data$Outgroup];
	} else	{
	ingroup_taxa <- taxa_subset[(1:length(taxa_subset))[!taxa_subset %in% outgroup_taxa]];
	}
# we now have all of the information that we need for the character-based part of FBD analyses.
# However, let's see if there are any taxa that belong to the ingroup-clade that are excluded!

#if (species_only)	{
#	taxon_names <- ingroup_taxa;
#	clade_members <- unique(sapply(taxon_names,diffindo_genus_names_from_species_names));
#	} else	{
clade_members <- ingroup_taxa;
#	}

fossil_intervals <- read.table(file.choose(),header=T);
max_age <- max(fossil_intervals$max);

revbayes_babble <- scribio_Rev_Bayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,coding_bias,outgroup_taxa,ingroup_taxa,max_age,write_scripts_directory,fbd_parameterization_script,character_rate_partitions,character_trend_partitions,fossil_interval_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);

filename <- paste(local_directory,analysis_name,sep="");
if (length(unique(character_rate_partitions))>1)
	filename <- paste(filename,"_Rate_Partitioned",sep="");
if (fbd_parameterization_script=="")	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_Analysis.Rev",sep="");
	filename <- paste(filename,"_Analysis.Rev",sep="");
	} else	{
#	filename <- paste(paste(write_scripts_directory,analysis_name,sep=""),"_FBD_Analysis.Rev",sep="");
	filename <- paste(filename,"_FBD_Analysis.Rev",sep="");
	}
write(revbayes_babble,file=filename);
#write(revbayes_babble,file=paste(paste(local_directory,analysis_name,sep=""),"_Partitioned_FBD_Analysis.Rev",sep=""));
#revbayes_babble <- scribio_Rev_Bayes_script_for_partitioned_character_data(analysis_name,initial_data,matrix_file_names,state_numbers,state_ordering,write_scripts_directory=write_scripts_directory,fbd_parameterization_script,extant_file,set_wdir,output_file_lead="output/",script_file_lead="scripts/",no_runs=4);
}

#### 
##### ROUTINES TO DOWNLOAD, CLEAN & ORGANIZE PALEODB DATA #######

# get basic stratigraphic information that RevBayes demands
accio_paleodb_data_for_Rev_Bayes <- function(otu_names,analysis_name,local_directory,control_taxon,zone_taxa,exclude_uncertain_taxa=T,taxon_level,onset,end,basic_environments=c("marine","unknown","terrestrial"),time_scale,zone_database,fossilworks_collections="",paleodb_rock_reidentifications="",paleodb_collection_edits="",lump_subgenera=F,species_only=T)	{
# otu_names: vector giving taxon names (matching those of the original nexus file)	
# study: name of the study (e.g., "Phacophida" or "Ordovician_Bucanids")	
# control_taxon: name of a taxonomic group that can be used as a control for sampling and diversification estimates
# zone_taxa: taxa such as conodonts, forams, trilobites, graptolites, ammonites, etc., that might co-occur with the study clade and that are used for biozonation;
# exclude_uncertain_taxa: if true (default), then questionable assignments are excluded
# onset: oldest collections & occurrences to download (should be older than the oldest members of you clade)
# end: youngest collections & occurrences to download
# basic_environment: environment type to download (defaults to all)
# fossilworks_collections: data.frame downloaded from Fossilworks.org
# paleodb_rock_reidentifications: data.frame providing corrections to the rock assignments of unedittable PaleoDB collections
# paleodb_collection_edits: data.frame providing other corrections to unedittable PaleoDB collections.

print("Getting occurrence & collection data for study taxa....")
data_compendium <- accio_occurrences_for_list_of_taxa(taxon_list=otu_names,lump_subgenera,species_only);
data_compendium$occurrences_compendium$flags <- simplify2array(data_compendium$occurrences_compendium$flags);
ingroup_finds <- evanesco_na_from_matrix(data=data_compendium$occurrences_compendium,replacement = "");
#ingroup_finds[match(991528,ingroup_finds$occurrence_no),]
if (species_only)	{
	no_species <- c();
	for (ot in 1:length(otu_names))	{
		if((sum(ingroup_finds$identified_rank[unique(which(ingroup_finds==otu_names[ot],arr.ind = T)[,1])]=="species")+sum(ingroup_finds$identified_rank[unique(which(ingroup_finds==otu_names[ot],arr.ind = T)[,1])]=="subspecies"))==0)	{
			no_species <- c(no_species,otu_names[ot]);
			}
		}
	if (length(no_species)>0)	{
		print("Attention: you requested species-level occurrences only, but the following taxa have only genus-level occurrences:");
		print(paste("    ",paste(no_species,collapse=","),sep=""));
		}
	}

if (exclude_uncertain_taxa)	{
	ingroup_finds <- subset(ingroup_finds,ingroup_finds$flags!= "uncertain genus, uncertain species");
	if (taxon_level=="genus" || taxon_level=="subgenus")	{
		ingroup_finds <- subset(ingroup_finds,ingroup_finds$flags!= "uncertain genus");
		} else if (taxon_level=="species" || taxon_level=="subspecies")	{
		ingroup_finds <- subset(ingroup_finds,ingroup_finds$flags!= "uncertain species");
		}
	}

ingroup_finds$taxon <- as.character(ingroup_finds$taxon);
nfinds <- nrow(ingroup_finds);

ingroup_collections <- data_compendium$collection_compendium[order(data_compendium$collection_compendium$collection_no),];
ingroup_collections <- ingroup_collections[match(unique(ingroup_finds$collection_no),ingroup_collections$collection_no),];
ingroup_collections <- ingroup_collections[order(ingroup_collections$collection_no),];
ncolls <- nrow(ingroup_collections);

#sum(ingroup_collections$max_ma<time_scale$ma_ub[match(end,time_scale$interval)])
#sum(ingroup_collections$min_ma>time_scale$ma_lb[match(onset,time_scale$interval)])
# sometimes it returns collections that are too young or too old because of the requested taxa; delete those collections
ingroup_collections <- subset(ingroup_collections,ingroup_collections$max_ma>=time_scale$ma_ub[match(end,time_scale$interval)]);
ingroup_collections <- subset(ingroup_collections,ingroup_collections$min_ma<=time_scale$ma_lb[match(onset,time_scale$interval)]);
ncolls <- nrow(ingroup_collections);
ingroup_finds <- ingroup_finds[(1:nrow(ingroup_finds))[ingroup_finds$collection_no %in% ingroup_collections$collection_no],]

# redate any collections that might be older or younger than the study interval: if we want only Jurassic collections, then Jurassic will be the oldest poss. early interval & youngest possible late interval
ingroup_collections$late_interval[ingroup_collections$min_ma<time_scale$ma_ub[match(end,time_scale$interval)]] <- end;
ingroup_collections$min_ma[ingroup_collections$min_ma<time_scale$ma_ub[match(end,time_scale$interval)]] <- time_scale$ma_ub[match(end,time_scale$interval)];
ingroup_collections$early_interval[ingroup_collections$max_ma>time_scale$ma_lb[match(onset,time_scale$interval)]] <- onset;
ingroup_collections$max_ma[ingroup_collections$min_ma<time_scale$ma_lb[match(onset,time_scale$interval)]] <- time_scale$ma_lb[match(onset,time_scale$interval)];

#write.csv(data_compendium$collection_compendium[order(data_compendium$collection_compendium$collection_no),],
write.csv(ingroup_collections[order(ingroup_collections$collection_no),],file=paste(local_directory,analysis_name,"_Collections.csv",sep=""),row.names=F);
#write.csv(data_compendium$collection_compendium[order(data_compendium$collection_compendium$collection_no),],
#write.csv(ingroup_finds,file="Fred_Finds.csv");
write.csv(ingroup_finds,file=paste(local_directory,analysis_name,"_Finds.csv",sep=""),row.names=F);

if (control_taxon[1]!="")	{
	print("Getting occurrence & collection data for control taxa....");
	control_data <- accio_data_for_control_groups_to_seed_FBD_analyses(control_taxon,onset,end,basic_environments,species_only = T);
	control_collections <- evanesco_na_from_matrix(control_data$control_collections,"");
	control_occurrences <- evanesco_na_from_matrix(control_data$control_occurrences,"");
	control_occurrences <- evanesco_indeterminate_species(paleodb_finds = control_occurrences);
	
	if (exclude_uncertain_taxa)	{
		control_occurrences <- subset(control_occurrences,control_occurrences$flags!= "uncertain genus, uncertain species");
		if (taxon_level=="genus" || taxon_level=="subgenus")	{
			control_occurrences <- subset(control_occurrences,control_occurrences$flags!= "uncertain genus");
			} else if (taxon_level=="species" || taxon_level=="subspecies")	{
			control_occurrences <- subset(control_occurrences,control_occurrences$flags!= "uncertain species");
			}
		}
	lost_finds <- ingroup_finds[(1:nrow(ingroup_finds))[!ingroup_finds$occurrence_no %in% control_occurrences$occurrence_no],];
	if (nrow(lost_finds)>0)	{
		lost_finds$taxon <- NULL;
		control_occurrences <- rbind(control_occurrences,lost_finds);
		control_occurrences <- control_occurrences[order(control_occurrences$collection_no,control_occurrences$occurrence_no),];
		}
	control_collections <- control_collections[control_collections$collection_no %in% control_occurrences$collection_no,];

 	mislaid_collections <- (1:ncolls)[!ingroup_collections$collection_no %in% control_collections$collection_no];
	if (length(mislaid_collections) > 0 )	{
		control_collections <- rbind(control_collections,ingroup_collections[mislaid_collections,]);
		control_collections <- control_collections[order(control_collections$collection_no),];
		ingroup_finds_temp <- ingroup_finds;
		ingroup_finds_temp$taxon <- NULL;
		control_occurrences <- rbind(control_occurrences,ingroup_finds_temp[ingroup_finds_temp$collection_no %in% ingroup_collections$collection_no[mislaid_collections],]);
		control_occurrences <- control_occurrences[order(control_occurrences$collection_no,control_occurrences$occurrence_no),];
		control_occurrences <- unique(control_occurrences);
		}
	if (length(control_taxon)>0)	{
		control_file_name <- paste(control_taxon,collapse="_&_");
		} else	{
		control_file_name <- control_taxon;
		}
	write.csv(evanesco_na_from_matrix(control_data$control_collections,""),file=paste(local_directory,control_file_name,"_Collections.csv",sep=""),row.names = F);
	write.csv(evanesco_na_from_matrix(control_data$control_occurrences,""),file=paste(local_directory,control_file_name,"_Finds.csv",sep=""),row.names = F);
	} else	{
	control_collections <- ingroup_collections;
	control_occurrences <- ingroup_finds;
	}

#if (!is.null(zone_taxa))	{
if (zone_taxa[1]!="")	{
	print("Getting occurrence data for zone taxa occupying the same collections....")
	zone_taxa_data <- accio_data_for_control_groups_to_seed_FBD_analyses(control_taxon=zone_taxa,onset,end,basic_environments,species_only = T)
	#zone_collections <- zone_taxa_data$control_collections[zone_taxa_data$control_collections$collection_no %in% control_collections$collection_no];
	zone_occurrences <- zone_taxa_data$control_occurrences[zone_taxa_data$control_occurrences$collection_no %in% control_collections$collection_no,];
	write.csv(zone_taxa_data$control_occurrences[zone_taxa_data$control_occurrences$collection_no %in% control_collections$collection_no,],
			  file=paste(local_directory,analysis_name,"_Zone_Taxa_Finds.csv",sep=""),
			  row.names = F);
	}

# this provides edits to biogeography due to old glitches.
print("Now some basic cleaning of the collections data....")
#if (!is.null(fossilworks_collections))	{
if (is.data.frame(fossilworks_collections))	{
	control_collections <- reparo_paleodb_paleogeography_with_fossilworks_data(paleodb_collections=control_collections,fossil_works_geography=fossilworks_collections);
	direct_dates <- data.frame(direct_ma=as.numeric(fossilworks_collections$direct_ma[match(control_collections$collection_no,fossilworks_collections$collection_no)]),
							   direct_ma_error=as.numeric(fossilworks_collections$direct_ma_error[match(control_collections$collection_no,fossilworks_collections$collection_no)]),
							   direct_ma_method=as.character(fossilworks_collections$direct_ma_method[match(control_collections$collection_no,fossilworks_collections$collection_no)]),
							   stringsAsFactors=hell_no);
	direct_dates <- evanesco_na_from_matrix(direct_dates,"");
	control_collections <- cbind(control_collections,direct_dates);
	}

# this provides en masse edits for rock units used in paleodb.
#if (!is.null(paleodb_rock_reidentifications))
if (is.data.frame(paleodb_rock_reidentifications))
	control_collections <- reparo_unedittable_paleodb_rock_identification(paleodb_collections=control_collections,paleodb_rock_reidentifications=paleodb_rock_reidentifications);

# this provides edits for paleodb collections that cannot currently be edited.
#if (!is.null(paleodb_collection_edits))
if (is.data.frame(paleodb_collection_edits))
	control_collections <- reparo_unedittable_paleodb_collections(paleodb_collections=control_collections,paleodb_collection_edits=paleodb_collection_edits);

# Correct the age ranges of collections using Gradstein 2012 + addenda (or another time scale)
control_collections <- redate_paleodb_collections_with_time_scale(paleodb_collections=control_collections,time_scale,zone_database);
if (zone_taxa[1]!="")	{
	output <- list(control_collections,control_occurrences,zone_occurrences);
	names(output) <- c("control_collections","control_occurrences","zone_occurrences");
	} else	{
	output <- list(control_collections,control_occurrences);
	names(output) <- c("control_collections","control_occurrences");
	}
return(output);
}

# get occurrences for a taxon from some span of time & environment
accio_occurrence_data <- function(taxa,onset="Cambrian",end="Holocene",basic_environments="terr,marine,unknown",species_only=TRUE,clean_entered_taxa=TRUE,directory="",save_files=TRUE,output_type=".csv") {
# Arguments:
# 	taxa: proper taxonomic name
# 	onset: onset geological interval from which you want new records
# 	end: end geological interval from which you want new records
# 	basic_environments: environments to download ("terr,marine,unknown" for "terrestrial, marine, unknown")
# 	species_only: if true, then eliminate Genus sp. identifications.
#	clean_entered_taxa: remove tags from entered taxa. (These always are removed from accepted names)
#	directory: where to send output files
#	save_files: if true, then output final results
# 	file_format: the end tag on the output files: '.xls' for Excel, '.txt', '.csv" for comma-delimited
taxa <- paste(taxa, collapse = ",");
if (!is.na(match("terrestrial",basic_environments)))
	basic_environments[match("terrestrial",basic_environments)] <- "terr";
basic_environments <- paste(basic_environments,collapse=",");
taxa <- gsub(" ","%20",taxa);
http <- paste("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&show=refattr,classext,rem,entname,abund&limit=all",sep = "");
fetch <- RCurl::getURL(http);
fetched <- gsub("\"","",simplify2array(strsplit(fetch,"\r\n"))[,1]);
if (!is.na(match("THIS REQUEST RETURNED NO RECORDS",fetched)))	{
	desired_finds <- c();
	} else	{
	all_finds <- utils::read.csv(text = fetch, header = FALSE, stringsAsFactors=hell_no);
	if (all_finds[1,1]=="Warning:")	{
		# this will happen only if something goes badly wrong.
		cc <- match("occurrence_no",all_finds[,1])
		kluge_mc_kluge_face <- character();
		for (mm in 1:ncol(all_finds))
			kluge_mc_kluge_face <- c(kluge_mc_kluge_face,as.character(all_finds[cc,mm]));
		colnames(all_finds) <- kluge_mc_kluge_face;
		all_finds <- all_finds[(cc+1):nrow(all_finds),];
		} else	{
		all_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=hell_no);
		}
	all_finds <- evanesco_na_from_matrix(data=all_finds,replacement = "");
	if (species_only)	{
		xxx <- (1:nrow(all_finds))[all_finds$identified_rank %in% c("species","subspecies")]
		desired_finds <- all_finds[xxx,];
#		desired_finds <- rbind(subset(all_finds,all_finds$identified_rank=="species"),subset(all_finds,all_finds$identified_rank=="subspecies"))
		}	else	{
		desired_finds <- all_finds;
		}
	if (nrow(desired_finds)==0)	{
		return(desired_finds);
		} else	{
		noccr <- nrow(desired_finds);
		taxon_name <- desired_finds$identified_name;
		flags1 <- sapply(taxon_name,revelio_uncertain_species_assignments);
		flags3 <- flags2 <- rep("",noccr);
		taxon_name <- desired_finds$identified_name[desired_finds$accepted_rank %in% c("genus","subgenus")];
		flags2[desired_finds$accepted_rank %in% c("genus","subgenus")] <- sapply(taxon_name,revelio_uncertain_genus_assignments);
		double <- (1:noccr)[flags1!=""][(1:noccr)[flags1!=""] %in% (1:noccr)[flags2!=""]];
		flags3[flags1!=""] <- flags1[flags1!=""];
		flags3[flags2!=""] <- flags2[flags2!=""];
		flags3[double] <- paste(unique(flags2[flags2!=""]),unique(flags1[flags1!=""]),sep=", ");
		desired_finds$flags <- simplify2array(flags3);
		
		
		# use flags field to note uncertain genus or species assignments.
#		desired_finds$flags <- sapply(taxon_name,identify_taxonomic_uncertainty);
#		for (tn in 1:nrow(desired_finds))	{
#			desired_finds$flags[tn] <- identify_taxonomic_uncertainty(taxon_name=desired_finds$identified_name[tn]);
#			}
		# some uncertain genus assignments will be clarified by updated generic assignments
#		uncertain_genera <- (1:noccr)[desired_finds$flags=="uncertain genus"];
#		recertain_genera <- uncertain_genera[desired_finds$difference[uncertain_genera]=="recombined as"]
#		desired_finds$flags[recertain_genera] <- "";
	
		if (clean_entered_taxa)	{
			cleaned_names <- sapply(as.character(desired_finds$identified_name),scourgify_taxon_names);
			desired_finds$identified_name <- cleaned_names;
			}	# removes tags such as "cf.", "?", "n. sp." from entered names
	
		if (species_only)	{
			entered_species <- sort(c((1:noccr)[desired_finds$accepted_rank=="species"],(1:noccr)[desired_finds$accepted_rank=="subspecies"]));
			unentered_species <- (1:noccr)[!(1:noccr) %in% entered_species];
			desired_finds$accepted_name[unentered_species] <- desired_finds$identified_name[unentered_species];
			noccr <- nrow(desired_finds);
			}
		cleaned_names <- sapply(as.character(desired_finds$accepted_name),scourgify_taxon_names);
		desired_finds$accepted_name <- cleaned_names;
	
		# make sure that type subgenera are consistently Lophospira (Lophospira)
		genus_name <- sort(unique(desired_finds$genus[desired_finds$genus!=""]));
		ngen <- length(genus_name);
		genus_subgenus <- base::t(sapply(genus_name,diffindo_subgenus_names_from_genus_names));
		subgenera <- genus_subgenus[genus_subgenus[,2]!="",];
		s_g <- nrow(subgenera);
		if (is.null(s_g))	{
			if (length(subgenera)==2)	{
				subgenera <- matrix(data=subgenera,nrow=1,ncol=2);
				s_g <- 1;
				} else	{
				subgenera <- matrix(data=0,nrow=0,ncol=2);
				s_g <- 0;
				}
			}
		type_subgenus <- subgenera[subgenera[,2]==subgenera[,1],];
		# first do non-type subgenera (e.g., Lophospira (Ruedemannia))
		nontype_subgenus <- subgenera[subgenera[,2]!=subgenera[,1],];
		nt_s_g <- nrow(nontype_subgenus);
		if (is.null(nt_s_g))	{
			if (length(nontype_subgenus)==2)	{
				nt_s_g <- 1;
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2);
				} else	{
				nontype_subgenus <- matrix(data=0,nrow=0,ncol=2);
				nt_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < nt_s_g)	{
			sg <- sg+1;
			if (!is.na(match(nontype_subgenus[sg,2],genus_subgenus[,1])))	{
	#			doubly_ranked <- genus_name[match(nontype_subgenus[sg,2],genus_subgenus[,1])]
	#			print(nontype_subgenus[sg,2]);
				doubly_ranked <- paste(nontype_subgenus[sg,1]," (",nontype_subgenus[sg,2],")",sep="");
				subgenus_finds <- (1:noccr)[desired_finds$genus==doubly_ranked];
				desired_finds$genus[subgenus_finds] <- nontype_subgenus[sg,2];
				nn <- match(nontype_subgenus[sg,2],genus_name);
				if (is.na(nn))	{
					# if name never is used alone, then add it to the genus lists
					nn <- match(nontype_subgenus[sg,2],genus_subgenus[,1]);
					genus_name <- insert_cell_into_vector_x(x=genus_name,new_value=as.character(genus_subgenus[nn,1]),cell_no=nn);
					genus_subgenus <- insert_row_into_matrix_x(x=genus_subgenus,new_row=c(genus_subgenus[nn,1],""),row_no=nn);
					}
				}
			}
	
		# now do type subgenera (e.g., Lophospira (Lophospira))
		t_s_g <- nrow(type_subgenus);
		if (is.null(t_s_g))	{
			if (length(type_subgenus)==2)	{
				t_s_g <- 1;
				type_subgenus <- matrix(data=type_subgenus,nrow=1,ncol=2)
				} else	{
				t_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < t_s_g)	{
			sg <- sg+1;
			subgenus_finds <- (1:noccr)[desired_finds$genus==type_subgenus[sg,1]];
			desired_finds$genus[subgenus_finds] <- paste(type_subgenus[sg,1]," (",type_subgenus[sg,2],")",sep="");
			}

		# now, make sure that subgenera are not listed in two different genera
		genus_name <- sort(unique(desired_finds$genus[desired_finds$genus!=""]));
		ngen <- length(genus_name);
		genus_subgenus <- base::t(sapply(genus_name,diffindo_subgenus_names_from_genus_names));
		subgenera <- genus_subgenus[genus_subgenus[,2]!="",];
		if (is.null(nrow(subgenera)))	{
			if (length(subgenera)==2)	{
				subgenera <- matrix(data=subgenera,nrow=1,ncol=2);
				} else	{
				subgenera <- matrix(data=0,nrow=1,ncol=2);
				}
			}
		nontype_subgenus <- subgenera[subgenera[,2]!=subgenera[,1],];
		if (is.null(nrow(nontype_subgenus)))	{
			if (length(nontype_subgenus)==2)	{
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2);
				} else	{
				nontype_subgenus <- matrix(0,0,2);
				}
			}
		if (nrow(nontype_subgenus)>0)
			nontype_subgenus <- nontype_subgenus[nontype_subgenus[,2]!="",];
		nt_s_g <- nrow(nontype_subgenus);
		if (is.null(nt_s_g))	{
			if (length(nontype_subgenus)==2)	{
				nt_s_g <- 1;
				nontype_subgenus <- matrix(data=nontype_subgenus,nrow=1,ncol=2)
				} else	{
				nt_s_g <- 0;
				}
			}
		sg <- 0;
		while (sg < nt_s_g)	{
			sg <- sg+1;
			if (sum(nontype_subgenus[sg,2]==nontype_subgenus[,2])>1)	{
				sgs <- (1:nt_s_g)[nontype_subgenus[,2] %in% nontype_subgenus[sg,2]]
				ssg <- 1;
				senior_entry <- paste(nontype_subgenus[sgs[1],1]," (",nontype_subgenus[sgs[1],2],")",sep="");
				while (ssg < sum(nontype_subgenus[sg,2]==nontype_subgenus[,2]))	{
					ssg <- ssg+1;
					double_entry <- paste(nontype_subgenus[sgs[ssg],1]," (",nontype_subgenus[sgs[ssg],2],")",sep="");
					desired_finds$genus[(1:noccr)[desired_finds$genus==double_entry]] <- senior_entry;
					}
				}
			}
	
		taxon_name <- desired_finds$accepted_name;
		species_epithet <- sapply(taxon_name,diffindo_species_epithets);
		desired_finds$accepted_name <- paste(desired_finds$genus,species_epithet);

		desired_finds <- evanesco_na_from_matrix(desired_finds,replacement="");
	
		if (save_files)	{
			taxa <- gsub(",","+",taxa);
			if (onset!=end)	{
				timespan <- paste(onset,"-",end,sep="");
				}	else	timespan <- onset;
			output <- paste(timespan,"_",taxa,"_Occurrences",output_type,sep="");
			if (directory!="")
				output <- paste(directory,output,sep="");
			if (output_type==".csv")	{
				output <- gsub("TRUE","",output);
				write.csv(desired_finds,file=output,row.names = FALSE);
				}	else	{
				output <- gsub("TRUE","",output);
				write.table(desired_finds,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
				}
			}
		}	
	return(desired_finds);
	}
}

# get collections for all taxa in vector taxon_list
accio_occurrences_for_list_of_taxa <- function(taxon_list,lump_subgenera=F,species_only=T)	{
# remove any funny symbols from taxon names
taxon_list <- sapply(taxon_list,scourgify_taxon_names);
ntaxa <- length(taxon_list);

#occurrences_compendium_list <- sapply(taxa,accio_occurrences_for_one_taxon,species_only);
#occurrences_compendium_list <- base::t(sapply(taxa,accio_occurrence_data,species_only=species_only,save_files=F));
occurrences_compendium <- c();
for (tx in 1:ntaxa)	{
#	taxon_finds <- occurrences_compendium_list[[tx]];
#	tx <- tx+1;
	taxa <- taxon_list[tx];
	taxon_finds <- accio_occurrence_data(taxa,species_only=species_only,save_files=F);
	if (!is.null(taxon_finds) && (nrow(taxon_finds)==0 && species_only))
		taxon_finds <- accio_occurrence_data(taxa,species_only=F,save_files=F);
	if (!lump_subgenera && !is.null(nrow(taxon_finds)))	{
		taxon_info <- accio_taxonomic_data_for_one_taxon(taxon=taxa);
		if (nrow(taxon_info)==1 || taxon_info[2,1]!="THIS REQUEST RETURNED NO RECORDS")	{
			backup_taxon_finds <- taxon_finds;
			if (taxon_info$taxon_rank=="genus")	{
				this_genus <- (1:nrow(taxon_finds))[taxon_finds$genus %in% c(taxa,paste(taxa," (",taxa,")",sep=""))];
				taxon_finds <- taxon_finds[this_genus,];
				# this kluge protects against wonky cases where the PaleoDB has conflicted information about genus/subgenus status
				if (nrow(taxon_finds)==0)	{
					taxon_finds <- backup_taxon_finds;
					genus_name <- taxon_finds$genus
					genus_subgenus <- sapply(genus_name,diffindo_subgenus_names_from_genus_names);
					taxon_finds$genus <- genus_subgenus[2,];
					this_genus <- (1:nrow(taxon_finds))[taxon_finds$genus %in% c(taxa,paste(taxa," (",taxa,")",sep=""))];
					taxon_finds <- taxon_finds[this_genus,];
					}
				}
			}
		}
	if (!is.null(nrow(taxon_finds)))	{
		taxon <- rep(taxa,nrow(taxon_finds));
		taxon_finds <- cbind(taxon=as.character(taxon),taxon_finds);
		if (tx==1)	{
			occurrences_compendium <- taxon_finds;
			} else	{
			occurrences_compendium <- rbind(occurrences_compendium,taxon_finds);
			}
		}
	}

collection_no <- sort(unique(occurrences_compendium$collection_no));
c <- 1;
collection_compendium <- c();
#while (c < length(collection_no))	{
#	c <- c+1;
for (c in 1:length(collection_no))	{
	if (is.null(collection_compendium))	{
#	if (c==1)	{
		collection_compendium <- accio_single_locality_info(collection_no=collection_no[c]);
		} else	{
		collection_compendium <- rbind(collection_compendium,accio_single_locality_info(collection_no=collection_no[c]));
		}
	}

named_rock_units <- collection_compendium$formation;
collection_compendium$formation <- sapply(named_rock_units,scourgify_rock_unit_names);
named_rock_units <- collection_compendium$member;
collection_compendium$member <- sapply(named_rock_units,scourgify_rock_unit_names);
named_rock_units <- collection_compendium$stratgroup;
collection_compendium$stratgroup <- sapply(named_rock_units,scourgify_rock_unit_names);
zone <- collection_compendium$zone[collection_compendium$zone!=""];
collection_compendium$zone[collection_compendium$zone!=""] <- sapply(zone,turgio_zone);
web_text <- collection_compendium$collection_name;
collection_compendium$collection_name <- sapply(web_text,scourgify_web_text_vanilla);
web_text <- collection_compendium$stratcomments;
collection_compendium$stratcomments <- sapply(web_text,scourgify_web_text_vanilla);

output <- list(collection_compendium,occurrences_compendium);
names(output) <- c("collection_compendium","occurrences_compendium");
return(output);
}

## Get General within Bin Occurrences
# count occurrences by subintervals, with finds attributable only to general intervals spread among subintervals
tally_presence_by_subinterval <- function(coll_no,early_interval,late_interval,hierarchical_chronostrat,lump_cooccr=T)	{
# coll_no: colletion number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# lump_cooccr: it true, the co-occurrences in the same collection are lumped into 1 find
if (lump_cooccr)	{
	orig_coll_no <- coll_no;
	coll_no <- unique(coll_no);
	early_interval <- early_interval[match(coll_no,orig_coll_no)];
	late_interval <- late_interval[match(coll_no,orig_coll_no)];
	}
ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)];
span <- 1+abs(bins_late-bins_early);
bin_finds <- vector(length=max(hierarchical_chronostrat$bin_last));
for (f in 1:ttl_finds)
	if (bins_early[f]==bins_late[f])
		bin_finds[bins_early[f]] <- 1;
return(bin_finds);
}

# count occurrences by subintervals, with finds attributable only to general intervals spread among subintervals
tally_presence_by_subinterval_minimum <- function(early_interval,late_interval,hierarchical_chronostrat,lump_cooccr=T,constrain=F)	{
# coll_no: colletion number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# lump_cooccr: it true, the co-occurrences in the same collection are lumped into 1 find

fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));

ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)];
bins_early[bins_early<fa_latest] <- fa_latest;
bins_late[bins_late>la_earliest] <- la_earliest;
span <- 1+abs(bins_late-bins_early);
bin_finds <- vector(length=max(hierarchical_chronostrat$bin_last));
for (f in 1:ttl_finds)	{
	bin_finds[bins_early[f]:bins_late[f]] <- bin_finds[bins_early[f]:bins_late[f]]+1/span[f];
	}
return(bin_finds);
}

# count occurrences by subintervals, with finds attributable only to general intervals spread among subintervals
tally_definite_occurrences_by_subinterval <- function(coll_no,early_interval,late_interval,hierarchical_chronostrat,lump_cooccr=T)	{
# coll_no: colletion number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# lump_cooccr: it true, the co-occurrences in the same collection are lumped into 1 find
if (lump_cooccr)	{
	orig_coll_no <- coll_no;
	coll_no <- unique(coll_no);
	early_interval <- early_interval[match(coll_no,orig_coll_no)];
	late_interval <- late_interval[match(coll_no,orig_coll_no)];
	}
ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)];
span <- 1+abs(bins_late-bins_early);
bin_finds <- vector(length=max(hierarchical_chronostrat$bin_last));
for (f in 1:ttl_finds)
	if (bins_early[f]==bins_late[f])
		bin_finds[bins_early[f]] <- bin_finds[bins_early[f]]+1;
return(bin_finds);
}

# count rock-units occupied by subintervals, with possible fractional counts
tally_collections_occupied_by_subinterval <- function(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=0.1)	{
# coll_no: collection number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# constrain: return only relevant time intervals; otherwise, return results for entire time scale.
finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
finest_chronostrat <- finest_chronostrat[unique(match(finest_chronostrat$bin_first,finest_chronostrat$bin_first)),];
if (constrain)	{
	early_interval <- as.character(taxon_collections$interval_lb);
	late_interval <- as.character(taxon_collections$interval_ub);
	fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
	la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));
	} else	{
	fa_latest <- min(hierarchical_chronostrat$bin_first);
	la_earliest <- max(hierarchical_chronostrat$bin_last);
	}
collection_finds <- data.frame(array(0,dim=c(nrow(taxon_collections),1+(la_earliest-fa_latest))));
colnames(collection_finds) <- finest_chronostrat$interval;
for (cn in 1:nrow(taxon_collections))	{
	relevant_collections <- taxon_collections[cn,];
	collection_finds[cn,] <- count_units_per_bin_fuzzily(relevant_collections,finest_chronostrat = finest_chronostrat,temporal_precision = temporal_precision);
	}
return(colSums(collection_finds));
}

# get collections from a certain span of time that include a particular taxon or set of taxa
accio_collection_data <- function(taxa,onset="Cambrian",end="Holocene",basic_environments="terr,marine,unknown",standardize_members=FALSE,directory="",save_files=TRUE,species_only=FALSE,output_type=".csv") {
taxa <- paste(taxa, collapse = ",");
if (!is.na(match("terrestrial",basic_environments)))
	basic_environments[match("terrestrial",basic_environments)] <- "terr";
basic_environments <- paste(basic_environments,collapse=",");
http <- paste("http://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,refattr,ent,entname,crmod",sep="");
#http <- paste("http://paleobiodb.org/data1.2/colls/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&show=loc,paleoloc,strat,stratext,refattr,entname,lith,env,crmod",sep="");
fetch <- RCurl::getURL(http);
collections <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=hell_no);
#	http <- paste("http://www.paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,",&show=full,etbasis,strat,lith,env,timebins,timecompare,ref,ent,entname,crmod",sep="")
if (species_only)	{
	http <- paste("http://www.paleobiodb.org/data1.2/occs/list.csv?base_name=",taxa,"&interval=",onset,",",end,"&envtype=",basic_environments,"&show=full",sep="");
	fetch <- RCurl::getURL(http);
	species_finds <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=hell_no);
	species_finds <- subset(species_finds,species_finds$identified_rank=="species");
	unique_collections <- unique(species_finds$collection_no);
	save_colls <- match(unique_collections,collections$collection_no);
	save_colls <- save_colls[!is.na(save_colls)];
	collections <- collections[save_colls,];
	}
ttl_coll <- nrow(collections);
# clean up rock unit names
#clean_groups <- scourgify_rock_unit_names(named_rock_unit=collections$stratgroup,delete_rock_type=TRUE)
named_rock_unit <- collections$stratgroup
clean_groups <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
collections$stratgroup <- clean_groups
named_rock_unit <- collections$formation;
clean_formations <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
collections$formation <- clean_formations;
named_rock_unit <- collections$member;
clean_members <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
collections$member <- clean_members;
collections$collection_subset <- evanesco_na_from_vector(collections$collection_subset,"")
# standardize member/formation ranks if possible
if (standardize_members)	{
	formations <- sort(unique(clean_formations))
	members <- sort(unique(clean_members))
	confusion <- sum(members %in% formations)
	if (confusion>0)	{
		member_or_formation <- members[(1:length(members))[members %in% formations]]
		for (c in 1:confusion)	{
			# Use latest opinion.  If latest opinion is "member," then reassign
			#	all collections to the latest formation/member combo
			# If the latest opinion is tied, then go with majority rule.  If that
			#	is tied, too, then just make the damned thing a formation....
			if (member_or_formation[c]!="")	{
				vote_formation <- (1:ttl_coll)[clean_formations %in% member_or_formation[c]]
				vote_member <- (1:ttl_coll)[clean_members %in% member_or_formation[c]]
				if (max(collections$ref_pubyr[vote_formation]) > max(collections$ref_pubyr[vote_member]))	{
					### elevate member to formation in appropriate collections
					collections$formation[vote_member] <- member_or_formation[c]
					collections$member[vote_member] <- ""
					}	else if (max(collections$ref_pubyr[vote_formation]) < max(collections$ref_pubyr[vote_member]))	{
#					for (cc in 1:length(vote_formation))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(collections$ref_pubyr[vote_member]),collections$ref_pubyr[vote_member])]
					collections$formation[vote_formation] <- collections$formation[latest_opinion]
					collections$member[vote_formation] <- member_or_formation[c]
#						}
					} else if (length(vote_formation) < length(vote_member))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(collections$ref_pubyr[vote_member]),collections$ref_pubyr[vote_member])]
					collections$formation[vote_formation] <- collections$formation[latest_opinion]
					collections$member[vote_formation] <- member_or_formation[c]
					} else	{
					collections$formation[vote_member] <- member_or_formation[c]
					collections$member[vote_member] <- ""
					}
				}
			}
		}
	}

## clean zones of question marks, aff.s, etc.
collections_w_zones <- (1:ttl_coll)[collections$zone!=""];
cwz <- length(collections_w_zones);
collections$zone <- as.character(collections$zone)
zone <- as.character(collections$zone[collections_w_zones]);
collections$zone[collections_w_zones] <- sapply(zone,turgio_zone);
#sort(unique(collections$zone[collections_w_zones]))

collections <- evanesco_na_from_matrix(collections,replacement = "");
collections$late_interval[collections$late_interval==""] <- collections$early_interval[collections$late_interval==""];
if (save_files)	{
	taxa <- gsub(",","+",taxa);
	if (onset!=end)	{
		timespan <- paste(onset,"-",end,sep="");
		}	else	timespan <- onset;
	output <- paste(timespan,"_",taxa,"_Collections",output_type,sep="")
	if (directory!="")
		output <- paste(directory,output,sep="");
	if (output_type==".csv")	{
		write.csv(collections,file=output,row.names = FALSE);
		}	else	{
		write.table(collections,file=output,sep = "\t",row.names = FALSE,col.names = TRUE);
		}
	}
return(collections)
}

# get collections data for one locality
accio_single_locality_info <- function(collection_no)	{
httpC <- paste("http://paleobiodb.org/data1.2/colls/list.csv?id=",collection_no,"&show=loc,paleoloc,strat,stratext,timebins,timecompare,lith,lithext,env,geo,methods,resgroup,ref,refattr,ent,entname,crmod",sep="");
accio <- RCurl::getURL(httpC)
coll_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=hell_no))
coll_info <- evanesco_na_from_vector(coll_info,"")
return(coll_info)
}

# organize paleodb data for basic diversification analyses.
accio_data_for_control_groups_to_seed_FBD_analyses <- function(control_taxon,onset="Cambrian",end="Holocene",basic_environments="terr,marine,unknown",species_only=T)	{
control_finds <- accio_occurrence_data(taxa=control_taxon,onset,end,basic_environments,species_only,clean_entered_taxa=TRUE,
					directory="",save_files=F,output_type=".csv");
control_collections <- accio_collection_data(taxa=control_taxon,onset,end,basic_environments,standardize_members=F,
					directory="",save_files=F,species_only=F,output_type=".csv");

control_collections <- control_collections[control_collections$collection_no %in% unique(control_finds$collection_no),];

output <- list(control_collections,control_finds);
names(output) <- c("control_collections","control_occurrences");
return(output);
}

# tally numbers of rock units
number_unique_rock_units <- function(paleodb_collections,zone_database,time_scale)	{
paleodb_rocks_info <- construct_stratigraphic_data_base_from_paleodb_collections(paleodb_collections=paleodb_collections,zone_database,time_scale);
paleodb_rocks <- paleodb_rocks_info$Rock_Database;
paleodb_rocks_to_zones <- paleodb_rocks_info$Rocks_to_Zones;
paleodb_rock_thesaurus <- accio_rock_unit_thesaurus_given_paleodb_collections(paleodb_rocks);
paleodb_collections <- evanesco_na_from_matrix(paleodb_collections,replacement="");
paleodb_collections$rock_no <- as.numeric(match_paleodb_collections_to_paleodb_rock_thesaurus(paleodb_collections=paleodb_collections,paleodb_rock_thesaurus));
paleodb_collections$rock_no_sr <- paleodb_rock_thesaurus$rock_no_sr[match(paleodb_collections$rock_no,paleodb_rock_thesaurus$rock_no)];
paleodb_collections$formation_no <- paleodb_rock_thesaurus$formation_no[match(paleodb_collections$rock_no,paleodb_rock_thesaurus$rock_no)];
paleodb_collections$formation_no_sr <- paleodb_rock_thesaurus$formation_no_sr[match(paleodb_collections$rock_no_sr,paleodb_rock_thesaurus$rock_no)];
paleodb_collections <- evanesco_na_from_matrix(paleodb_collections,replacement=0);
return(paleodb_collections);
}

## Convert PaleoDB data to FBD information ##
# create basic summaries of first and last appearances that RevBayes needs for FBD
accio_stratigraphic_information_for_Rev_Bayes <- function(taxa,paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank,sampling_unit="collection",faux_recent=T,lump_cooccr=T,constrain=F,end_FBD="",temporal_precision=0.1)	{
notu <- length(taxa);
n_intervals <- max(hierarchical_chronostrat$bin_last);
if (taxon_rank=="species")	{
	taxon_rank <- "accepted_name";
	} else if (taxon_rank=="subgenus")	{
	genus_name <- taxa;
	broken_up_names <- sapply(genus_name,diffindo_subgenus_names_from_genus_names);
	taxa[broken_up_names[2,]!=""] <- broken_up_names[2,broken_up_names[2,]!=""];
	}
if (sampling_unit=="collections" || sampling_unit=="localities" || sampling_unit=="locality" || sampling_unit=="site" || sampling_unit=="sites")
	sampling_unit <- "collection";
taxon_col <- match(taxon_rank,colnames(paleodb_finds));	# different taxonomic ranks are in different columns
if (constrain)	{
	earliest_bin <- min(hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)]);
	latest_bin <- max(hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)]);
	} else	{
	earliest_bin <- min(hierarchical_chronostrat$bin_first);
	latest_bin <- max(hierarchical_chronostrat$bin_last);
	}
if (end_FBD!="")	{
	latest_bin <- min(latest_bin,hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)]);
	}

if (is.null(paleodb_finds$interval_lb))	{
	paleodb_finds$interval_lb <- as.character(paleodb_collections$interval_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]);
	paleodb_finds$interval_ub <- as.character(paleodb_collections$interval_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]);
	}
if (is.null(paleodb_finds$bin_lb))	{
	paleodb_finds$bin_lb <- hierarchical_chronostrat$bin_first[match(paleodb_finds$interval_lb,hierarchical_chronostrat$interval)];
	paleodb_finds$bin_ub <- hierarchical_chronostrat$bin_last[match(paleodb_finds$interval_ub,hierarchical_chronostrat$interval)];
	}
if (is.null(paleodb_finds$ma_lb))	{
	paleodb_finds$ma_lb <- as.numeric(paleodb_collections$ma_lb[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]);
	paleodb_finds$ma_ub <- as.numeric(paleodb_collections$ma_ub[match(paleodb_finds$collection_no,paleodb_collections$collection_no)]);
	}

fossil_intervals <- data.frame(taxon=as.character(rep("",notu)),min=as.numeric(rep(0,notu)),
	max=as.numeric(rep(0,notu)),stringsAsFactors=hell_no);
fossil_intervals_fuzzy <- data.frame(taxon=as.character(rep("",notu)),
	latest_poss_la=as.numeric(rep(0,notu)),earliest_poss_la=as.numeric(rep(0,notu)),
	latest_poss_fa=as.numeric(rep(0,notu)),earliest_poss_fa=as.numeric(rep(0,notu)),stringsAsFactors=hell_no);
fossil_information_detailed <- data.frame(taxon=as.character(rep("",notu)),
	latest_poss_la=as.numeric(rep(0,notu)),earliest_poss_la=as.numeric(rep(0,notu)),
	latest_poss_fa=as.numeric(rep(0,notu)),earliest_poss_fa=as.numeric(rep(0,notu)),
	total_finds=as.numeric(rep(0,notu)),finds_per_bin=as.character(rep("",notu)),stringsAsFactors=hell_no);
finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];

for (tx in 1:notu)	{
	fossil_information_detailed$taxon[tx] <- fossil_intervals_fuzzy$taxon[tx] <- fossil_intervals$taxon[tx] <- taxa[tx];
	# get PaleoDB collection numbers
	#taxon_colls <- as.numeric(paleodb_finds$collection_no[paleodb_finds[,taxon_col]==taxa[tx]]);
	taxon_finds <- subset(paleodb_finds,paleodb_finds[,taxon_col]==taxa[tx]);
	taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain genus, uncertain species");
	
	if (nrow(taxon_finds)==0)	{
		poss_finds <- which(paleodb_finds==taxa[tx],arr.ind = T);
		if (length(poss_finds)>0)	{
			colls_w_finds <- poss_finds[,1];
			taxon_finds <- paleodb_finds[colls_w_finds,];
			}
		}
	if (taxon_level=="genus" || taxon_level=="subgenus")	{
		taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain genus");
		} else if (taxon_level=="species" || taxon_level=="subspecies")	{
		taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain species");
		}
#	taxon_finds[order(taxon_finds$ma_lb,decreasing = T),];
	taxon_finds_set <- subset(taxon_finds,taxon_finds$interval_lb==taxon_finds$interval_ub);
	taxon_finds_fzy <- subset(taxon_finds,taxon_finds$interval_lb!=taxon_finds$interval_ub);
	
	# unfixed but definitely older collections
	if (nrow(taxon_finds_set)>0)	{
		if (nrow(taxon_finds_fzy)>0 && (min(taxon_finds_fzy$bin_ub) < min(taxon_finds_set$bin_lb)))	{
	#		print(c(tx,"l"));
			taxon_finds_older <- subset(taxon_finds_fzy,taxon_finds_fzy$bin_ub < min(taxon_finds_set$bin_lb));
			taxon_finds_older$bin_lb <- taxon_finds_older$bin_ub;
			taxon_finds_older$interval_lb <- taxon_finds_older$interval_ub;
			taxon_finds_older$ma_lb <- hierarchical_chronostrat$ma_lb[match(taxon_finds_older$interval_ub,hierarchical_chronostrat$interval)];
			taxon_finds_set <- rbind(taxon_finds_set,taxon_finds_older);
			taxon_finds_set <- taxon_finds_set[order(taxon_finds_set$collection_no),];
			}
		# unfixed but definitely younger collections
		if (nrow(taxon_finds_fzy)>0 && (max(taxon_finds_fzy$bin_lb) > max(taxon_finds_set$bin_ub)))	{
	#		print(c(tx,"u"));
			taxon_finds_younger <- subset(taxon_finds_fzy,taxon_finds_fzy$bin_lb > max(taxon_finds_set$bin_ub));
			taxon_finds_younger$bin_ub <- taxon_finds_younger$bin_lb;
			taxon_finds_younger$interval_ub <- taxon_finds_younger$interval_lb;
			taxon_finds_younger$ma_ub <- hierarchical_chronostrat$ma_ub[match(taxon_finds_younger$interval_lb,hierarchical_chronostrat$interval)];
			taxon_finds_set <- rbind(taxon_finds_set,taxon_finds_younger);
			taxon_finds_set <- taxon_finds_set[order(taxon_finds_set$collection_no),];
			}
		if (lump_cooccr || sampling_unit=="rock" || taxon_rank=="species")	{
			taxon_coll_nos <- sort(unique(taxon_finds_set$collection_no));
			} else	{
			taxon_coll_nos <- sort(taxon_finds_set$collection_no);
			}
		fossil_information_detailed$earliest_poss_fa[tx] <- fossil_intervals_fuzzy$earliest_poss_fa[tx] <- fossil_intervals$max[tx] <- max(taxon_finds_set$ma_lb);
		fossil_information_detailed$latest_poss_fa[tx] <- fossil_intervals_fuzzy$latest_poss_fa[tx] <- max(taxon_finds_set$ma_ub);
		fossil_information_detailed$earliest_poss_la[tx] <- fossil_intervals_fuzzy$earliest_poss_la[tx] <- min(taxon_finds_set$ma_lb);
		fossil_information_detailed$latest_poss_la[tx] <- fossil_intervals_fuzzy$latest_poss_la[tx] <- fossil_intervals$min[tx] <- min(taxon_finds_set$ma_ub);
		} else if (nrow(taxon_finds_fzy)>0)	{
		fossil_information_detailed$earliest_poss_fa[tx] <- fossil_intervals_fuzzy$earliest_poss_fa[tx] <- fossil_intervals$max[tx] <- max(taxon_finds$ma_lb);
		fossil_information_detailed$latest_poss_fa[tx] <- fossil_intervals_fuzzy$latest_poss_fa[tx] <- max(taxon_finds$ma_ub);
		fossil_information_detailed$earliest_poss_la[tx] <- fossil_intervals_fuzzy$earliest_poss_la[tx] <- min(taxon_finds$ma_lb);
		fossil_information_detailed$latest_poss_la[tx] <- fossil_intervals_fuzzy$latest_poss_la[tx] <- fossil_intervals$min[tx] <- min(taxon_finds$ma_ub);
		}
	
	if (nrow(taxon_finds)>0)	{
		coll_no <- match(taxon_finds$collection_no,paleodb_collections$collection_no);
		if (lump_cooccr || sampling_unit=="rock")
			coll_no <- unique(coll_no);
		fossil_information_detailed$total_finds[tx] <- length(coll_no);
		taxon_collections <- paleodb_collections[coll_no,];
		if (sampling_unit=="collection")	{
#			early_interval <- as.character(paleodb_collections$interval_lb[coll_no]);
#			late_interval <- as.character(paleodb_collections$interval_ub[coll_no]);
			f_p_b <- tally_collections_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision);
			} else	{
#			taxon_collections <- name_unnamed_rock_units(paleodb_collections=taxon_collections,finest_chronostrat,constrain=constrain,temporal_precision=temporal_precision);
			taxon_collections <- name_unnamed_rock_units(paleodb_collections=taxon_collections,finest_chronostrat);
			f_p_b <- tally_rock_units_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision);
			}
		names(f_p_b) <- finest_chronostrat$interval[earliest_bin:latest_bin];
		fossil_information_detailed$finds_per_bin[tx] <- paste(sprintf("%5.2f", f_p_b[earliest_bin:latest_bin]),collapse=";");
		} else	{
		f_p_b <- rep(0,1+(latest_bin-earliest_bin));
		fossil_information_detailed$finds_per_bin[tx] <- paste(sprintf("%5.2f", f_p_b),collapse=";");
		}
#	print(paleodb_collections$ma_lb[123])
	}
#fossil_intervals[31,]
#fossil_intervals_fuzzy[31,]
if (faux_recent)	{
	ingroup_ranges <- sepkoskify_paleodb_data(paleodb_finds,taxon_names=taxa,interval_names=finest_chronostrat$interval);
	if (end_FBD=="")	{
		youngest <- min(ingroup_ranges$ma_min);
		end_FBD <- rebin_collection_with_time_scale(age=youngest,onset_or_end = onset,fine_time_scale = finest_chronostrat);
		faux_recent <- min(fossil_intervals$min);
		fossil_intervals$min[fossil_intervals$max<youngest] <- youngest-0.01;
		fossil_intervals$min <- round(fossil_intervals$min-faux_recent,3);
		fossil_intervals$max <- round(fossil_intervals$max-faux_recent,3);

		fossil_intervals_fuzzy$latest_poss_la[fossil_intervals_fuzzy$latest_poss_la<youngest] <- youngest-0.01;
		fossil_intervals_fuzzy$earliest_poss_la[fossil_intervals_fuzzy$earliest_poss_la<youngest] <- youngest-0.01;
		fossil_intervals_fuzzy$latest_poss_la <- fossil_intervals_fuzzy$latest_poss_la-faux_recent;
		fossil_intervals_fuzzy$earliest_poss_la <- fossil_intervals_fuzzy$earliest_poss_la-faux_recent;
		fossil_intervals_fuzzy$latest_poss_fa <- fossil_intervals_fuzzy$latest_poss_fa-faux_recent;
		fossil_intervals_fuzzy$earliest_poss_fa <- fossil_intervals_fuzzy$earliest_poss_fa-faux_recent;
		} else	{
		if (min(fossil_intervals_fuzzy$latest_poss_fa)<finest_chronostrat$ma_lb[latest_bin])	{
			latest_bin <- 1+finest_chronostrat$bin_first[match(rebin_collection_with_time_scale(age=min(fossil_intervals_fuzzy$latest_poss_fa),onset_or_end = "onset",fine_time_scale = finest_chronostrat),finest_chronostrat$interval)];
			if (latest_bin>max(finest_chronostrat$bin_first))	{
				latest_bin <- max(finest_chronostrat$bin_first);
				end_FBD <- finest_chronostrat$interval[latest_bin]; 
				} else {
				end_FBD <- finest_chronostrat$interval[latest_bin]; 
				}
			}
		faux_recent <- youngest <- finest_chronostrat$ma_lb[match(end_FBD,finest_chronostrat$interval)];
		fossil_intervals$min <- fossil_intervals$min-youngest;
		fossil_intervals$min[fossil_intervals$min<0] <- 0;
		fossil_intervals$max <- fossil_intervals$max-youngest;
		fossil_intervals_fuzzy$latest_poss_la[fossil_intervals_fuzzy$latest_poss_la<youngest] <- youngest;
		fossil_intervals_fuzzy$earliest_poss_la[fossil_intervals_fuzzy$earliest_poss_la<youngest] <- youngest;
		fossil_intervals_fuzzy$latest_poss_la <- fossil_intervals_fuzzy$latest_poss_la-faux_recent;
		fossil_intervals_fuzzy$earliest_poss_la <- fossil_intervals_fuzzy$earliest_poss_la-faux_recent;
		fossil_intervals_fuzzy$latest_poss_fa <- fossil_intervals_fuzzy$latest_poss_fa-faux_recent;
		fossil_intervals_fuzzy$earliest_poss_fa <- fossil_intervals_fuzzy$earliest_poss_fa-faux_recent;
		}
	}

fossil_intervals$min <- temporal_precision*round(fossil_intervals$min/temporal_precision,0);
fossil_intervals$max <- temporal_precision*round(fossil_intervals$max/temporal_precision,0);
fossil_intervals_fuzzy$latest_poss_la <- temporal_precision*round(fossil_intervals_fuzzy$latest_poss_la/temporal_precision,0);
fossil_intervals_fuzzy$latest_poss_fa <- temporal_precision*round(fossil_intervals_fuzzy$latest_poss_fa/temporal_precision,0);
fossil_intervals_fuzzy$earliest_poss_la <- temporal_precision*round(fossil_intervals_fuzzy$earliest_poss_la/temporal_precision,0);
fossil_intervals_fuzzy$earliest_poss_fa <- temporal_precision*round(fossil_intervals_fuzzy$earliest_poss_fa/temporal_precision,0);

output <- list(fossil_intervals,fossil_intervals_fuzzy,fossil_information_detailed);
names(output) <- c("fossil_intervals","fossil_intervals_fuzzy","fossil_information_detailed");
return(output);
}

# matrices of finds per bin in different ways
accio_per_stratigraphic_interval_sampling_information_for_Rev_Bayes <- function(taxa,paleodb_finds,paleodb_collections,hierarchical_chronostrat,taxon_rank,sampling_unit="collection",lump_cooccr=T,constrain=F,end_FBD="",temporal_precision=0.1)	{
# taxa: vector of taxon names
# paleodb_finds: data.frame of Paleobiology Database finds
# paleodb_collections: data.frame of collections from the PaleoDB corresponding to finds & linked by collection_no
# hierarchical_chronostrat: chronostratigraphic scheme with information aboust subintervals
# taxon_rank: rank (genus, species, etc.) at which analysis is conducted
# sampling_unit: "collection" to count collections occupied, "rock" to count rock-units occupied.
# lump_cooccr: if TRUE, the only one occurrence per colleciton
# constrain: if TRUE, then do not extend possible range of find counts beyond the latest possible FA & earliest possible LA
if (constrain)	{
	earliest_bin <- min(hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)]);
	latest_bin <- max(hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)]);
	} else	{
	earliest_bin <- min(hierarchical_chronostrat$bin_first);
	latest_bin <- max(hierarchical_chronostrat$bin_last);
	}
if (end_FBD!="")
	latest_bin <- min(latest_bin,hierarchical_chronostrat$bin_first[match(end_FBD,hierarchical_chronostrat$interval)]);
notu <- length(taxa);
n_intervals <- max(hierarchical_chronostrat$bin_last);
if (taxon_rank=="species")	{
	taxon_rank <- "accepted_name";
	} else if (taxon_rank=="subgenus")	{
	genus_name <- taxa;
	broken_up_names <- sapply(genus_name,diffindo_subgenus_names_from_genus_names);
	taxa[broken_up_names[2,]!=""] <- broken_up_names[2,broken_up_names[2,]!=""];
	}
taxon_col <- match(taxon_rank,colnames(paleodb_finds));
if (sampling_unit=="collections" || sampling_unit=="localities" || sampling_unit=="locality" || sampling_unit=="site" || sampling_unit=="sites")
	sampling_unit <- "collection";
finds_per_bin <- array(0,dim=c(notu,n_intervals));
colnames(finds_per_bin) <- hierarchical_chronostrat$interval[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last];
rownames(finds_per_bin) <- taxa;
sampled_in_bin <- definite_finds_per_bin <- finds_per_bin <- data.frame(finds_per_bin);
if (sampling_unit=="collection")	{
	sample_col_coll <- match("collection_no",colnames(paleodb_collections));
	sample_col_occr <- match("collection_no",colnames(paleodb_finds));
	} else	{
	sample_col_coll <- match("rock_no_sr",colnames(paleodb_collections));
	if (is.null(paleodb_finds$rock_no_sr))
		paleodb_finds$rock_no_sr <- paleodb_collections$rock_no_sr[match(paleodb_finds$collection_no,paleodb_collections$collection_no)];
	sample_col_occr <- match("rock_no_sr",colnames(paleodb_finds));
	}

tx <- 0;
while (tx < notu)	{			# for debugging.....
	tx <- tx+1;					# for debugging.....
#for (tx in 1:notu)	{
	# get PaleoDB collection numbers
#	print(taxa[tx]);
	# note: we want to use collections initially just because some collections within the same rock
	#	unit can have different ages;
	taxon_finds <- subset(paleodb_finds,paleodb_finds[,taxon_col]==taxa[tx]);
	taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain genus, uncertain species");
	if (taxon_level=="genus" || taxon_level=="subgenus")	{
		taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain genus");
		} else if (taxon_level=="species" || taxon_level=="subspecies")	{
		taxon_finds <- subset(taxon_finds,taxon_finds$flags != "uncertain species");
		}
	if (nrow(taxon_finds)>0)	{
		occupation <- taxon_finds$collection_no;
		if (lump_cooccr || sampling_unit!="collection")
			occupation <- unique(occupation);
#	taxon_colls <- as.numeric(paleodb_finds$collection_no[paleodb_finds[,taxon_col]==taxa[tx]]);
	# get which collection in this data set the PaleoDB colleciton is
		coll_no <- match(occupation,paleodb_collections$collection_no);
		early_interval <- as.character(paleodb_collections$interval_lb[coll_no]);
		late_interval <- as.character(paleodb_collections$interval_ub[coll_no]);
		taxon_collections <- paleodb_collections[coll_no,];
		if (sampling_unit=="collection")	{
			finds_per_bin[tx,] <- tally_collections_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision);
#			finds_per_bin[tx,] <- tally_occurrences_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat,lump_cooccr,constrain);
			definite_finds_per_bin[tx,] <- tally_definite_occurrences_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat);
			} else	{
			## START HERE!!!!
			taxon_collections <- name_unnamed_rock_units(paleodb_collections=taxon_collections,finest_chronostrat=subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last));
			finds_per_bin[tx,] <- tally_rock_units_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=temporal_precision);
			definite_finds_per_bin[tx,] <- tally_rock_units_definitely_occupied_by_subinterval(taxon_collections,hierarchical_chronostrat,constrain=F);
			}
		sampled_in_bin[tx,] <- tally_presence_by_subinterval(coll_no,early_interval,late_interval,hierarchical_chronostrat);
		}
	# latest possible first appearance would be the earliest late interval
	}
colnames(sampled_in_bin) <- colnames(definite_finds_per_bin) <- colnames(finds_per_bin) <- gsub("\\."," ",colnames(finds_per_bin));
output <- list(finds_per_bin[earliest_bin:latest_bin],definite_finds_per_bin[earliest_bin:latest_bin],sampled_in_bin[earliest_bin:latest_bin]);
names(output) <- c("finds_per_bin","definite_finds_per_bin","sampled_in_bin");
return(output);
}

# use Fossilworks information to fill in blanks in PaleoDB data
# there are some cases where the PaleoDB does not provide paleogeographic info but Fossil Works does. 
reparo_paleodb_paleogeography_with_fossilworks_data <- function(paleodb_collections,fossil_works_geography)	{
ncolls <- nrow(paleodb_collections);
l_p_1 <- (1:ncolls)[is.na(paleodb_collections$paleolat)];
l_p_2 <- (1:ncolls)[paleodb_collections$paleolat==""];
lost_paleogeography <- sort(unique(c(l_p_1,l_p_2)));
retrieved_paleogeography <- match(paleodb_collections$collection_no[lost_paleogeography],fossil_works_geography$collection_no);
lost_paleogeography <- lost_paleogeography[!is.na(retrieved_paleogeography)];
retrieved_paleogeography <- retrieved_paleogeography[!is.na(retrieved_paleogeography)];
still_broke <- sort(unique(c(l_p_1,l_p_2)))[!sort(unique(c(l_p_1,l_p_2))) %in% lost_paleogeography];
paleodb_collections$paleolat[still_broke] <- paleodb_collections$paleolng[still_broke] <- paleodb_collections$geoplate[still_broke] <- 0;
if (length(retrieved_paleogeography)>0)	{
	paleodb_collections$geoplate[lost_paleogeography] <- as.numeric(fossil_works_geography$plate[retrieved_paleogeography]);
	paleodb_collections$paleolat[lost_paleogeography] <- as.numeric(fossil_works_geography$paleolatdec[retrieved_paleogeography]);
	paleodb_collections$paleolng[lost_paleogeography] <- as.numeric(fossil_works_geography$paleolngdec[retrieved_paleogeography]);
	paleodb_collections$geoplate <- evanesco_na_from_vector(data=paleodb_collections$geoplate,replacement = 0);
	paleodb_collections$paleolat <- evanesco_na_from_vector(data=paleodb_collections$paleolat,replacement = 0);
	paleodb_collections$paleolng <- evanesco_na_from_vector(data=paleodb_collections$paleolng,replacement = 0);

	lost_geoplate <- (1:ncolls)[paleodb_collections$geoplate==0];
	}

return(paleodb_collections);
}


# this provides edits for corrections of rock units that cannot currently be edited online
reparo_unedittable_paleodb_rock_identification <- function(paleodb_collections,paleodb_rock_reidentifications)	{
ncolls <- nrow(paleodb_collections);
ttl_re_id <- c();
for (rr in 1:nrow(paleodb_rock_reidentifications))	{
	m1 <- (1:ncolls)[paleodb_collections$formation %in% paleodb_rock_reidentifications$Formation[rr]];
	m2 <- (1:ncolls)[paleodb_collections$zone %in% paleodb_rock_reidentifications$Member[rr]];
	m3 <- (1:ncolls)[paleodb_collections$zone %in% paleodb_rock_reidentifications$Zone[rr]];
	re_id <- m1[m1 %in% m2];
	re_id <- re_id[re_id %in% m3];
	paleodb_collections$formation[re_id] <- as.character(paleodb_rock_reidentifications$New_Formation[rr]);
	paleodb_collections$member[re_id] <- as.character(paleodb_rock_reidentifications$New_Member[rr]);
	paleodb_collections$zone[re_id] <- as.character(paleodb_rock_reidentifications$New_Zone[rr]);
	ttl_re_id <- c(ttl_re_id,re_id);
	}
return(paleodb_collections);
}

# edits for paleodb collections that cannot be editted online
reparo_unedittable_paleodb_collections <- function(paleodb_collections,paleodb_collection_edits)	{
editted_fields <- colnames(paleodb_collection_edits);

editted_colls <- match(paleodb_collections$collection_no,paleodb_collection_edits$collection_no);
editted_colls <- editted_colls[!is.na(editted_colls)];
colls_to_edit <- match(paleodb_collection_edits$collection_no,paleodb_collections$collection_no);
colls_to_edit <- colls_to_edit[!is.na(colls_to_edit)];

for (ef in 1:length(editted_fields))	{
	if (editted_fields[ef]!="collection_no")	{
		field_to_edit <- match(editted_fields[ef],colnames(paleodb_collections));
		paleodb_collections[colls_to_edit,field_to_edit] <- paleodb_collection_edits[editted_colls,ef];
#		cbind(paleodb_collections[colls_to_edit,field_to_edit],paleodb_collection_edits[editted_colls,2]);
		}
	}
return(paleodb_collections);
}

### clean up & organize formation, member & group data
# routine to cleanup rock unit names
scourgify_rock_unit_names <- function(named_rock_unit,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)	{
# was: clean_rock_unit_names
# named_rock_unit: string giving the name of a formation, member or group
# delete_rock_type: if true, the "Burgess Shale" becomes "Burgess" This is here because workers are
#	inconsistent about including rock-types in formation names
if (is.na(named_rock_unit))	named_rock_unit <- "";

named_rock_unit <- transmogrify_to_title_case(named_rock_unit);

if (named_rock_unit=="")	{
	return(named_rock_unit);
	} else	{
	nru <- named_rock_unit;		# for debugging
#	named_rock_unit <- str_lowercase(named_rock_unit);		# function no longer works
#	named_rock_unit <- str_ucfirst(named_rock_unit);		# function no longer works
	named_rock_unit <- gsub(" \\(\\?\\)","",named_rock_unit);
	named_rock_unit <- gsub("\\(\\?\\) ","",named_rock_unit);
	named_rock_unit <- gsub("\\(\\?\\)","",named_rock_unit);
	named_rock_unit <- gsub("\\?","",named_rock_unit);
	named_rock_unit <- gsub("\"", "",named_rock_unit);
	named_rock_unit <- gsub("\'s ", "s ",named_rock_unit);
	named_rock_unit <- gsub("’s ", "s ",named_rock_unit);
	named_rock_unit <- gsub("\'", "’",named_rock_unit);
	named_rock_unit <- gsub("\u009d","",named_rock_unit);
	named_rock_unit <- gsub("\035","",named_rock_unit);
	named_rock_unit <- gsub("\x8","",named_rock_unit);
	named_rock_unit <- gsub(" - ","-",named_rock_unit);
	named_rock_unit <- gsub(" — ","-",named_rock_unit);
	named_rock_unit <- gsub("Ste.-","Ste. ",named_rock_unit);
	named_rock_unit <- gsub("Ste-","Ste. ",named_rock_unit);
	named_rock_unit <- gsub("St.-","St. ",named_rock_unit);
	named_rock_unit <- gsub("St-","St. ",named_rock_unit);
	named_rock_unit <- gsub("lower part","",named_rock_unit);
	named_rock_unit <- gsub("middle part","",named_rock_unit);
	named_rock_unit <- gsub("upper part","",named_rock_unit);
	named_rock_unit <- gsub("-bearing","",named_rock_unit);
	named_rock_unit <- gsub("Ã„","A",named_rock_unit);
	named_rock_unit <- gsub("ã„","a",named_rock_unit);
	named_rock_unit <- gsub("Á","A",named_rock_unit);
	named_rock_unit <- gsub("á","a",named_rock_unit);
	named_rock_unit <- gsub("Ã¡","a",named_rock_unit);
	named_rock_unit <- gsub("ã¡","a",named_rock_unit);
	named_rock_unit <- gsub("Ã¤","a",named_rock_unit);
	named_rock_unit <- gsub("ã¤","a",named_rock_unit);
	named_rock_unit <- gsub("Ã¥","a",named_rock_unit);
	named_rock_unit <- gsub("ã¥","a",named_rock_unit);
	named_rock_unit <- gsub("√†","a",named_rock_unit);
	named_rock_unit <- gsub("√§","a",named_rock_unit);
	named_rock_unit <- gsub("Ã§","c",named_rock_unit);
	named_rock_unit <- gsub("ã§","c",named_rock_unit);
	named_rock_unit <- gsub("Ã©","e",named_rock_unit);
	named_rock_unit <- gsub("ã©","e",named_rock_unit);
	named_rock_unit <- gsub("Ã¨","e",named_rock_unit);
	named_rock_unit <- gsub("Ã±","n",named_rock_unit);
	named_rock_unit <- gsub("Ã–","O",named_rock_unit);
	named_rock_unit <- gsub("Ã¸","o",named_rock_unit);
	named_rock_unit <- gsub("√≥","o",named_rock_unit);
	named_rock_unit <- gsub("√∂","o",named_rock_unit);
	named_rock_unit <- gsub("Ã¶","o",named_rock_unit);
	named_rock_unit <- gsub("Ãµ","o",named_rock_unit);
	named_rock_unit <- gsub("Ã´","o",named_rock_unit);
	named_rock_unit <- gsub("ã´","o",named_rock_unit);
	named_rock_unit <- gsub("Ã¼","u",named_rock_unit);
	named_rock_unit <- gsub("ã¼","u",named_rock_unit);
	named_rock_unit <- gsub("≈´","u",named_rock_unit);
	named_rock_unit <- gsub("√Æ","i",named_rock_unit);
	named_rock_unit <- gsub("√æ","i",named_rock_unit);
	named_rock_unit <- gsub("≈†","S",named_rock_unit);
	named_rock_unit <- gsub("≈°","s",named_rock_unit);
	named_rock_unit <- gsub("√•","a",named_rock_unit);
	named_rock_unit <- gsub("&#367;","u",named_rock_unit);
	named_rock_unit <- gsub("&#945;","α",named_rock_unit);

	named_rock_unit <- gsub(" Part Of "," ",named_rock_unit);
	named_rock_unit <- transmogrify_diacritics(named_rock_unit);

	# now, remove particular words
#	named_rock_unit_2 <- gsub("-"," ",named_rock_unit);
	n_r_u <- strsplit(named_rock_unit," ")[[1]];
	rock_names <- length(n_r_u);
	if (n_r_u[rock_names]=="Part" && rock_names>1)	{
		rock_names <- rock_names-1;
		n_r_u <- n_r_u[1:rock_names]
		}
	
	thesaurus <- cbind(c("ft","ft.","mt","mt.","ste.","ste","st","st.","ls","lst","limeston","limstone","limestonee","qzt.","sh","claystones","limestones","slates","shales","siltstones","sandstones"),
	c("Fort","Fort","Mountain","Mountain","Sainte","Sainte","Saint","Saint","Limestone","Limestone","Limestone","Limestone","Limestone","Quartzite","Shale","Claystone","Limestone","Slate","Shale","Siltstone","Sandstone"));
	th_length <- nrow(thesaurus);	

	edit <- (1:rock_names)[tolower(n_r_u) %in% thesaurus[,1]];
	if (length(edit)>0)	{
		n_r_u[edit] <- thesaurus[match(tolower(n_r_u[edit]), thesaurus[,1]),2]
		if (n_r_u[1]=="Mountain")	n_r_u[1] <- "Mount";
		}

	bad_words <- c("basal","bed","beds","biofacies","biozone","couches","facies","fm.","fm","formacion","formation","horizons","horizon","layer","level","member","mbr","mb","miembro","niveau","section","series","shelly","stage","suite","subst.","subsuite","subzone","tongue","unit","units","zone","bottom","top","(lower)","(middle)","(upper)","(bottom)","(top)","regional");
	uncensored <- (1:rock_names)[!tolower(n_r_u) %in% bad_words];
	
	named_rock_unit <- paste(n_r_u[uncensored],collapse = " ");
	
	if (dehyphenate)	{
		named_rock_unit <- gsub("-"," ",named_rock_unit);
		named_rock_unit <- gsub("—"," ",named_rock_unit);
		}

	if (delete_rock_type)	{
		named_rock_unit_2 <- gsub("-"," ",named_rock_unit);
		named_rock_unit_2 <- gsub("—"," ",named_rock_unit_2);
		named_rock_unit_2 <- gsub("/"," ",named_rock_unit_2);
		named_rock_unit_2 <- gsub("&"," ",named_rock_unit_2);
		n_r_u <- strsplit(named_rock_unit_2," ")[[1]];
		rock_names <- length(n_r_u);
		bad_words_2 <- c("argillaceous","ashes","ash","calcaerous","calcaire","carbonate","chalk","cherts","chert","clay","claystone","claystones","conglomerates","conglomerate","coquina","coquinas","dolomites","dolomite","dolostones","dolostone","flags","glauconites","glauconite","glauconitics","glauconitic","gres","grauwacke","greywacke","greywackes","grits","grit","kalk","limestone","limestones","limeston","limstone","ls.","ls","lst","lst.","marlstones","marlstone","marl","marls","marly","micrites","micrite","mudstones","mudstone","ooid","ooids","phosphatics","phosphatic","phosphorite","phosphorites","qzt.","quartzite","quartzites","sandstone","sandstones","shales","schichten","schistes","shale","shaly","siltstones","siltstone","tillite","tillites","tuff","tuffs","volcanic","volcanics");
		uncensored <- (1:rock_names)[!tolower(n_r_u) %in% bad_words_2];
		if (length(uncensored) < (rock_names-1))	{
			censored <- (1:rock_names)[tolower(n_r_u) %in% bad_words_2];
			if (length(uncensored[tolower(n_r_u[uncensored]) %in% "and"])>0)	{
				et <- uncensored[tolower(n_r_u[uncensored]) %in% "and"];
				if (et > min(censored) && et < max(censored))
					uncensored <- (1:length(uncensored))[tolower(n_r_u[uncensored])!="and"];
				}
			if (length(uncensored[tolower(n_r_u[uncensored]) %in% "et"])>0)	{
				et <- uncensored[tolower(n_r_u[uncensored]) %in% "et"];
				if (et > min(censored) && et < max(censored))
					uncensored <- (1:length(uncensored))[tolower(n_r_u[uncensored])!="et"];
				}
			if (length(uncensored[tolower(n_r_u[uncensored]) %in% "or"])>0)	{
				et <- uncensored[tolower(n_r_u[uncensored]) %in% "or"];
				if (et > min(censored) && et < max(censored))
					uncensored <- (1:length(uncensored))[n_r_u[uncensored]!="or"];
				}
			if (length(censored)==1 && censored<rock_names)	{
				french <- c("du","de","a");
				if (!is.na(match(tolower(n_r_u[censored+1]),french)))	{
					uncensored <- uncensored[uncensored!=censored+1];
					}
				}
			}
		named_rock_unit <- paste(n_r_u[uncensored],collapse = " ");
#		print(named_rock_unit);
		}

	if (delete_informal)	{
		informals <- c("basal","inferieur","lower","middle","upper","uppermost","superieur","informal");
		n_r_u <- strsplit(named_rock_unit," ")[[1]];
		rock_names <- length(n_r_u);
		formals <- (1:rock_names)[!tolower(n_r_u) %in% informals];
		named_rock_unit <- paste(n_r_u[formals],collapse = " ");
		named_rock_unit <- gsub("1", "",named_rock_unit);
		named_rock_unit <- gsub("2", "",named_rock_unit);
		named_rock_unit <- gsub("3", "",named_rock_unit);
		named_rock_unit <- gsub("4", "",named_rock_unit);
		named_rock_unit <- gsub("5", "",named_rock_unit);
		named_rock_unit <- gsub("6", "",named_rock_unit);
		named_rock_unit <- gsub("7", "",named_rock_unit);
		named_rock_unit <- gsub("8", "",named_rock_unit);
		named_rock_unit <- gsub("9", "",named_rock_unit);
		named_rock_unit <- gsub("0", "",named_rock_unit);
		
		if (sum(letters %in% tolower(named_rock_unit))==1)
			named_rock_unit <- "";
		}
#	print(named_rock_unit);
	if (named_rock_unit=="")	{
		return(named_rock_unit);
		} else	{
		named_rock_unit <- gsub(" )",")",named_rock_unit);
		j <- strsplit(as.character(named_rock_unit),split="",fixed=TRUE)[[1]];
		if (length(j)>0 && j[length(j)]==" ")	{
			if (length(j)>1)	{
				while (j[length(j)]==" ")	j <- j[1:(length(j)-1)];
				rnr <- c();
				for (i in 1:length(j))	{
					rnr<- paste(rnr,j[i],sep="");
					}
				named_rock_unit <- rnr;
				} else {
				named_rock_unit <- "";
				} 
			} else if (length(j)>0 && j[1]==" ")	{
			while (j[1]==" ")	j <- j[2:length(j)];
			rnr <- c();
			for (i in 1:length(j))	{
				rnr<- paste(rnr,j[i],sep="");
				}
			named_rock_unit <- rnr;
			}
		named_rock_unit <- gsub("  "," ",named_rock_unit);
		detritus <- c(",","/","."," and ","and"," or ","or"," ");
		if (sum(detritus==named_rock_unit))
			named_rock_unit <- "";
#		print(named_rock_unit);

#		print(c(nru,named_rock_unit))		# for debugging
		return(named_rock_unit);
		}
	}
#return(named_rock_unit);
}

# count rock-units occupied by subintervals, with possible fractional counts
tally_rock_units_definitely_occupied_by_subinterval <- function(taxon_collections,hierarchical_chronostrat,constrain=F)	{
# coll_no: collection number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# constrain: return only relevant time intervals; otherwise, return results for entire time scale.
finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
if (constrain)	{
	early_interval <- as.character(taxon_collections$interval_lb);
	late_interval <- as.character(taxon_collections$interval_ub);
	fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
	la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));
	} else	{
	fa_latest <- min(hierarchical_chronostrat$bin_first);
	la_earliest <- max(hierarchical_chronostrat$bin_last);
	}
unique_rocks <- sort(unique(taxon_collections$rock_no_sr[taxon_collections$rock_no_sr>0]));
rock_finds <- data.frame(array(0,dim=c(length(unique_rocks),1+(la_earliest-fa_latest))));
rownames(rock_finds) <- unique_rocks;
colnames(rock_finds) <- finest_chronostrat$interval;
# get rock-units
# this needs to be inside time loop!
#tx_colls <- nrow(taxon_finds);
set_collections <- subset(taxon_collections,taxon_collections$interval_lb==taxon_collections$interval_ub);
unset_collections <- subset(taxon_collections,taxon_collections$interval_lb!=taxon_collections$interval_ub);
for (rn in 1:length(unique_rocks))	{
	this_rock <- subset(set_collections,set_collections$rock_no_sr==unique_rocks[rn]);
	if (nrow(this_rock)>0)	{
		this_rock_bins <- unique(hierarchical_chronostrat$bin_first[match(this_rock$interval_lb,hierarchical_chronostrat$interval)])
		rock_finds[rn,this_rock_bins] <- 1;
		}
	}
return(colSums(rock_finds[fa_latest:la_earliest]));
}

# sometimes 2+ rock units are entered; this will split them
diffindo_lumped_rock_units <- function(named_rock_unit)	{
named_rock_unit <- gsub("-"," - ",named_rock_unit);
named_rock_unit <- gsub("—"," — ",named_rock_unit);
named_rock_unit <- gsub("/"," / ",named_rock_unit);
named_rock_unit <- gsub("&"," and ",named_rock_unit);
named_rock_unit <- gsub(","," , ",named_rock_unit);
named_rock_unit <- gsub("  "," ",named_rock_unit);
splitters <- c("-","—","/","and","or",",");
n_r_u <- strsplit(named_rock_unit," ")[[1]];
#n_r_u <- strsplit(named_rock_unit,splitters)[[1]];
rock_names <- length(n_r_u);
splits <- (1:rock_names)[tolower(n_r_u) %in% splitters]
if (length(splits)>0)	{
	splits <- c(0,splits,rock_names+1);
	multi_rock_units <- c()
	for (i in 1:(length(splits)-1))	{
		multi_rock_units <- c(multi_rock_units,paste(n_r_u[(splits[i]+1):(splits[i+1]-1)],collapse=" "));
		}
	return(multi_rock_units);
	}	else	{
	return(named_rock_unit);
	}
}

deformalize_rock_unit_names <- function(rock_unit_name,informals=c("lower","middle","upper","uppermost","Lower","Middle","Upper","Uppermost"))	{
split_formation <- stringr::str_split_fixed(rock_unit_name,pattern=" ",n=2);
if (!is.na(match(tolower(split_formation[[1]]),informals)))	{
	if (length(split_formation[1,])==2)	{
		rock_unit_name <- split_formation[[2]];
		} else	{
		rock_unit_name <- ""
		}
	}
return(rock_unit_name)
}

# get rock names sans "lower" or "upper"
transmogrify_informal_formation_addendum_to_member <- function(rock_unit_name,informals=c("lower","middle","upper","Lower","Middle","Upper"))	{
split_formation <- stringr::str_split_fixed(rock_unit_name,pattern=" ",n=2);
if (!is.na(match(tolower(split_formation[[1]]),informals)))	{
	xxx <- strsplit(split_formation[[1]],"")[[1]]
	while (xxx[1]=="")	xxx <- xxx[2:length(xxx)]
	xxx[1] <- toupper(xxx[1]);
	informal <- c();
	for (i in 1:length(xxx))	informal <- paste(informal,xxx[i],sep="");
	rock_unit_name <- paste(split_formation[[2]]," (",informal,")",sep="");
	}
return(rock_unit_name)
}

# count rock-units occupied by subintervals, with possible fractional counts
tally_rock_units_occupied_by_subinterval <- function(taxon_collections,hierarchical_chronostrat,constrain=F,temporal_precision=0.1)	{
# coll_no: collection number
# early_interval: earliest possible chronostratigraphic interval for corresponding coll_no
# late_interval: lateest possible chronostratigraphic interval for corresponding coll_no
# hierarchical_chronostrat: table denoting which chronostratigraphic units are subunits of others
# constrain: return only relevant time intervals; otherwise, return results for entire time scale.
finest_chronostrat <- subset(hierarchical_chronostrat,hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last);
finest_chronostrat <- finest_chronostrat[unique(match(finest_chronostrat$bin_first,finest_chronostrat$bin_first)),];
if (constrain)	{
	early_interval <- as.character(taxon_collections$interval_lb);
	late_interval <- as.character(taxon_collections$interval_ub);
	fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(late_interval,hierarchical_chronostrat$interval)]));
	la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(early_interval,hierarchical_chronostrat$interval)]));
	} else	{
	fa_latest <- min(hierarchical_chronostrat$bin_first);
	la_earliest <- max(hierarchical_chronostrat$bin_last);
	}
taxon_collections <- name_unnamed_rock_units(paleodb_collections=taxon_collections,finest_chronostrat);
unique_rocks <- sort(unique(taxon_collections$rock_no_sr[taxon_collections$rock_no_sr>0]));
rock_finds <- data.frame(array(0,dim=c(length(unique_rocks),1+(la_earliest-fa_latest))));
colnames(rock_finds) <- finest_chronostrat$interval;
# get rock-units
rn <- 0;
while (rn < length(unique_rocks))	{
	rn <- rn+1;
	relevant_collections <- subset(taxon_collections,taxon_collections$rock_no_sr==unique_rocks[rn]);
	rock_finds[rn,] <- count_units_per_bin_fuzzily(relevant_collections,finest_chronostrat = finest_chronostrat,temporal_precision = temporal_precision);
	}
return(colSums(rock_finds))
}

# routine to put unique numbers on unnamed rock units unique to their geoplate & time
name_unnamed_rock_units <- function(paleodb_collections,finest_chronostrat)	{
# paleodb_collections <- refined_collections;
intervals <- finest_chronostrat$interval;
ncolls <- nrow(paleodb_collections);
if (is.null(paleodb_collections$bin_lb))	{
	paleodb_collections$bin_lb <- finest_chronostrat$bin_first[match(paleodb_collections$interval_lb,finest_chronostrat$interval)];
	paleodb_collections$bin_ub <- finest_chronostrat$bin_last[match(paleodb_collections$interval_ub,finest_chronostrat$interval)];
	}
for (i in min(paleodb_collections$bin_lb):max(paleodb_collections$bin_ub))	{
	interval_collections <- subset(paleodb_collections,paleodb_collections$interval_lb==intervals[i]);
	interval_collections <- subset(interval_collections,interval_collections$interval_ub==intervals[i]);
	nloc <- nrow(interval_collections);
	if (nrow(interval_collections)>0)	{
		interval_collections_named_rocks <- subset(interval_collections,interval_collections$rock_no_sr>0);
		interval_collections_unnamed_rocks <- subset(interval_collections,interval_collections$rock_no_sr==0);
		# assumed that unnamed rocks on geoplates with no other rocks are unique rock units
		rockless_plates <- unique(as.numeric(interval_collections_unnamed_rocks$geoplate));
		rockless_plates <- rockless_plates[!is.na(rockless_plates)];
		rock_ided_plates <- unique(as.numeric(interval_collections_named_rocks$geoplate));
		rock_ided_plates <- rock_ided_plates[!is.na(rock_ided_plates)];
		plate_needs_rock <- rockless_plates[!rockless_plates %in% rock_ided_plates];
		if (length(plate_needs_rock)>0)	{
			# add inferred rock units to collections
			implied_rocks <- paste("Geoplate",plate_needs_rock,intervals[i]);
			coll_w_nameable_rocks <- (1:nloc)[interval_collections$geoplate %in% rockless_plates[!rockless_plates %in% rock_ided_plates]];
			interval_collections$formation[coll_w_nameable_rocks] <- implied_rocks[match(as.numeric(interval_collections$geoplate[coll_w_nameable_rocks]),rockless_plates[!rockless_plates %in% rock_ided_plates])];
			mx_rock_no <- max(paleodb_collections$rock_no_sr);
			interval_collections$rock_no_sr[coll_w_nameable_rocks] <- mx_rock_no+match(as.numeric(interval_collections$geoplate[coll_w_nameable_rocks]),rockless_plates[!rockless_plates %in% rock_ided_plates]);
			update_these <- match(interval_collections$collection_no[coll_w_nameable_rocks],paleodb_collections$collection_no);
			paleodb_collections$formation[update_these] <- interval_collections$formation[coll_w_nameable_rocks];
			paleodb_collections$formation_sr[update_these] <- interval_collections$formation_sr[coll_w_nameable_rocks];
			paleodb_collections$rock_no_sr[update_these] <- interval_collections$rock_no_sr[coll_w_nameable_rocks];
			}
		} else	{
		mx_rock_no <- 0;
		}
	}

if (sum(paleodb_collections$rock_no_sr)==0)	{
	geoplates <- as.numeric(unique(paleodb_collections$geoplate));
	gp <- 0;
	while (gp < length(geoplates))	{
		gp <- gp+1;
		this_plate <- (1:nrow(paleodb_collections))[as.numeric(paleodb_collections$geoplate)==geoplates[gp]];
		paleodb_collections$rock_no_sr[this_plate] <- paleodb_collections$rock_no[this_plate] <- gp+mx_rock_no;
		paleodb_collections$formation_no_sr[this_plate] <- paleodb_collections$formation_no[this_plate] <- gp+mx_rock_no;
		}
	}

return(paleodb_collections);
}

# summarize rock unit data from a set of paleodb collections
accio_rock_unit_data_from_paleodb_collections <- function(paleodb_collections,standardize_members=TRUE) {
ttl_coll <- nrow(paleodb_collections);

# clean up rock unit names
#clean_groups <- scourgify_rock_unit_names(named_rock_unit=collections$stratgroup,delete_rock_type=TRUE)
named_rock_unit <- paleodb_collections$stratgroup;
clean_groups <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
paleodb_collections$stratgroup <- clean_groups;
named_rock_unit <- paleodb_collections$formation;
clean_formations <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
paleodb_collections$formation <- clean_formations;
named_rock_unit <- paleodb_collections$member;
clean_members <- sapply(as.character(named_rock_unit),scourgify_rock_unit_names,dehyphenate=FALSE,delete_rock_type=FALSE,delete_informal=FALSE)
paleodb_collections$member <- clean_members;

# standardize member/formation ranks if possible
if (standardize_members)	{
	formations <- sort(unique(clean_formations))
	members <- sort(unique(clean_members))
	confusion <- sum(members %in% formations)
	if (confusion>0)	{
		member_or_formation <- members[(1:length(members))[members %in% formations]]
		for (c in 1:confusion)	{
			# Use latest opinion.  If latest opinion is "member," then reassign
			#	all collections to the latest formation/member combo
			# If the latest opinion is tied, then go with majority rule.  If that
			#	is tied, too, then just make the damned thing a formation....
			if (member_or_formation[c]!="")	{
				vote_formation <- (1:ttl_coll)[clean_formations %in% member_or_formation[c]]
				vote_member <- (1:ttl_coll)[clean_members %in% member_or_formation[c]]
				if (max(paleodb_collections$ref_pubyr[vote_formation]) > max(paleodb_collections$ref_pubyr[vote_member]))	{;
					### elevate member to formation in appropriate collections
					paleodb_collections$formation[vote_member] <- member_or_formation[c];
					paleodb_collections$member[vote_member] <- "";
					}	else if (max(paleodb_collections$ref_pubyr[vote_formation]) < max(paleodb_collections$ref_pubyr[vote_member]))	{;
#					for (cc in 1:length(vote_formation))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(paleodb_collections$ref_pubyr[vote_member]),paleodb_collections$ref_pubyr[vote_member])];
					paleodb_collections$formation[vote_formation] <- paleodb_collections$formation[latest_opinion];
					paleodb_collections$member[vote_formation] <- member_or_formation[c];
#						}
					} else if (length(vote_formation) < length(vote_member))	{
					## get the latest opinion, and assign the rock unit as a member to that formation
					latest_opinion <- vote_member[match(max(paleodb_collections$ref_pubyr[vote_member]),paleodb_collections$ref_pubyr[vote_member])];
					paleodb_collections$formation[vote_formation] <- paleodb_collections$formation[latest_opinion];
					paleodb_collections$member[vote_formation] <- member_or_formation[c];
					} else	{
					paleodb_collections$formation[vote_member] <- member_or_formation[c];
					paleodb_collections$member[vote_member] <- "";
					}
				}
			}
		}
	}

keep <- match(c("formation","member","stratgroup","zone","early_interval","late_interval","max_ma","min_ma"),colnames(paleodb_collections));
if (!is.null(paleodb_collections$ma_lb))	{
	keep <- c(keep,match(c("ma_lb","ma_ub"),colnames(paleodb_collections)));
	}
if (!is.null(paleodb_collections$interval_lb))	{
	keep <- c(keep,match(c("interval_lb","interval_ub"),colnames(paleodb_collections)));
	}
if (!is.null(paleodb_collections$direct_ma))	{
	keep <- c(keep,match(c("direct_ma","direct_ma_error"),colnames(paleodb_collections)));
	}
remove <- (1:ncol(paleodb_collections))[!(1:ncol(paleodb_collections)) %in% keep];
rock_info <- paleodb_collections;
rock_info <- rock_info[,-remove];
rock_info <- rock_info[,order(keep)];
#ncolls <- nrow(rock_info);
#rrr <- (1:ncolls)[rock_info$formation %in% "Shuayb"]
#rock_info[rrr,]
rock_info <- unique(rock_info);

return(rock_info)
}

# routine to organize Paleobiologyy Database rock units to facilitate matching to external rock unit database
accio_rock_unit_thesaurus_given_paleodb_collections	<- function(paleodb_rocks)	{
# standardize group usage (these often are not entered)
formation_names <- unique(paleodb_rocks$formation[paleodb_rocks$formation!=""]);
for (fn in 1:length(formation_names))	{
	named_rock_unit <- paleodb_rocks$group[paleodb_rocks$formation==formation_names[fn]];
	named_rock_unit <- unique(named_rock_unit[named_rock_unit!=""]);
	if (length(named_rock_unit)==1)	{
		paleodb_rocks$group[paleodb_rocks$formation==formation_names[fn]] <- named_rock_unit;
		} else if (length(named_rock_unit)>1)	{
		formation_group <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=T,delete_rock_type=T);
		examples <- c();
		groupies <- paleodb_rocks$group[paleodb_rocks$formation==formation_names[fn]];
		for (fg in 1:length(formation_group))
			examples <- c(examples,sum(groupies==formation_group[fg]));
		paleodb_rocks$group[paleodb_rocks$formation==formation_names[fn]] <- unique(named_rock_unit[named_rock_unit!=""])[match(max(examples),examples)]
		}
	}
paleodb_rocks <- unique(paleodb_rocks);

n_rocks <- nrow(paleodb_rocks);
rock_no <- 1:n_rocks;
formation_no <- match(paleodb_rocks$formation,formation_names);
formation_no[is.na(formation_no)] <- 0;

#formation_no[(1:n_rocks)[is.na(formation_no)]] <- rock_no[(1:n_rocks)[is.na(formation_no)]];

#### shift formation & rock unit entries for group-only records to "•"
has_formation <- (1:n_rocks)[paleodb_rocks$formation!=""];
has_member <- (1:n_rocks)[paleodb_rocks$member!=""];
has_group <- (1:n_rocks)[paleodb_rocks$group!=""];
group_only <- has_group[!has_group %in% c(has_formation,has_member)];
formation_only <- has_formation[!has_formation %in% has_member];
member_only <- has_member[!has_member %in% has_formation];
formation_and_member <- has_formation[has_formation %in% has_member];

paleodb_rocks$group <- as.character(paleodb_rocks$group);
paleodb_rocks$formation <- as.character(paleodb_rocks$formation);
paleodb_rocks$member <- as.character(paleodb_rocks$member);
full_name <- paleodb_rocks$formation;
full_name[formation_and_member] <- paste(full_name[formation_and_member]," (",paleodb_rocks$member[formation_and_member],")",sep="");
full_name[member_only] <- paleodb_rocks$member[member_only];
full_name[group_only] <- paleodb_rocks$group[group_only];
formation_no[group_only] <- rock_no[group_only];

#paleodb_rocks$full_name <- as.character(paleodb_rocks$full_name);
#paleodb_rocks$formation[group_only] <- paleodb_rocks$member[group_only] <- paleodb_rocks$group[(1:n_rocks)[paleodb_rocks$group==""]] <- "";
#wagner_rocks$formation[groups_only_nos] <- "•"
#wagner_rocks$full_name[groups_only_nos] <- "§";
#wagner_rocks$member[groups_only_nos] <- "¢";
#wagner_rocks$group[(1:n_rocks)[wagner_rocks$group==""]] <- "£";

formation_clean_basic <- named_rock_unit <- as.character(paleodb_rocks$formation);
formation_clean_basic[has_formation] <- sapply(named_rock_unit[has_formation],scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
formation_clean_no_rock_formal <- formation_clean_no_rock <- formation_clean_basic;
formation_clean_no_rock[has_formation] <- sapply(named_rock_unit[has_formation],scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
formation_clean_no_rock_formal[has_formation] <- sapply(named_rock_unit[has_formation],scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);
member_clean_basic <- as.character(paleodb_rocks$member);
named_rock_unit <- as.character(member_clean_basic[has_member]);
member_clean_basic[has_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
member_clean_no_rock <- member_clean_no_rock_formal <- member_clean_basic;
member_clean_no_rock[has_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
member_clean_no_rock_formal[has_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);

group_clean_basic <- paleodb_rocks$group;
named_rock_unit <- paleodb_rocks$group[has_group];
group_clean_basic[has_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal = FALSE);
group_clean_no_rock <- group_clean_no_rock_formal <- group_clean_basic;
group_clean_no_rock[has_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = FALSE);
group_clean_no_rock_formal[has_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal = TRUE);

### suppose we starte with formation = "Burgess Shale" & member = "Lower)
####  rock_unit_clean_basic: "Burgess Shale (Lower)"
####  rock_unit_clean_no_rock: "Burgess (Lower)" 
####  rock_unit_clean_no_rock_formal: "Burgess" 
rock_unit_clean_basic <- formation_clean_basic;
has_member <- (1:n_rocks)[member_clean_basic!=""];	# some members are erased
member_only <- has_member[!has_member %in% has_formation];
formation_and_member <- has_formation[has_formation %in% has_member];
rock_unit_clean_basic[formation_and_member] <- paste(formation_clean_basic[formation_and_member]," (",member_clean_basic[formation_and_member],")",sep="");
rock_unit_clean_basic[group_only] <- group_clean_basic[group_only];

rock_unit_clean_no_rock <- formation_clean_no_rock;
has_member <- (1:n_rocks)[member_clean_no_rock!=""];	# some members are erased
member_only <- has_member[!has_member %in% has_formation];
formation_and_member <- has_formation[has_formation %in% has_member];
rock_unit_clean_no_rock[formation_and_member] <- paste(formation_clean_no_rock[formation_and_member]," (",member_clean_no_rock[formation_and_member],")",sep="");
rock_unit_clean_no_rock[group_only] <- group_clean_no_rock[group_only];

rock_unit_clean_no_rock_formal <- formation_clean_no_rock_formal;
has_formation <- (1:n_rocks)[formation_clean_no_rock_formal!=""]; # some formations lost
has_member <- (1:n_rocks)[member_clean_no_rock_formal!=""]; 	# some members lost
member_only <- has_member[!has_member %in% has_formation];
formation_and_member <- has_formation[has_formation %in% has_member];
rock_unit_clean_no_rock_formal[formation_and_member] <- paste(formation_clean_no_rock_formal[formation_and_member]," (",member_clean_no_rock_formal[formation_and_member],")",sep="");

rock_no_sr <- rock_no;
m_a_f <- has_member[member_clean_no_rock_formal[has_member] %in% formation_clean_no_rock_formal];
members_as_formations <- member_clean_no_rock_formal[m_a_f];
rock_no_sr[m_a_f] <- match(members_as_formations,formation_clean_no_rock_formal);
g_a_f <- group_only[group_clean_no_rock[group_only] %in% formation_clean_no_rock];
groups_as_formations <- group_clean_no_rock[g_a_f];
rock_no_sr[g_a_f] <- match(groups_as_formations,formation_clean_no_rock);
aaa <- formation_no[match(rock_no_sr,rock_no)]
formation_no_sr <- formation_no[match(rock_no_sr,rock_no)]

#senior_rock_units <- unique(rock_unit_clean_no_rock[rock_unit_clean_no_rock!=""]);
#r_u_c_n_r <- length(rock_unit_clean_no_rock[rock_unit_clean_no_rock!=""]);
#if (length(senior_rock_units) < r_u_c_n_r)	{
#	for (sr in 1:length(length(senior_rock_units)))	{
#		this_senior <- (1:r_u_c_n_r)[senior_rock_units[sr]==rock_unit_clean_no_rock];
#		if (length(this_senior)>1)	print(sr)
#		}
#	}
rock_unit_senior <- full_name[match(rock_no_sr,rock_no)];
organized_rocks <- data.frame(rock_no_sr=as.numeric(rock_no_sr),rock_no=as.numeric(rock_no),formation_no=as.numeric(formation_no),formation_no_sr=as.numeric(formation_no_sr),
			formation=as.character(paleodb_rocks$formation),member=as.character(paleodb_rocks$member),
			full_name=as.character(full_name),group=as.character(paleodb_rocks$group),
			interval_lb=as.character(paleodb_rocks$interval_lb),interval_ub=as.character(paleodb_rocks$interval_ub),
			ma_lb=as.numeric(paleodb_rocks$ma_lb),ma_ub=as.numeric(paleodb_rocks$ma_ub),
			rock_unit_senior=as.character(rock_unit_senior),formation_clean_basic=as.character(formation_clean_basic),
			formation_clean_no_rock=as.character(formation_clean_no_rock),formation_clean_no_rock_formal=as.character(formation_clean_no_rock_formal),
			member_clean_basic=as.character(member_clean_basic),member_clean_no_rock=as.character(member_clean_no_rock),
			member_clean_no_rock_formal=as.character(member_clean_no_rock_formal),group_clean_basic=as.character(group_clean_basic),
			group_clean_no_rock=as.character(group_clean_no_rock),group_clean_no_rock_formal=as.character(group_clean_no_rock_formal),
			rock_unit_clean_basic=as.character(rock_unit_clean_basic),rock_unit_clean_no_rock=as.character(rock_unit_clean_no_rock),
			rock_unit_clean_no_rock_formal=as.character(rock_unit_clean_no_rock_formal),
			stringsAsFactors=hell_no);
unique_rocks <- sort(unique(rock_no_sr));
for (ur in 1:length(unique_rocks))	{
#	this_rock <- subset(organized_rocks,organized_rocks$rock_no_sr==41);
	syns <- (1:nrow(organized_rocks))[organized_rocks$rock_no_sr==ur];
	if (length(syns)>1)	{
#		print(ur);
		oldest <- match(max(organized_rocks$ma_lb[syns]),organized_rocks$ma_lb[syns]);
		youngest <- match(min(organized_rocks$ma_ub[syns]),organized_rocks$ma_ub[syns]);
		organized_rocks$interval_lb[syns] <- as.character(organized_rocks$interval_lb[syns[oldest]]);
		organized_rocks$interval_ub[syns] <- as.character(organized_rocks$interval_ub[syns[oldest]]);
		organized_rocks$ma_lb[syns] <- as.numeric(organized_rocks$ma_lb[syns[oldest]]);
		organized_rocks$ma_ub[syns] <- as.numeric(organized_rocks$ma_ub[syns[oldest]]);
		
		if (length(unique(organized_rocks$group[syns]))>1)	{
			stratgroup <- unique(organized_rocks$group[syns]);
			stratgroup <- stratgroup[stratgroup!=""];
			if (length(stratgroup)==1)
				organized_rocks$group[syns] <- stratgroup;
			}
		}
	}

#organized_rocks <- data.frame(cbind(rock_unit_senior,formation_clean_basic,formation_clean_no_rock,formation_clean_no_rock_formal,member_clean_basic,member_clean_no_rock,member_clean_no_rock_formal,group_clean_basic,group_clean_no_rock,group_clean_no_rock_formal,rock_unit_clean_basic,rock_unit_clean_no_rock,rock_unit_clean_no_rock_formal),stringsAsFactors=hell_no);
#paleodb_rocks <- cbind(paleodb_rocks,organized_rocks);
#wagner_rocks <- cbind(wagner_rocks,as.character(rock_unit_senior),formation_clean_basic,formation_clean_no_rock,formation_clean_no_rock_formal,member_clean_basic,member_clean_no_rock,member_clean_no_rock_formal,group_clean_basic,group_clean_no_rock,group_clean_no_rock_formal,rock_unit_clean_basic,rock_unit_clean_no_rock,rock_unit_clean_no_rock_formal);
return(organized_rocks);
}

# construct a database of rock units from the PaleoDB
#	returns two tables: one for rock units, and another for zones to which rock units are assigned.
construct_stratigraphic_data_base_from_paleodb_collections <- function(paleodb_collections,zone_database="",time_scale)	{
# paleodb_collections: data.frame of collections data from PaleoDB
rocks <- accio_rock_unit_data_from_paleodb_collections(paleodb_collections,standardize_members = T);
nrock_combos <- nrow(rocks);
rocks <- rocks[order(rocks$formation,rocks$member,rocks$zone,rocks$stratgroup,-rocks$max_ma),];
has_formation <- (1:nrock_combos)[rocks$formation!=""];
has_group <- (1:nrock_combos)[rocks$stratgroup!=""];
has_member <- (1:nrock_combos)[rocks$member!=""];
keepers <- sort(unique(c(has_formation,has_group,has_member)));	# keep only information with rock names
rocks <- rocks[keepers,];
rocks <- rocks[order(rocks$formation,rocks$member,rocks$stratgroup,-rocks$max_ma),];
# redo numbering without unnamed rocks
nrock_combos <- nrow(rocks);
has_formation <- (1:nrock_combos)[rocks$formation!=""];
has_group <- (1:nrock_combos)[rocks$stratgroup!=""];
has_member <- (1:nrock_combos)[rocks$member!=""];
has_zone <- (1:nrock_combos)[rocks$zone!=""];
has_group_only <- has_group[!has_group %in% has_formation];
has_group_only <- has_group_only[!has_group_only %in% has_member];

formation_entered <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_entered[has_formation] <- sapply(named_rock_unit,scourgify_rock_unit_names);
member_entered <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_entered[has_member] <- sapply(named_rock_unit,scourgify_rock_unit_names);
group_entered <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_entered[has_group] <- sapply(named_rock_unit,scourgify_rock_unit_names);

zone <- rocks$zone;
zones_cleaned <-  sapply(zone,turgio_zone);

formation_dehyph <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_dehyph[has_formation] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE);
member_dehyph <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_dehyph[has_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE);
group_dehyph <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_dehyph[has_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE);

formation_rockless <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_rockless[has_formation] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);
member_rockless <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_rockless[has_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);
group_rockless <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_rockless[has_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE);

formation_rockless_formal <- rocks$formation;
named_rock_unit <- rocks$formation[has_formation];
formation_rockless_formal[has_formation] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);
member_rockless_formal <- rocks$member;
named_rock_unit <- rocks$member[has_member];
member_rockless_formal[has_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);
group_rockless_formal <- rocks$stratgroup;
named_rock_unit <- rocks$stratgroup[has_group];
group_rockless_formal[has_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate = TRUE,delete_rock_type = TRUE,delete_informal=TRUE);

# rock_base has all of the combinations of rock units, intervals and zones.
if (!is.null(rocks$interval_lb))	{
	rocks$early_interval <- rocks$interval_lb;
	rocks$late_interval <- rocks$interval_ub;
	}
if (!is.null(rocks$ma_lb))	{
	rocks$max_ma <- rocks$ma_lb;
	rocks$min_ma <- rocks$ma_ub;
	}
rock_base <- data.frame(formation=as.character(formation_entered),member=as.character(member_entered),group=as.character(group_entered),
			  early_interval=as.character(rocks$early_interval),late_interval=as.character(rocks$late_interval),zone = as.character(zones_cleaned),
			  max_ma=as.numeric(rocks$max_ma),min_ma=as.numeric(rocks$min_ma),
			  formation_dehyph=as.character(formation_dehyph),member_dehyph=as.character(member_dehyph),group_dehyph=as.character(group_dehyph),
			  formation_rockless=as.character(formation_rockless),member_rockless=as.character(member_rockless),group_rockless=as.character(group_rockless),
			  formation_rockless_formal=as.character(formation_rockless_formal),member_rockless_formal=as.character(member_rockless_formal),group_rockless_formal=as.character(group_rockless_formal),
			  stringsAsFactors=hell_no);
nrock_combos <- nrow(rock_base);
formation_names <- unique(formation_rockless[formation_rockless!=""]);

### START HERE!!!!
for (ff in 1:length(formation_names))	{
	ff_records <- subset(rock_base,rock_base$formation_rockless==formation_names[ff]);
	ff_r <- length(ff_records$group);
#	rb_r <- (1:nrocks)[rock_base$formation_rockless==formation_names[ff]];
	rb_r <- as.numeric(rownames(ff_records))
	if (sum(ff_records$group!="")>0 && sum(ff_records$group!="")<ff_r)	{
		ff_g <- (1:ff_r)[ff_records$group!=""][1];
		rock_base$group[rb_r] <- as.character(ff_records$group[ff_g]);
		rock_base$group_dehyph[rb_r] <- as.character(ff_records$group_dehyph[ff_g]);
		rock_base$group_rockless[rb_r] <- as.character(ff_records$group_rockless[ff_g]);
		rock_base$group_rockless_formal[rb_r] <- as.character(ff_records$group_rockless_formal[ff_g]);
		}
	formation_names_entered <- unique(ff_records$formation);
	if (length(formation_names_entered)==2)	{
#		print(ff);
		if (length(formation_names_entered[formation_names_entered!=formation_names[ff]])==2)	{
#			print(ff);
			# compare dehyphenated names to rockless names
			rock_base$formation[rb_r] <- sort(formation_names_entered[unique(formation_dehyph[rb_r])!=formation_names[ff]],decreasing = TRUE)[1];
#			sort(ff_records$formation[ff_records$formation_dehyph==ff_records$formation_rockless[1]],decreasing = TRUE)[1];
			} else	{
			rock_base$formation[rb_r] <- formation_names_entered[formation_names_entered!=formation_names[ff]];
			}
		} else if (length(formation_names_entered)>2)	{
		print(ff);
		}

	unique_members_entered <- sort(unique(ff_records$member[ff_records$member!=""]));
	unique_members_rockless <- sort(unique(ff_records$member_rockless[ff_records$member_rockless!=""]));
	
	if (length(unique_members_entered) > length(unique_members_rockless))	{
		for (mm in 1:length(unique_members_rockless))	{
			mm_records <- subset(ff_records,ff_records$member_rockless==unique_members_rockless[mm]);
			unique_members <- unique(mm_records$member);
#			if (length(unique(mm_records$member_dehyph))==1)	{
			if (length(unique_members)==2)	{
				rb_m <- rb_r[rock_base$member_rockless[rb_r]==unique_members_rockless[mm]];
#				rock_base$member[rb_m] <- unique(mm_records$member[mm_records$member!=unique_members_rockless[mm]]);
				if (length(unique_members[unique_members!=unique_members_rockless[mm]])==1)	{
					rock_base$member[rb_m] <- unique_members[unique_members!=unique_members_rockless[mm]];
					} else	{
					unique_members_dehyph <- unique(mm_records$member_dehyph);
					if (length(unique_members_dehyph[unique_members_dehyph!=unique_members_rockless[mm]])==1)	{
						rock_base$member[rb_m] <- unique_members[unique_members!=unique_members_rockless[mm]];
						}
					}
				}
			}
		}
	}

rock_base <- unique(rock_base);	# all remaining rock unit + stratigraphic assignment combinations

nrock_combos <- nrow(rock_base);
formation_member_names <- rock_base$formation_rockless_formal;
has_formation <- (1:nrock_combos)[rock_base$formation_rockless!=""];
formation_member_names[has_formation] <- rock_base$formation_rockless[has_formation];
w_members_rockless <- (1:nrock_combos)[rock_base$member_rockless!=""];
w_member_and_formation <- w_members_rockless[w_members_rockless %in% has_formation];
formation_member_names[w_member_and_formation] <- paste(as.character(rock_base$formation_rockless[w_member_and_formation])," (",as.character(rock_base$member_rockless[w_member_and_formation]),")",sep="");
rock_base <- cbind(rock_base,formation_member_names);

rocks_to_time_scale <- rocks_to_zones <- c();
#unique_rock_units <- unique(rock_base$formation[rock_base$formation!=""]);
formation_member_names <- formation_member_names[formation_member_names!=""];
unique_rock_units <- unique(sort(c(formation_names,formation_member_names)));

for (uru in 1:length(unique_rock_units))	{
	ru <- match(unique_rock_units[uru],rock_base$formation_member_names);
	if (is.na(ru))	{
		ff <- match(unique_rock_units[uru],rock_base$formation_rockless);
		frm <- rock_base$formation[ff];
		grp <- rock_base$group[ff];
		mmb <- "";
		} else	{
		grp <- rock_base$group[ru];
		frm <- rock_base$formation[ru];
		mmb <- rock_base$member[ru];
		}
	if (is.na(match(unique_rock_units[uru],formation_names)))	{
		# case where rock unit is not simply a formation
		ur_records <- subset(rock_base,rock_base$formation_member_names==unique_rock_units[uru]);
		if(ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]!="")	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			} else	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$early_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			}
		unique_zones <- sort(unique(ur_records$zone[ur_records$zone!=""]));
		nz <- length(unique_zones);
		if (nz>0)
			rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
		} else	{
		# case where rock unit is simply a formation
		ur_records <- subset(rock_base,rock_base$formation_rockless==unique_rock_units[uru]);
		if(ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]!="")	{
			# add rock info to rocks_to_time_scale
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$late_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			} else	{
			rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,unique_rock_units[uru],max(ur_records$max_ma),min(ur_records$min_ma),ur_records$early_interval[match(max(ur_records$max_ma),ur_records$max_ma)],ur_records$early_interval[match(min(ur_records$min_ma),ur_records$min_ma)]));
			}
		unique_zones <- sort(unique(ur_records$zone[ur_records$zone!=""]));
		nz <- length(unique_zones);
		if (nz>0)
			# fill out table linking rocks to zones.
			rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
		}
	}

unique_groups <- unique(rock_base$group_rockless[has_group_only]);
gu <- 0;
while (gu < length(unique_groups))	{
	gu <- gu+1;
	grp <- unique_groups[gu];
	frm <- mmb <- "";
	gr_records <- subset(rock_base,rock_base$group_rockless==grp);
	if(gr_records$late_interval[match(min(gr_records$min_ma),gr_records$min_ma)]!="")	{
		rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,"",max(gr_records$max_ma),min(gr_records$min_ma),gr_records$early_interval[match(max(gr_records$max_ma),gr_records$max_ma)],gr_records$late_interval[match(min(gr_records$min_ma),gr_records$min_ma)]));
		} else	{
		rocks_to_time_scale <- rbind(rocks_to_time_scale,c(frm,mmb,grp,"",max(gr_records$max_ma),min(gr_records$min_ma),gr_records$early_interval[match(max(gr_records$max_ma),gr_records$max_ma)],gr_records$early_interval[match(min(gr_records$min_ma),gr_records$min_ma)]));
		}
	unique_zones <- sort(unique(gr_records$zone[gr_records$zone!=""]));
	nz <- length(unique_zones);
	if (nz>0)
		rocks_to_zones <- rbind(rocks_to_zones,cbind(rep(frm,nz),rep(mmb,nz),rep(grp,nz),unique_zones));
	}

rock_database <- data.frame(formation=as.character(rocks_to_time_scale[,1]),
							member=as.character(rocks_to_time_scale[,2]),
							group=as.character(rocks_to_time_scale[,3]),
							formation_member=as.character(rocks_to_time_scale[,4]),
							interval_lb=as.character(rocks_to_time_scale[,7]),
							interval_ub=as.character(rocks_to_time_scale[,8]),
							ma_lb=as.numeric(rocks_to_time_scale[,5]),
							ma_ub=as.numeric(rocks_to_time_scale[,6]),
							stringsAsFactors=hell_no
							);
rocks_to_zones <- data.frame(formation=as.character(rocks_to_zones[,1]),member=as.character(rocks_to_zones[,2]),group=as.character(rocks_to_zones[,3]),zone=as.character(rocks_to_zones[,4]),stringsAsFactors=hell_no);
nrocks <- nrow(rock_database);
r_t_z <- nrow(rocks_to_zones);

#if (!is.null(zone_database))	{
if (is.data.frame(zone_database))	{
	formations_to_zones <- (1:nrow(rocks_to_zones))[rocks_to_zones$formation!=""];
	members_to_zones <- (1:nrow(rocks_to_zones))[rocks_to_zones$member!=""];
	formations_and_members_to_zones <- members_to_zones[members_to_zones %in% formations_to_zones];
	members_only_to_zones <- members_to_zones[!members_to_zones %in% formations_to_zones];
	rock_units_w_zones <- rocks_to_zones$formation[formations_to_zones];
	rock_units_w_zones[formations_and_members_to_zones] <- paste(rocks_to_zones$formation[formations_and_members_to_zones]," (",rocks_to_zones$member[formations_and_members_to_zones],")",sep="");
	zoned_formations <- unique(rocks_to_zones$formation[rocks_to_zones$formation!=""]);
	for (zf in 1:length(zoned_formations))	{
		this_form <- subset(rocks_to_zones,rocks_to_zones$formation==zoned_formations[zf]);
		# first collapse all of the zones into one string separated by '; '; then splitusing '; '; because multizone collections use '; '
		#		to separate zones, this will put all zones in a string 'a; b; c; d' which we then parts to c(a,b,c,d)
		form_zones <- unique(strsplit(paste(this_form$zone,collapse="; "),split="; ")[[1]]);
#		rn <- match(rocks_to_zones$formation[zf],rock_database$formation_member);
		rn <- match(zoned_formations[zf],rock_database$formation_member);
		if (is.na(rn))	{
#			rocks_to_zones$formation[zf] <- scourgify_rock_unit_names(named_rock_unit=rocks_to_zones$formation[zf],dehyphenate=T,delete_rock_type=T);
#			rn <- match(rocks_to_zones$formation[zf],rock_database$formation_member);
			try_this <- scourgify_rock_unit_names(named_rock_unit=zoned_formations[zf],dehyphenate=T,delete_rock_type=T);
			rn <- match(try_this,rock_database$formation_member);
			}
		if (is.na(rn))	{
#			rocks_to_zones$formation[zf] <- scourgify_rock_unit_names(named_rock_unit=rocks_to_zones$formation[zf],dehyphenate=T,delete_rock_type=T,delete_informal = T);
#			rn <- match(rocks_to_zones$formation[zf],rock_database$formation_member);
			try_this <- scourgify_rock_unit_names(named_rock_unit=zoned_formations[zf],dehyphenate=T,delete_rock_type=T,delete_informal = T);
			rn <- match(try_this,rock_database$formation_member);
			}
		if (!is.na(rn))	{
			zone_nos <- c();
			for (zn in 1:length(form_zones))
				zone_nos <- c(zone_nos,unique(which(zone_database==form_zones[zn],arr.ind = T)[,1]));
			zone_nos <- sort(unique(zone_nos));
			if (length(zone_nos)>0)	{
				# if zones make rock younger than oldest possible age, then adjust accordingly
				this_rock_base <- subset(rock_base,rock_base$formation==zoned_formations[zf]);
				if (nrow(this_rock_base)==0)
					this_rock_base <- subset(rock_base,rock_base$formation_rockless==zoned_formations[zf]);
				if (nrow(this_rock_base)==0)
					this_rock_base <- subset(rock_base,rock_base$formation_rockless_formal==zoned_formations[zf]);
				zoneless_formation <- subset(rock_base,rock_base$formation==zoned_formations[zf]);
				zoneless_formation <- subset(zoneless_formation,zoneless_formation$zone=="");
				if (rock_database$ma_lb[rn] > max(zone_database$ma_lb[zone_nos]) && (nrow(zoneless_formation)==0 || rock_database$ma_lb[rn] > max(zoneless_formation$max_ma)))	{
					if (nrow(zoneless_formation)==0 || max(zone_database$ma_lb[zone_nos]) > max(zoneless_formation$max_ma))	{
						oldest_zone <- zone_nos[match(max(zone_database$ma_lb[zone_nos]),zone_database$ma_lb[zone_nos])];
						rock_database$ma_lb[rn] <- zone_database$ma_lb[oldest_zone];
						rock_database$interval_lb[rn] <- zone_database$interval_lb[oldest_zone];
						} else	{
						oldest_zone <- match(max(zoneless_formation$max_ma),zoneless_formation$max_ma);
						rock_database$ma_lb[rn] <- zoneless_formation$max_ma[oldest_zone];
						rock_database$interval_lb[rn] <- zoneless_formation$early_interval[oldest_zone];
						}
					}
				# if zones make rock older than youngest possible age, then adjust accordingly
				if (rock_database$ma_ub[rn] < min(zone_database$ma_ub[zone_nos]) && (nrow(zoneless_formation)==0 || rock_database$ma_ub[rn] < min(zoneless_formation$min_ma)))	{
					if (nrow(zoneless_formation)==0 || min(zone_database$ma_ub[zone_nos]) < min(zoneless_formation$min_ma))	{
						youngest_zone <- zone_nos[match(min(zone_database$ma_ub[zone_nos]),zone_database$ma_ub[zone_nos])];
						rock_database$ma_ub[rn] <- zone_database$ma_ub[youngest_zone];
						rock_database$interval_ub[rn] <- zone_database$interval_ub[youngest_zone];
						} else	{
						youngest_zone <- match(min(zoneless_formation$min_ma),zoneless_formation$min_ma);
						rock_database$ma_ub[rn] <- zoneless_formation$min_ma[youngest_zone];
						rock_database$interval_ub[rn] <- zoneless_formation$late_interval[youngest_zone];
						}
					}
				}
			formation_members <- unique(this_form$member[this_form$member!=""]);
			fm <- 0;
			while (fm < length(formation_members))	{
				fm <- fm+1;
				whole_name <- paste(zoned_formations[zf]," (",formation_members[fm],")",sep="");
				this_rock_base <- subset(rock_base,rock_base$formation==zoned_formations[zf]);
				this_member_base <- subset(this_rock_base,this_rock_base$member==formation_members[fm]);
				memb_zones <- unique(strsplit(paste(this_member_base$zone,collapse="; "),split="; ")[[1]]);
				zoneless_member <- subset(this_member_base,this_member_base$zone=="");
				rnm <- match(whole_name,rock_database$formation_member);
				if (is.na(rnm))	{
					whole_name <- paste(zoned_formations[zf]," (",scourgify_rock_unit_names(named_rock_unit=formation_members[fm],dehyphenate=T,delete_rock_type=T),")",sep="");
					rnm <- match(whole_name,rock_database$formation_member);
					}
				if (is.na(rnm))	{
					whole_name <- paste(zoned_formations[zf]," (",scourgify_rock_unit_names(named_rock_unit=formation_members[fm],dehyphenate=T,delete_rock_type=T,delete_informal=T),")",sep="");
					rnm <- match(whole_name,rock_database$formation_member);
					}
				if (!is.na(rnm))	{
					zone_nos <- c();
					memb_zones <- memb_zones[memb_zones!=""];
					zn <- 0;
					while (zn < length(memb_zones))	{
						zn <- zn+1;
						zone_nos <- c(zone_nos,unique(which(zone_database==memb_zones[zn],arr.ind = T)[,1]));
						}
					zone_nos <- sort(unique(zone_nos));
					if (length(zone_nos)>0)	{
						# if zones make rock younger than oldest possible age, then adjust accordingly
						if (rock_database$ma_lb[rnm] > max(zone_database$ma_lb[zone_nos]) && (nrow(zoneless_member)==0 || rock_database$ma_lb[rnm] > max(zoneless_member$max_ma)))	{
							if (nrow(zoneless_member)==0 || max(zone_database$ma_lb[zone_nos]) < max(zoneless_member$max_ma))	{
								oldest_zone <- zone_nos[match(max(zone_database$ma_lb[zone_nos]),zone_database$ma_lb[zone_nos])];
								rock_database$ma_lb[rnm] <- zone_database$ma_lb[oldest_zone];
								rock_database$interval_lb[rnm] <- zone_database$interval_lb[oldest_zone];
								} else	{
								oldest_zone <- match(max(zoneless_member$max_ma),zoneless_member$max_ma);
								rock_database$ma_lb[rnm] <- zoneless_member$max_ma[oldest_zone];
								rock_database$interval_lb[rnm] <- zoneless_member$early_interval[oldest_zone];
								}
							}
						# if zones make rock older than youngest possible age, then adjust accordingly
						if (rock_database$ma_ub[rnm] < min(zone_database$ma_ub[zone_nos]) && (nrow(zoneless_member)==0 || rock_database$ma_ub[rnm] < min(zoneless_member$min_ma)))	{
							if (nrow(zoneless_member)==0 || min(zone_database$ma_ub[zone_nos]) > min(zoneless_member$min_ma))	{
								youngest_zone <- zone_nos[match(min(zone_database$ma_ub[zone_nos]),zone_database$ma_ub[zone_nos])];
								rock_database$ma_ub[rnm] <- zone_database$ma_ub[youngest_zone];
								rock_database$interval_ub[rnm] <- zone_database$interval_ub[youngest_zone];
								} else	{
								youngest_zone <- match(min(zoneless_member$min_ma),zoneless_member$min_ma);
								rock_database$ma_ub[rnm] <- zoneless_member$min_ma[youngest_zone];
								rock_database$interval_ub[rnm] <- zoneless_member$late_interval[youngest_zone];
								}
							}
						}
					}
				}
			}
#		print(zf);
		}

	group_names <- unique(rock_database$group[rock_database$formation_member==""]);
#	group_names <- group_names[group_names!=""];
	for (gg in 1:length(group_names))	{
	#	ff_records <- subset(rock_base,rock_base$formation_rockless==formation_names[ff]);
		gg_records <- unique(rbind(subset(rock_base,rock_base$group==group_names[gg]),subset(rock_base,rock_base$group_dehyph==group_names[gg]),subset(rock_base,rock_base$group_rockless==group_names[gg])));
		gb_r <- as.numeric(rownames(gg_records))
		gg_records <- subset(gg_records,gg_records$zone!="");
		gn <- (1:nrow(rock_base))[rock_base$group %in% group_names]
		group_zones <- unique(strsplit(paste(rock_base$zone[gn],collapse="; "),split="; ")[[1]]);
		group_zones <- group_zones[group_zones!=""];
		if (length(group_zones)>0)	{
			zone_nos <- c();
			for (zn in 1:length(group_zones))
				zone_nos <- c(zone_nos,unique(which(zone_database==group_zones[zn],arr.ind = T)[,1]));
			zone_nos <- sort(unique(zone_nos));
			if (length(zone_nos)>0)	{
				rn <- (1:nrow(rock_database))[rock_database$group %in% group_names[gg]];
				rn <- rn[rock_database$formation_member[rn]==""];
				if (!is.na(rn))	{
					# if zones make rock younger than oldest possible age, then adjust accordingly
					if (rock_database$ma_lb[rn] > max(zone_database$ma_lb[zone_nos]))	{
						oldest_zone <- zone_nos[match(max(zone_database$ma_lb[zone_nos]),zone_database$ma_lb[zone_nos])];
						# if this is a younger stratigraphic interval than the PaleoDB originally allowed, then change it IF there are not
						#	unzoned collections that make it older
						if (time_scale$ma_ub[match(rock_database$interval_lb[rn],time_scale$interval)] > zone_database$ma_lb[oldest_zone])	{
							other_colls <- (1:nrow(paleodb_collections))[paleodb_collections$stratgroup==group_names[gg]]
							zoneless <- other_colls[paleodb_collections$zone[other_colls]==""];
							if (length(zoneless)>0 && max(time_scale$ma_ub[match(paleodb_collections$early_interval[zoneless],time_scale$interval)]) < zone_database$ma_lb[oldest_zone])	{
								rock_database$ma_lb[rn] <- zone_database$ma_lb[oldest_zone];
								rock_database$interval_lb[rn] <- zone_database$interval_lb[oldest_zone];
								}
							}
						rock_database$interval_lb[rn] <- zone_database$interval_lb[oldest_zone];
						}
					# if zones make rock older than youngest possible age, then adjust accordingly
					if (rock_database$ma_ub[rn] < min(zone_database$ma_ub[zone_nos]))	{
						youngest_zone <- zone_nos[match(min(zone_database$ma_ub[zone_nos]),zone_database$ma_ub[zone_nos])];
						# if this is an older stratigraphic interval than the PaleoDB originally allowed, then change it IF there are not
						#	unzoned collections that make it younger
						rock_database$ma_ub[rn] <- zone_database$ma_ub[youngest_zone];
						if (time_scale$ma_lb[match(rock_database$interval_ub[rn],time_scale$interval)] < zone_database$ma_ub[youngest_zone])	{
							other_colls <- (1:nrow(paleodb_collections))[paleodb_collections$stratgroup==group_names[gg]]
							zoneless <- other_colls[paleodb_collections$zone[other_colls]==""];
							if (length(zoneless)>0 && min(time_scale$ma_lb[match(paleodb_collections$late_interval[zoneless],time_scale$interval)]) < zone_database$ma_lb[youngest_zone])	{
								rock_database$ma_ub[rn] <- zone_database$ma_ub[youngest_zone];
								rock_database$interval_ub[rn] <- zone_database$interval_ub[youngest_zone];
								}
							}
						}
					}
				}
			}
		}
	}
output <- list(rock_database,rocks_to_zones);
names(output) <- c("Rock_Database","Rocks_to_Zones");
return(output);
}

# routine to take rock unit thesaurus and provide numbers for unique rock units
match_paleodb_collections_to_paleodb_rock_thesaurus <- function(paleodb_collections,paleodb_rock_thesaurus)	{
ncolls <- nrow(paleodb_collections);
colls_w_formations <- (1:ncolls)[paleodb_collections$formation!=""];
colls_w_members <- (1:ncolls)[paleodb_collections$member!=""];
colls_w_formations_and_members <- colls_w_formations[colls_w_formations %in% colls_w_members];
colls_w_members_only <- colls_w_members[!colls_w_members %in% colls_w_formations];
colls_w_groups <- (1:ncolls)[paleodb_collections$stratgroup!=""];
colls_w_groups_only <- colls_w_groups[!colls_w_groups %in% c(colls_w_formations,colls_w_members_only)];
colls_w_rocks <- sort(unique(c(colls_w_formations,colls_w_members,colls_w_groups)));

paleodb_clean_formation_basic <- paleodb_collections$formation;
named_rock_unit <- paleodb_collections$formation[colls_w_formations];
paleodb_clean_formation_basic[colls_w_formations] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=T);

paleodb_clean_member_basic <- paleodb_collections$member;
named_rock_unit <- paleodb_collections$member[colls_w_members];
paleodb_clean_member_basic[colls_w_members] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=T);

paleodb_clean_rock_unit_basic <- paleodb_clean_formation_basic;
paleodb_clean_rock_unit_basic[colls_w_formations_and_members] <- paste(paleodb_clean_formation_basic[colls_w_formations_and_members]," (",paleodb_clean_member_basic[colls_w_formations_and_members],")",sep="");
paleodb_clean_rock_unit_basic[colls_w_members_only] <- paleodb_clean_member_basic[colls_w_members_only];

paleodb_clean_group_basic <- paleodb_collections$stratgroup;
named_rock_unit <- paleodb_collections$stratgroup[colls_w_groups];
paleodb_clean_group_basic[colls_w_groups] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=T);

paleodb_clean_formation_no_rock <- paleodb_collections$formation;
named_rock_unit <- paleodb_collections$formation[colls_w_formations];
paleodb_clean_formation_no_rock[colls_w_formations] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

paleodb_clean_member_no_rock <- paleodb_collections$member;
named_rock_unit <- paleodb_collections$member[colls_w_members];
paleodb_clean_member_no_rock[colls_w_members] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

paleodb_clean_rock_unit_no_rock <- paleodb_clean_formation_no_rock;
qqq <- colls_w_formations_and_members[paleodb_clean_member_no_rock[colls_w_formations_and_members]!=""];
paleodb_clean_rock_unit_no_rock[qqq] <- paste(paleodb_clean_formation_no_rock[qqq]," (",paleodb_clean_member_no_rock[qqq],")",sep="");
paleodb_clean_rock_unit_no_rock[colls_w_members_only] <- paleodb_clean_member_no_rock[colls_w_members_only];

paleodb_clean_group_no_rock <- paleodb_collections$stratgroup;
named_rock_unit <- paleodb_collections$stratgroup[colls_w_groups];
paleodb_clean_group_no_rock[colls_w_groups] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=T,delete_rock_type=T);

cleaned_rocks <- data.frame(cbind(paleodb_clean_formation_basic,paleodb_clean_member_basic,paleodb_clean_group_basic,paleodb_clean_formation_no_rock,paleodb_clean_member_no_rock,paleodb_clean_group_no_rock),stringsAsFactors=hell_no);

unique_rock_units <- sort(unique(paleodb_clean_rock_unit_no_rock[paleodb_clean_rock_unit_no_rock!=""]));
rock_nos <- vector(length=ncolls);
nrocks <- nrow(paleodb_rock_thesaurus);
no_finds <- c();
for (uru in 1:length(unique_rock_units))	{
	rn <- (1:nrocks)[paleodb_rock_thesaurus$rock_unit_clean_no_rock==unique_rock_units[uru]];
#	rn <- unique(which(paleodb_rock_thesaurus$rock_unit_clean_no_rock==paleodb_clean_rock_unit_no_rock[uru],arr.ind = T)[,1]);
	if (length(rn)==1)	{
		rock_nos[(1:ncolls)[paleodb_clean_rock_unit_no_rock==unique_rock_units[uru]]] <- rn;
		} else	{
		no_finds <- c(no_finds,uru);
		}
	}

unided_colls <- colls_w_rocks[rock_nos[colls_w_rocks]==0];
unided_colls_w_formations <- unided_colls[paleodb_clean_rock_unit_basic[unided_colls]!=""];
#for (ucwf in 1:length(unided_colls_w_formations))	{
#for (uc in 1:length(unided_colls))	{
uc <- 0;
while (uc < length(unided_colls))	{
	uc <- uc+1;
	coll_no <- unided_colls[uc];
	if (cleaned_rocks$paleodb_clean_member_basic[coll_no]!="")	{
		rn <- match(cleaned_rocks$paleodb_clean_member_basic[coll_no],paleodb_rock_thesaurus$formation_clean_basic);
		if (is.na(rn))
			rn <- match(cleaned_rocks$paleodb_clean_member_no_rock[coll_no],paleodb_rock_thesaurus$formation_clean_no_rock);
#		rns <- unique(which(paleodb_rock_thesaurus==cleaned_rocks$paleodb_clean_member_basic[coll_no],arr.ind = T)[,1]);
		if (!is.na(rn))	{
			# assign collections with same member to this rock number
			# clear them out of the que, too.
			assigned <- unided_colls[cleaned_rocks$paleodb_clean_member_basic[unided_colls] %in% cleaned_rocks$paleodb_clean_member_basic[coll_no]];
			rock_nos[assigned] <- rn;
			unided_colls <- unided_colls[!unided_colls %in% assigned];
			uc <- uc-1;
			}
		}
	if (rock_nos[coll_no]==0 && cleaned_rocks$paleodb_clean_group_basic[coll_no]!="")	{
		rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_basic %in% cleaned_rocks$paleodb_clean_group_basic[coll_no]];
		if (!is.na(rn) && length(rn)>1)
			rn <- rn[paleodb_rock_thesaurus$formation[rn]==""];
		if (is.na(rn))
			rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_basic %in% cleaned_rocks$paleodb_clean_rock_unit_basic[coll_no]];
		if (is.na(rn))
			rn <- (1:nrocks)[paleodb_rock_thesaurus$group_clean_no_rock %in% cleaned_rocks$paleodb_clean_rock_unit_no_rock[coll_no]];
		if (!is.na(rn) && length(rn)==1)	{
			assigned <- unided_colls[cleaned_rocks$paleodb_clean_group_basic[unided_colls] %in% cleaned_rocks$paleodb_clean_group_basic[coll_no]];
			assigned <- assigned[paleodb_collections$formation[assigned]==""];
			if (length(assigned)>0)	{
				rock_nos[assigned] <- rn;
				unided_colls <- unided_colls[!unided_colls %in% assigned];
				uc <- uc-1;
				}
			}
		}
#	print(c(uc,length(unided_colls)));
	}
return(rock_nos)
}

redate_collections_with_direct_dates <- function(collections,finest_chronostrat)	{
ncolls <- nrow(collections);
collections$direct_ma <- as.numeric(collections$direct_ma);
collections$direct_ma_error <- as.numeric(collections$direct_ma_error);
collections$direct_ma <- evanesco_na_from_vector(collections$direct_ma,0);
collections$direct_ma_error <- evanesco_na_from_vector(collections$direct_ma_error,0);
beakerheads <- (1:ncolls)[collections$direct_ma>0];
if (length(beakerheads)>0)	{
	if (!is.null(collections$ma_lb))	{
		collections$ma_lb[beakerheads] <- collections$direct_ma[beakerheads]+collections$direct_ma_error[beakerheads];
		collections$ma_ub[beakerheads] <- collections$direct_ma[beakerheads]-collections$direct_ma_error[beakerheads];
		age <- collections$ma_lb[beakerheads];
		collections$interval_lb[beakerheads] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
		age <- collections$ma_ub[beakerheads];
		collections$interval_ub[beakerheads] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
#		collections$interval_ub[beakerheads] <- rebin_collection_with_time_scale(age=collections$ma_ub[beakerheads],onset_or_end="end",fine_time_scale=finest_chronostrat);
		} else	{
		collections$max_ma[beakerheads] <- collections$direct_ma[beakerheads]+collections$direct_ma_error[beakerheads];
		collections$min_ma[beakerheads] <- collections$direct_ma[beakerheads]-collections$direct_ma_error[beakerheads];
		age <- collections$max_ma[beakerheads];
		collections$early_interval[beakerheads] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=finest_chronostrat);
		age <- collections$min_ma[beakerheads];
		collections$late_interval[beakerheads] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=finest_chronostrat);
#		collections$late_interval[beakerheads] <- rebin_collection_with_time_scale(age=collections$ma_lb[beakerheads],onset_or_end="end",fine_time_scale=finest_chronostrat);
		}
	}
return(collections);
}

## clean up & organize zone data
# construct a thesaurus for zones that provides several ways to find a zone's name.
accio_zone_thesaurus <- function(zone_data)	{
### construct a zone thesaurus.
zone_thesaurus <- cbind(zone_data$zone,zone_data$zone_sr,zone_data$zone_species,zone_data$zone_species_sr,zone_data$non_taxon_zone,zone_data$non_taxon_zone_label,zone_data$non_taxon_zone_sr,zone_data$non_taxon_zone_label_sr);
zone_thesaurus <- zone_data;
zone_thesaurus$zone_type <- zone_thesaurus$Regional_Scale <- zone_thesaurus$ma_lb <- zone_thesaurus$ma_ub <- zone_thesaurus$interval_lb <- zone_thesaurus$interval_ub <- NULL;
zone_thesaurus <- unique(zone_thesaurus);
return(zone_thesaurus);
}

# more general routine for counting things per interval
count_units_per_bin_fuzzily <- function(relevant_collections,finest_chronostrat,temporal_precision=0.1)	{
#	prob_find_this_rock_bin <- array(0,dim=c(nrow(this_rock),length(finest_chronostrat$interval)));
prob_find <- array(0,dim=c(length(finest_chronostrat$interval)));
prob_find_this_set <- c();
r_c <- 0;
while (r_c < nrow(relevant_collections))	{
	r_c <- r_c+1;
	prob_find_this_case <- array(0,dim=c(length(finest_chronostrat$interval)));
	aa <- min(finest_chronostrat$bin_first[match(relevant_collections$interval_lb[r_c],finest_chronostrat$interval)]);
	zz <- max(finest_chronostrat$bin_last[match(relevant_collections$interval_ub[r_c],finest_chronostrat$interval)]);
	if (aa==zz)	{
		prob_find_this_case[aa] <- 1;
		} else	{
		if (relevant_collections$ma_lb[r_c]!=relevant_collections$ma_ub[r_c])	{
			ma_range <- abs(relevant_collections$ma_lb[r_c]-relevant_collections$ma_ub[r_c]);
			range_start <- relevant_collections$ma_lb[r_c];
			range_end <- relevant_collections$ma_ub[r_c];
			} else	{
			ma_range <- abs(mean(finest_chronostrat$ma_lb[aa],finest_chronostrat$ma_ub[aa])-mean(finest_chronostrat$ma_lb[zz],finest_chronostrat$ma_ub[zz]));
			range_start <- mean(finest_chronostrat$ma_lb[aa],finest_chronostrat$ma_ub[aa]);
			range_end <- mean(finest_chronostrat$ma_lb[zz],finest_chronostrat$ma_ub[zz]);
			}
		round_level <- ceiling(-log10(temporal_precision));
		ma_range_finds <- array(1/(ma_range/temporal_precision),dim=c(ma_range/temporal_precision));
		ma_range_subbin_finds <- temporal_precision/ma_range;
		ma_range_find_ages <- round(seq(range_start-temporal_precision,range_end,by=-temporal_precision),round_level);
		for (bn in aa:zz)	{
			prob_bin <- ma_range_subbin_finds*sum(ma_range_find_ages[ma_range_find_ages<round(finest_chronostrat$ma_lb[bn],round_level)] %in% ma_range_find_ages[ma_range_find_ages>=round(finest_chronostrat$ma_ub[bn],round_level)]);
		#		rock_finds[rn,bn] <- max(rock_finds[rn,bn],prob_bin);
			prob_find_this_case[bn] <- prob_bin;
			}
		}
	prob_find_this_set <- rbind(prob_find_this_set,prob_find_this_case);
#	prob_find_this_set;
	}
if (length(prob_find_this_set)==nrow(finest_chronostrat))	{
	prob_find <- 1-exp(log(1-prob_find_this_set));
	}	else	{
	prob_find <- 1-exp(colSums(log(1-prob_find_this_set)));
	}
return(prob_find);
}

# clean taxon entries of ?, aff., etc.
turgio_zone <- function(zone,dbug=F)	{
if (dbug)	print(zone);
zone <- scourgify_web_text_vanilla(web_text = zone);
zone <- gsub("<sub>","",zone);
zone <- gsub("<-sub>","",zone);

zone <- gsub(", ","-",zone);
zone <- gsub("—","-",zone);
zone <- gsub("-","-",zone);
zone <- gsub("–","-",zone);
zone <- gsub("-"," - ",zone);
zone <- gsub("  "," ",zone);
zone <- gsub("  "," ",zone);
zone_molecularized <- strsplit(zone," ")[[1]];
zone_molecularized <- zone_molecularized[zone_molecularized!=""];
zone <- paste(zone_molecularized, collapse = " ");
	
zone <- gsub("  "," ",zone);
zone <- gsub("  "," ",zone);
zone <- gsub(" (\\?)","",zone);
zone <- gsub("(\\?)","",zone);
zone <- gsub(" \\?"," ",zone);
zone <- gsub("\\? ","",zone);
zone <- gsub("\\?","",zone);
zone <- gsub("—","-",zone);
zone <- gsub("-","-",zone);
zone <- gsub("–","-",zone);
zone <- gsub("/" ,"-",zone);
zone <- gsub("‚Äò,Äò","\"",zone);
zone <- gsub("“","\"",zone);
zone <- gsub("”","\"",zone);
zone <- gsub("‘","\"",zone);
zone <- gsub("’","\"",zone);
zone <- gsub("zone of ","",zone);
zone <- gsub("Zone of ","",zone);
zone <- gsub(" of the "," ",zone);
zone <- gsub("prob. ","",zone);
zone <- gsub("probably ","",zone);
zone <- gsub("MP","MP ",zone);
zone <- gsub("MP  ","MP ",zone);
zone <- gsub("MN","MN ",zone);
zone <- gsub("MN  ","MN ",zone);
zone <- gsub(" \\(Graptolite\\)" ,"",zone);
zone <- gsub(" \\(Conodont\\)" ,"",zone);
zone <- gsub(" \\(Trilobite\\)" ,"",zone);
zone <- gsub("graptolites" ,"",zone);
zone <- gsub("lower - upper ","",zone);
#zone <- gsub(" -" ,"-",zone);
#zone <- gsub("- " ,"-",zone);

if (zone != "")	{
	zone_detritus <- c("biozone","zone","zones","subzone","subzones","level","levels","bed","beds","layer","fauna","interval","local","total","range","ammonite","trilobite","conodont","coral","graptolite","reference","base","max:","min:","close","between","of","the","spore","assemblage","part","s.l.","s.l.\\)","(s.l.\\)");
	zone_molecularized <- strsplit(zone," ")[[1]];
	zone_molecularized[zone_molecularized %in% c("or","and","to","through")] <- "-";
	molecules <- (1:length(zone_molecularized))[!tolower(zone_molecularized) %in% zone_detritus];
	if (length(molecules)>0)	{
		wrong_informals <- c("earliest","early","late","latest");
		informals <- c("lowermost","lower","upper","uppermost","middle");
		if (length(molecules[tolower(zone_molecularized) %in% wrong_informals])>0)	{
			fix_these <- molecules[tolower(zone_molecularized) %in% wrong_informals];
			zone_molecularized[fix_these] <- informals[match(tolower(zone_molecularized[fix_these]),wrong_informals)];
			}
		subsz <- molecules[tolower(zone_molecularized) %in% informals];
		if (length(subsz)>0)	{
			zone_molecularized[subsz] <- tolower(zone_molecularized[subsz]);
		# if zone is entered as Polygnathus costatus upper, rewrite to upper Polygnathus costatus
		# if zone is entered as Polygnathus costatus upper - Polygnathus ensensis lower, rewrite to
		#		upper Polygnathus costatus - lower Polygnathus ensensis
			if (sum(zone_molecularized=="-")==0)	{
				zone_molecularized <- c(zone_molecularized[subsz],zone_molecularized[!molecules %in% subsz]);
				} else	{
				separators <- c(0,molecules[zone_molecularized=="-"],1+max(length(zone_molecularized)));
				for (sz in 1:length(subsz))	{
					sz_sep <- max(separators[separators<subsz[sz]]);
					if (!zone_molecularized[sz_sep+1] %in% informals)	{
						sz_sep_n <- separators[match(sz_sep,separators)+1];
						rearranged <- zone_molecularized[(sz_sep+1):(sz_sep_n-1)]
						informal_term <- rearranged[rearranged %in% informals];
						zone_molecularized[(sz_sep+1):(sz_sep_n-1)] <- rearranged <- c(informal_term,rearranged[!rearranged %in% informals]);
						}
					}
				}
			}
		zone <- paste(zone_molecularized[molecules], collapse = " ");
		if (zone_molecularized[1]=="MN" || zone_molecularized[1]=="MP")	{
			breakpts <- (1:length(zone_molecularized))[zone_molecularized %in% "-"];
			if (length(breakpts)>0)	{
				for (bp in length(breakpts):1)	{
					if (zone_molecularized[breakpts[bp]+1] != zone_molecularized[1])	{
						zone_molecularized <- c(zone_molecularized[1:breakpts[bp]],zone_molecularized[1],zone_molecularized[(breakpts[bp]+1):length(zone_molecularized)]);
						}
					}
				zone <- paste(zone_molecularized,collapse=" ");
				}
			}
		} else	{
		zone <- "";
		}
	}
zone <- gsub(" - 0" ,"0",zone);
zone <- gsub(" - 1" ,"1",zone);
zone <- gsub(" - 2" ,"2",zone);
zone <- gsub(" - 3" ,"3",zone);
zone <- gsub(" - 4" ,"4",zone);
zone <- gsub(" - 5" ,"5",zone);
zone <- gsub(" - 6" ,"6",zone);
zone <- gsub(" - 7" ,"7",zone);
zone <- gsub(" - 8" ,"8",zone);
zone <- gsub(" - 9" ,"9",zone);
zone <- zone[!zone %in% 1:9];
if (length(zone)==0)
	zone <- "";
return(zone);
}

# Separate "Redlichia chinensis - Kootenia gimmelfarbi" into "Redlichia chinensis" & "Kootenia gimmelfarbi"
diffindo_zone <- function(zone)	{
zone <- gsub("–","-",zone);
zone <- gsub(" -" ,"-",zone);
zone <- gsub("- " ,"-",zone);
multizones <- stringr::str_split(zone,pattern="-")[[1]];
return(multizones);
}

# turn Rossodus manitouensis zone to manitouensis zone
transmogrify_full_zone_names_to_species_names_only <- function(zone)	{
zone <- turgio_zone(zone);
multizones <- diffindo_zone(zone);
#for (z in 1:length(multizones))	{
z <- 0;
poss_species <- array("",dim=length(multizones));
qualifier <- c("lowermost","lower","middle","upper","uppermost");
while (z < length(multizones))	{
	z <- z+1
	this_zone <- strsplit(multizones[z],split=" ")[[1]];
	zed <- names <- length(this_zone);
	v <- "";
	while (tolower(this_zone[zed])==this_zone[zed] && zed >= 1)	{
		this_zone[zed] <- gsub(")","",this_zone[zed]);
		if (!is.na(match(this_zone[1],qualifier)))	{
			qlf <- match(this_zone[1],qualifier);
			if (v=="") {
				poss_species[z] <- paste(qualifier[qlf]," ",this_zone[zed],poss_species[z],sep=v);	# work backwards to capture subspecies names
				} else	{
				poss_species[z] <- paste(qualifier[qlf],this_zone[zed],poss_species[z],sep=v);	# work backwards to capture subspecies names
				}
			} else	{
			poss_species[z] <- paste(this_zone[zed],poss_species[z],sep=v);	# work backwards to capture subspecies names
			}
		zed <- zed-1;
		v <- " ";
		}
	}
zone_species <- "";
if (sum(poss_species!="") > 0)	{
	to_link <- (1:length(poss_species))[!poss_species %in% ""];	# do this for "Agnostus smithi - Nevadella - Redlichia goofyi" zones
	zz <- 1;
	zone_species <- poss_species[to_link[zz]];					# it will now return "smithi-goofyi"
	while (zz < length(to_link))	{
		zz <- zz+1;
		zone_species <- paste(zone_species,"-",poss_species[to_link[zz]],sep="");
		}
	}
return(zone_species)
}

# turn Rossodus manitouensis zone to Rossodus zone
transmogrify_full_zone_names_to_genus_names_only <- function(zone)	{
zone <- turgio_zone(zone);
multizones <- diffindo_zone(zone);
if (length(multizones)>0)	{
	taxon_name <- multizones;
	genus_name <- sapply(taxon_name,diffindo_genus_names_from_species_names);
	return(paste(genus_name,collapse="-"));
	} else	{
	return("");
	}
}

accio_genus_subgenus_combinations_for_zones <- function(zone)	{
zone <- turgio_zone(zone);
genus_name <- transmogrify_full_zone_names_to_genus_names_only(zone);
output <- c("","");
if (length(genus_name)>0)	{
	species_epithet <- transmogrify_full_zone_names_to_species_names_only(zone);
	if (length(species_epithet)>0)	{
		gsg <- strsplit(genus_name,split=" ")[[1]];
		if (length(gsg)==2)	{
			g_s_g <- strsplit(gsg[2],split="")[[1]];
			if (g_s_g[1]=="(")	{
				gsg[2] <- paste(g_s_g[2:(length(g_s_g)-1)],collapse="");
				output <- c(paste(gsg[1],species_epithet),paste(gsg[2],species_epithet));
				} else	{
				output <- c(zone,zone);
				}
			} else if (length(gsg)==1)	{
			output <- c(zone,zone);
			}
		}
	}
return(output);
}

# for zones called things like "Zone G1" or "Ashgill Shelly Zone 1"
### returns vector that turns zone="Zone G-1 (Hintzeia celsaora)" into:
###		zone_name[1] = "Zone G1"
###		zone_name[2] = "G1"
aparecium_nontaxon_zone <- function(zone)	{
this_zone <- stringr::str_split(zone,pattern=" ")[[1]];
zone_title <- match("zone",tolower(this_zone));
zone_name <- array("",dim=2);
names(zone_name) <- c("non_taxon_zone","non_taxon_zone_label");
#zone_name <- "";
if (!is.na(zone_title))	{
	if (zone_title < length(this_zone))	{
		zone_name[1] <- this_zone[1];
		for (zt in 2:(zone_title+1))	{
			zone_name[1] <- paste(zone_name[1],this_zone[zt],sep=" ");
			}
		zone_name[1] <- gsub("-", "",zone_name[1]);
		zone_name[1] <- gsub("—", "",zone_name[1]);
		zone_name[2] <- this_zone[zone_title+1];
		zone_name[2] <- gsub("-", "",zone_name[2]);
		zone_name[2] <- gsub("—", "",zone_name[2]);
		}
	}
#zone_name <- data_frame("F")
return(zone_name);
}

### Routines to deal with vagaries of occurrence ids
# routine to find and eliminate informal species designations that the PaleoDB considers species IDs
evanesco_indeterminate_species <- function(paleodb_finds)	{
taxon_name <- sort(unique(paleodb_finds$accepted_name));
species_epithets <- sapply(taxon_name,diffindo_species_epithets);
indet_species <- c("sp.");
indet_species <- c(indet_species,paste("sp.",LETTERS));
indet_species <- c(indet_species,paste("sp.",letters));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp. nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("indet.",i));
species_epithets <- sapply(species_epithets,accio_embedded_informal_names);
taxon_names <- taxon_name[!species_epithets %in% indet_species];
echino_species <- taxon_names <- unique(taxon_names);
taxon_names <- sapply(echino_species,echinoscrub);
taxon_names <- taxon_names[taxon_names!=""];
paleodb_finds <- paleodb_finds[paleodb_finds$accepted_name %in% taxon_names,];
return(paleodb_finds);
}

# get those hard to find informals!
accio_embedded_informal_names <- function(species_epithets)	{
for (se in 1:length(species_epithets))	{
	if (length(strsplit(species_epithets[se],split=" ")[[1]])>1)	{
		if (sum(strsplit(species_epithets[se],split=" ")[[1]] %in% c(0:9,letters,LETTERS,"genus","sp.","indet."))>0)
			species_epithets[se] <- "sp. 1";
		}
	}
return(species_epithets)
}

# pentameral funkiness begone!
echinoscrub <- function(echino_species)	{
echinobabble <- c("columnals","debris","stem","stems","holdfast","holdfasts","miscellanea","miscellaneus","ossicle","ossicles","plate","plates");
j <- simplify2array(strsplit(echino_species," ")[[1]]);
if (sum(tolower(j) %in% echinobabble)>0)	{
	return("");
	} else	{
	return(echino_species);
	}
}

# remove duplicate occurrences due to synonomies and/or lumping collections
evanesco_duplicate_occurrences_paleodb <- function(occurrences)	{
ttl_finds_cnd <- ttl_finds <- nrow(occurrences);
species <- sort(unique(occurrences$accepted_name));
otu <- length(species)

reduced_finds <- data.frame(collection_no=as.numeric(occurrences$collection_no),
							accepted_name=as.character(occurrences$accepted_name));
reduced_finds <- unique(reduced_finds);
keepers <- as.numeric(rownames(reduced_finds));
return(occurrences[keepers,]);
}

# update species names from newer opinions
transmogrify_accepted_species_name <- function(identified_name,accepted_genus)	{
#transmogrify_accepted_species_name <- function(name_to_fix)	{
#print(name_to_fix[1:2])
#accepted_genus <- name_to_fix[1];
species_epithet <- diffindo_species_epithets(identified_name);
return(paste(accepted_genus,species_epithet))
}

# PaleoDB occurrences report Genus (Subgenus); this provides Subgenus as a separate field
add_subgenus_names_to_paleodb_finds <- function(paleodb_finds)	{
# routine to provide a separate column giving subgenus names to paleodb occurence records
genus_name <- paleodb_finds$genus;
genus_subgenus <- sapply(genus_name,diffindo_subgenus_names_from_genus_names);
paleodb_finds$genus <- as.character(genus_subgenus[1,]);
subgenus <- as.character(genus_subgenus[2,]);
subgenus[subgenus==""] <- paleodb_finds$genus[subgenus==""];
if (is.na(match("subgenus",colnames(paleodb_finds))))	{
	paleodb_finds <- tibble::add_column(paleodb_finds, subgenus=as.character(subgenus), .after = match("genus",colnames(paleodb_finds)));
	} else	{
	paleodb_finds$subgenus <- as.character(subgenus);
	}
return(paleodb_finds);
}

## Routines to get classic first & last appearance data
sepkoskify_paleodb_data_one_taxon <- function(taxon,pbdb_finds,transpose=FALSE)	{
taxon_record <- subset(pbdb_finds,pbdb_finds$accepted_name==taxon);
if (nrow(taxon_record)==0)	{
	# fill vectors with nothing if no finds.
	taxon_found <- which(pbdb_finds==taxon,arr.ind = T);
	if (length(taxon_found)>0)	{
		# taxon found under another name!
		taxon_columns <- match(c("phylum","class","order","family","genus","subgenus"),colnames(pbdb_finds));
		if (sum(taxon_columns %in% taxon_found[1,2])>0)	{
			coll_w_taxon <- which(pbdb_finds==taxon,arr.ind = T)[,1];
			taxon_record <- pbdb_finds[coll_w_taxon,];
			} else if (sum(unique(colnames(paleodb_finds)[taxon_found[,2]]) %in% "identified_name")>0)	{
			coll_w_taxon <- which(pbdb_finds==taxon,arr.ind = T)[,1];
			taxon_record <- pbdb_finds[coll_w_taxon,];
			}
		}
	}

if (nrow(taxon_record)==0)	{
	# fill vectors with nothing if no finds.
	ma_max <- ma_min <- bin_lb <- bin_ub <- 0;
	} else {
	taxon_record_set <- subset(taxon_record,taxon_record$bin_lb==taxon_record$bin_ub);
	taxon_record_fuzzy <- subset(taxon_record,taxon_record$bin_lb!=taxon_record$bin_ub);
	if (nrow(taxon_record_set)>0)	{
		ma_max <- max(taxon_record_set$ma_lb);
		ma_min <- min(taxon_record_set$ma_ub);
		bin_lb <- min(taxon_record_set$bin_lb);
		bin_ub <- max(taxon_record_set$bin_ub);
		if (nrow(taxon_record_fuzzy)>0)	{
			if (min(taxon_record_fuzzy$bin_ub) < bin_lb)	{
				bin_lb <- min(taxon_record_fuzzy$bin_ub);
				if (!is.na(match(bin_lb,taxon_record$bin_lb)))
					taxon_record$ma_lb[match(bin_lb,taxon_record$bin_lb)];
				} else if (max(taxon_record_fuzzy$bin_lb) > bin_ub)	{
				bin_ub <- max(taxon_record_fuzzy$bin_lb);
				if (!is.na(match(bin_ub,taxon_record$bin_ub)))
					taxon_record$ma_ub[match(bin_ub,taxon_record$bin_ub)];
				}
			}
		} else if (nrow(taxon_record_fuzzy)>0)	{
		ma_max <- max(taxon_record$ma_lb[taxon_record$ma_lb>=max(taxon_record$ma_ub)]);
		ma_min <- min(taxon_record$ma_ub[taxon_record$ma_ub<=min(taxon_record$ma_lb)]);
		bin_lb <- min(taxon_record$bin_lb[taxon_record$bin_lb<=min(taxon_record$bin_ub)]);
		bin_ub <- max(taxon_record$bin_ub[taxon_record$bin_ub>=max(taxon_record$bin_lb)]);
		}
	}
if (transpose)	{
	output <- base::t(data.frame(taxon=as.character(taxon),bin_lb=as.numeric(bin_lb),bin_ub=as.numeric(bin_ub),ma_max=as.numeric(ma_max),ma_min=as.numeric(ma_min)));
	colnames(output) <- taxon;
	} else	{
	output <- data.frame(taxon=as.character(taxon),bin_lb=as.numeric(bin_lb),bin_ub=as.numeric(bin_ub),ma_max=as.numeric(ma_max),ma_min=as.numeric(ma_min));
	rownames(output) <- taxon;
	}
return(output);
}
	
sepkoskify_paleodb_data <- function(pbdb_finds,taxon_names,interval_names="")  {
#taxon <- taxon_names;
tranpose <- FALSE;
compendium <- c();
#compend <- sapply(taxon,sepkoskify_paleodb_data_one_taxon,pbdb_finds,tranpose);
for (i in 1:length(taxon_names))
	compendium <- rbind(compendium,sepkoskify_paleodb_data_one_taxon(taxon=taxon_names[i],pbdb_finds))
#	print(c(i,nrow(compendium)));
#compend <- data.frame(base::t(sapply(notu,sepkoskify_one_taxon_from_occurrence_data,taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub)));
return(compendium);
}

sepkoskify_occurrence_data <- function(taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub,sampled_taxa,interval_names="",condense=TRUE)  {
# taxon_no: vector giving taxon number from PaleoDB or elsewhere for occurrences
# stage_no: vector given the rank of the stratigraphic unit, with 1 the oldest
# sample_no: the collection number from PaleoDB or elswhere; can also be rock unit number
# sample_age_lb: lower bound on sample age in millions of years, with 400 older than 300
# sample_age_ub: upper bound on sample age in millions of years, with 400 older than 300
# interval_names: names of the stratigraphic units, with interval_names[1] being the name going to stage_no[n]=1
# condense: if TRUE, then only information for taxa actually present is returned;
#	if FALSE, then absent numbers (e.g., 2 if 1 and 3 present) are inserted with no information
#	FALSE is useful if you are building the compendium over increasingly refined stratigraphic data
fossil_record <- data.frame(taxon_no=as.numeric(taxon_no),stage_no=as.numeric(stage_no),sample_no=as.numeric(sample_no),sample_age_lb=as.numeric(sample_age_lb),sample_age_ub=as.numeric(sample_age_ub));
taxonomic_dictionary <- unique(data.frame(taxon_no=as.numeric(taxon_no),sampled_taxa=as.character(sampled_taxa)));
taxonomic_dictionary <- taxonomic_dictionary[order(taxonomic_dictionary$taxon_no),]
ntaxa <- length(taxonomic_dictionary$sampled_taxa);
notu <- taxonomic_dictionary$taxon_no;
compend <- base::t(sapply(notu,sepkoskify_one_taxon_from_occurrence_data,taxon_no,stage_no,sample_no,sample_age_lb,sample_age_ub));
colnames(compend) <- c("bin_lb","bin_ub","ma_max","ma_min");
if (length(interval_names)>0)	{
	compendium <- data.frame(bin_lb=as.numeric(frak_it(compend[,1])),
							 bin_ub=as.numeric(frak_it(compend[,2])),
							 ma_max=as.numeric(frak_it(compend[,3])),
							 ma_min=as.numeric(frak_it(compend[,4])),
							 stage_lb=as.character(as.character(interval_names[frak_it(compend[,1])])),
							 stage_ub=as.character(as.character(interval_names[frak_it(compend[,2])])),
							 stringsAsFactors=hell_no);
	} else	{
	compendium <- data.frame(bin_lb=as.numeric(frak_it(compend[,1])),
							 bin_ub=as.numeric(frak_it(compend[,2])),
							 ma_max=as.numeric(frak_it(compend[,3])),
							 ma_min=as.numeric(frak_it(compend[,4])));
	}
rownames(compendium) <- taxonomic_dictionary$sampled_taxa;
if (!condense)	{
	absent_taxa <- (1:max(notu))[!(1:max(notu)) %in% notu];
	ataxa <- length(absent_taxa);
	all_taxa <- c(notu,absent_taxa);
	if (length(interval_names)>0)	{
		dummy <- data.frame(
		bin_lb=as.numeric(rep(0,ataxa)),
		bin_ub=as.numeric(rep(0,ataxa)),
		ma_max=as.numeric(rep(0,ataxa)),
		ma_min=as.numeric(rep(0,ataxa)),
		stage_lb=as.character(rep("",ataxa)),
		stage_ub=as.character(rep("",ataxa)),
		stringsAsFactors=hell_no);
		} else	{
		dummy <- data.frame(
		bin_lb=as.numeric(rep(0,ataxa)),
		bin_ub=as.numeric(rep(0,ataxa)),
		ma_max=as.numeric(rep(0,ataxa)),
		ma_min=as.numeric(rep(0,ataxa)));
		}
	rownames(dummy) <- absent_taxa;
	compendium <- rbind(compendium,dummy);
	compendium <- compendium[order(all_taxa),];
	}
return(compendium);
}

accio_synoptic_richness <- function(taxon_ranges,interval_names="")	{
# taxon_ranges: a taxon x 2 matrix where column 1 gives 1st appearance bin & column 2
#	gives 2nd appearance bin. Bin 10 must precede Bin 11.
#ntaxa <- nrow(taxon_ranges);
bins <- max(taxon_ranges[,2]);
mnbin <- min(taxon_ranges[taxon_ranges[,1]!=0,1]);
jacks <- c();

for (b in mnbin:bins)	jacks <- c(jacks,sum((taxon_ranges[,1]<=b)*(taxon_ranges[,2]>=b)));

if (length(interval_names)>0)	names(jacks) <- interval_names[mnbin:bins];
return(jacks)
}

##### ROUTINES TO ORGANIZE PALEODB DATA WITH EXTERNAL DATABASE ######
# routine to convert intervals to some standardized chronostratigraphic scale (e.g., international unts)
reset_paleodb_intervals_to_desired_time_scale <- function(collections,finest_chronostrat,time_scale)	{
ncolls <- nrow(collections);
collections$late_interval[collections$late_interval==""] <- collections$early_interval[collections$late_interval==""];
problem_early_intervals <- (1:ncolls)[is.na(match(collections$early_interval,finest_chronostrat$interval))];
problem_late_intervals <- (1:ncolls)[is.na(match(collections$late_interval,finest_chronostrat$interval))];
pei <- 0;
while (pei < length(problem_early_intervals))	{
	pei <- pei+1;
	coll_no <- problem_early_intervals[pei];
	int_no <- match(collections$early_interval[coll_no],time_scale$interval);
	collections$early_interval[coll_no] <- finest_chronostrat$interval[max(1,sum(time_scale$ma_lb[int_no]<=finest_chronostrat$ma_lb))];
	if (is.na(match(collections$early_interval[coll_no],finest_chronostrat$interval)))
		collections$early_interval[coll_no] <- rebin_collection_with_time_scale(age=collections$max_ma[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat)
	}

pli <- 0;
while (pli < length(problem_late_intervals))	{
	pli <- pli+1;
	coll_no <- problem_late_intervals[pli];
	int_no <- match(collections$late_interval[coll_no],time_scale$interval);

	collections$late_interval[coll_no] <- finest_chronostrat$interval[max(1,sum(time_scale$ma_ub[int_no]<finest_chronostrat$ma_lb))];
	if (is.na(collections$late_interval[coll_no]))
		collections$late_interval[coll_no] <- rebin_collection_with_time_scale(age=collections$min_ma[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
	}
return(collections);
}

# refine PaleoDB dates given finer time scale information than the PaleoDB uses.
redate_paleodb_collections_with_time_scale <- function(paleodb_collections,time_scale,zone_database)	{
# paleodb_collections: dataframe where:
#	paleodb_collections$early_interval gives oldest possible interval
#	paleodb_collections$late_interval gives youngest possible interval (this is filled if blank)
#	paleodb_collections$max_ma gives early_interal onset (initially from PaleoDB, then from Gradstein)
#	paleodb_collections$min_ma gives late_interal end (initially from PaleoDB, then from Gradstein)
# time_scale: dataframe where:
#	time_scale$interval gives interval name
#	time_scale$ma_lb gives interval onset (given Gradstein et al. 2012)
#	time_scale$ma_ub gives interval end (given Gradstein et al. 2012)
# zone_database: dataframe where:
#	we search for the entered interval names in one of several versions of the zone name
#	zone_database$interval_lb gives onset interval
#	zone_database$interval_ub gives end interval
#	zone_database$ma_lb gives interval onset (given Gradstein et al. 2012)
#	zone_database$ma_ub gives interval end (given Gradstein et al. 2012)

ncolls <- nrow(paleodb_collections);
no_late <- (1:ncolls)[paleodb_collections$late_interval==""];
paleodb_collections$late_interval[no_late] <- paleodb_collections$early_interval[no_late];
early_intervals <- match(paleodb_collections$early_interval,time_scale$interval);
late_intervals <- match(paleodb_collections$late_interval,time_scale$interval);
#paste(paleodb_collections$collection_no[paleodb_collections$early_interval %in% "Freboldi"],collapse=",");
if (sum(is.na(early_intervals))>0 || sum(is.na(late_intervals))>0)	{
	poss_zones <- sort(unique(c(paleodb_collections$early_interval[is.na(early_intervals)],paleodb_collections$late_intervals[is.na(late_intervals)])));
	trouble <- trouble_pz <- c();
	for (pz in 1:length(poss_zones))	{
#		pz <- pz+1;
		xxx <- unique(c(which(zone_database==poss_zones[pz],arr.ind = T)[,1],which(zone_database==tolower(poss_zones[pz]),arr.ind = T)[,1]));
		this_trouble <- c();
		if (length(xxx)==1)	{
			paleodb_collections$early_interval[paleodb_collections$early_interval %in% poss_zones[pz]] <- zone_database$interval_lb[xxx];
			paleodb_collections$late_interval[paleodb_collections$late_interval %in% poss_zones[pz]] <- zone_database$interval_ub[xxx];
			} else if (length(xxx)>1)	{
			if(length(unique(zone_database$interval_lb[xxx]))==1)	{
				temp_zones <- paleodb_collections$zone[paleodb_collections$early_interval %in% poss_zones[pz]]
				temp_zones[temp_zones==""] <- zone_database$zone_sr[xxx[1]];
				paleodb_collections$zone[paleodb_collections$early_interval %in% poss_zones[pz]] <- temp_zones;
				paleodb_collections$early_interval[paleodb_collections$early_interval %in% poss_zones[pz]] <- zone_database$interval_lb[xxx[1]];
				} else	{
#				print(paste("eff me down",pz));
				#zone_database$zone[xxx];
				this_trouble <- c(this_trouble,xxx);
				}
			if(length(unique(zone_database$interval_ub[xxx]))==1)	{
				paleodb_collections$late_interval[paleodb_collections$late_interval %in% poss_zones[pz]] <- zone_database$interval_ub[xxx[1]];
#				paleodb_collections$min_ma[paleodb_collections$late_interval %in% poss_zones[pz]];
				} else	{
#				print(paste("eff me up",pz));
				#print(zone_database[xxx,]);
				this_trouble <- c(this_trouble,xxx);
				}
			}
		trouble<- c(trouble,(unique(this_trouble)));
		trouble_pz <- c(trouble_pz,rep(pz,length(unique(this_trouble))));
		}
#	zone_database[trouble,];
	}
early_intervals <- match(paleodb_collections$early_interval,time_scale$interval);
late_intervals <- match(paleodb_collections$late_interval,time_scale$interval);
#old_maxes <- collections$max_ma
#old_mins <- collections$min_ma
# sum(is.na(late_intervals))
#(1:nrow(paleodb_collections))[is.na(late_intervals)]
#paleodb_collections$early_interval[(1:nrow(paleodb_collections))[is.na(late_intervals)]]
#paleodb_collections$late_interval[(1:nrow(paleodb_collections))[is.na(late_intervals)]]
paleodb_collections$max_ma[!is.na(early_intervals)] <- time_scale$ma_lb[early_intervals[!is.na(early_intervals)]];
paleodb_collections$min_ma[!is.na(late_intervals)] <- time_scale$ma_ub[late_intervals[!is.na(late_intervals)]];
misentered <- (1:ncolls)[paleodb_collections$max_ma==paleodb_collections$min_ma];
if (length(misentered)>0)	{
	dummy <- paleodb_collections$early_interval[misentered];
	paleodb_collections$early_interval[misentered] <- paleodb_collections$late_interval[misentered];
	paleodb_collections$late_interval[misentered] <- dummy;
	paleodb_collections$max_ma[misentered] <- time_scale$ma_lb[match(paleodb_collections$early_interval[misentered],time_scale$interval)];
	paleodb_collections$min_ma[misentered] <- time_scale$ma_ub[match(paleodb_collections$late_interval[misentered],time_scale$interval)];
	}
return(paleodb_collections);
}

# refine PaleoDB dates given the zone information
redate_paleodb_collections_with_zones <- function(paleodb_collections,zone_matches,zone_database,time_scale,emend_paleodb=TRUE)	{
# paleodb_collections: dataframe of collections data downloaded from PaleoDB
# zone_matches: vector giving zones (as numbers) for collections with ';' separating different zones
# zone_database: dataframe with
#	dataframe$ma_lb: onset of zone (with 485.4 = 485.4 million years ago)
#	dataframe$ma_ub: end of zone (with 443.4 = 443.4 million years ago)
# chronostrat: dataframe with chronostratigraphic information.
#	note: it helps to make this include onlyy the finest intervals
# emend_paleodb: if T, then max_ma, min_ma, early_interval & late_interval are edited
#	if F, then ma_lb, ma_ub, interval_lb & interval_ub are appended to paleodb_collections
chronostrat_units <- unique(c(unique(zone_database$interval_lb),unique(zone_database$interval_ub),unique(paleodb_collections$early_interval),unique(paleodb_collections$late_interval)));
chronostrat_units <- chronostrat_units[chronostrat_units!=""];
time_scale <- subset(time_scale,time_scale$ma_ub<1.25*max(paleodb_collections$max_ma));
if (length(unique(time_scale$scale[match(chronostrat_units,time_scale$interval)]))==1)	{
	chronostrat <- accio_hierarchical_timescale(chronostrat_units,time_scale,regional_scale=time_scale$scale[match(chronostrat_units[1],time_scale$interval)]);
	finest_chronostrat <- chronostrat[chronostrat$bin_first==chronostrat$bin_last,];
	} else	{
	chronostrat <- accio_hierarchical_timescale(chronostrat_units,time_scale);
	chronostrat_b <- accio_hierarchical_timescale(chronostrat_units,time_scale=subset(time_scale,time_scale$scale=="International"));
	finest_chronostrat <- chronostrat_b[chronostrat_b$bin_first==chronostrat_b$bin_last,];
	}
finest_chronostrat$ma_lb <- 0.001*round(finest_chronostrat$ma_lb/0.001,0);
finest_chronostrat$ma_ub <- 0.001*round(finest_chronostrat$ma_ub/0.001,0);

ncolls <- nrow(paleodb_collections);
colls_w_zones <- (1:ncolls)[zone_matches!=""];
c_w_z <- length(colls_w_zones);
ma_lb <- 0.001*round(paleodb_collections$max_ma/0.001,0);
ma_ub <- 0.001*round(paleodb_collections$min_ma/0.001,0);
interval_lb <- paleodb_collections$early_interval;
redate_these <- (1:ncolls)[is.na(match(interval_lb,chronostrat$interval))];
new_intervals <- chronostrat$interval[match(paleodb_collections$max_ma[redate_these],chronostrat$ma_lb)];
interval_lb[redate_these[!is.na(new_intervals)]] <- new_intervals[!is.na(new_intervals)];
# cases where max_ma is within interval
if (sum(is.na(new_intervals))>0)	{
	re_redate_these <- redate_these[is.na(new_intervals)];
	for (rdt in 1:length(re_redate_these))	{
		bin <- sum(paleodb_collections$max_ma[re_redate_these[rdt]]<=chronostrat$ma_lb)
		interval_lb[redate_these[match(re_redate_these[rdt],redate_these)]] <- chronostrat$interval[bin];
		}
	}
# fill in blank late intervals with early_interval
paleodb_collections$late_interval[paleodb_collections$late_interval==""] <- paleodb_collections$early_interval[paleodb_collections$late_interval==""];
interval_ub <- paleodb_collections$late_interval;
redate_these <- (1:ncolls)[is.na(match(interval_ub,chronostrat$interval))];
new_intervals <- chronostrat$interval[match(paleodb_collections$min_ma[redate_these],chronostrat$ma_ub)];
interval_ub[redate_these[!is.na(new_intervals)]] <- new_intervals[!is.na(new_intervals)];
# cases where min_ma is within interval
if (sum(is.na(new_intervals))>0)	{
	re_redate_these <- redate_these[is.na(new_intervals)];
	for (rdt in 1:length(re_redate_these))	{
		bin <- sum(paleodb_collections$min_ma[re_redate_these[rdt]]<=chronostrat$ma_ub)
		interval_ub[redate_these[match(re_redate_these[rdt],redate_these)]] <- chronostrat$interval[bin];
		}
	}

# make sure to retain original interval_ub if new one is adjacent
#badness <- c();
#ages <- cbind(ma_lb,ma_ub);
for (cz in 1:c_w_z)	{
	coll_no <- colls_w_zones[cz];
	zone_nos <- as.numeric(strsplit(zone_matches[coll_no],";")[[1]]);
	entered_zone <- paleodb_collections$zone[coll_no];
	relv_zones <- zone_database[zone_nos,];
	zone_types <- unique(relv_zones$zone_type);

	group_zones <- relv_zones[1:length(zone_types),];
#	group_zones$zone <- group_zones$zone_sr <- ""
	for (zt in 1:length(zone_types))	{
		group_zones$zone_type[zt] <- zone_types[zt];
		type_zones <- subset(relv_zones,relv_zones$zone_type==zone_types[zt]);
		onset <- match(max(type_zones$ma_lb),type_zones$ma_lb);
		end <- match(min(type_zones$ma_ub),type_zones$ma_ub);
		group_zones$ma_lb[zt] <- type_zones$ma_lb[onset];
		group_zones$interval_lb[zt] <- type_zones$interval_lb[onset];
		group_zones$ma_ub[zt] <- type_zones$ma_ub[end];
		group_zones$interval_ub[zt] <- type_zones$interval_ub[end];
		}
#	zone_overlap <- accio_minimum_joint_ranges(lb=c(ma_lb[coll_no],group_zones$ma_lb),ub=c(ma_ub[coll_no],group_zones$ma_ub));
	# get overlap between zones and original range
	for (gz in 1:nrow(group_zones))	{
		if (gz==1)	{
			zone_overlap <- accio_temporal_overlap(lb1=ma_lb[coll_no],ub1=ma_ub[coll_no],lb2=group_zones$ma_lb[gz],ub2=group_zones$ma_ub[gz]);
			} else	{
			zone_overlap <- rbind(zone_overlap,accio_temporal_overlap(lb1=ma_lb[coll_no],ub1=ma_ub[coll_no],lb2=group_zones$ma_lb[gz],ub2=group_zones$ma_ub[gz]));
			}
		}
	if (sum(zone_overlap$ma_lb>0)>0)	{
		zone_overlap <- subset(zone_overlap,zone_overlap$ma_lb>0);
		# if zones are disjunct within alloted range, then the latest ma_lb will be younger than the oldest ma_ub
		ma_lb[coll_no] <- round(min(zone_overlap$ma_lb[zone_overlap$ma_lb>max(zone_overlap$ma_ub)]),3);
		# if zones are disjunct within alloted range, then the oldest ma_ub will be older than the youngest ma_lb
		ma_ub[coll_no] <- round(max(zone_overlap$ma_ub[zone_overlap$ma_ub<min(zone_overlap$ma_lb)]),3);
		interval_lb[coll_no] <- rebin_collection_with_time_scale(age=ma_lb[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat);
		interval_ub[coll_no] <- rebin_collection_with_time_scale(age=ma_ub[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
		} else if (sum(zone_overlap$ma_lb==0)==nrow(zone_overlap) && entered_zone!="")	{
#		print(paste("totally effed",cz,coll_no));
		if (nrow(group_zones)==1)	{
			ma_lb[coll_no] <- round(group_zones$ma_lb[1],3);
			ma_ub[coll_no] <- round(group_zones$ma_ub[1],3);
#			interval_lb[coll_no] <- group_zones$interval_lb[1];
#			interval_ub[coll_no] <- group_zones$interval_ub[1];
			interval_lb[coll_no] <- rebin_collection_with_time_scale(age=ma_lb[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat);
			interval_ub[coll_no] <- rebin_collection_with_time_scale(age=ma_ub[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
			} else	{
#			joint_ranges <- accio_minimum_joint_ranges(lb=group_zones$ma_lb,ub=group_zones$ma_ub);
			ma_lb[coll_no] <- round(min(group_zones$ma_lb[group_zones$ma_lb>max(group_zones$ma_ub)]),3);
			ma_ub[coll_no] <- round(max(group_zones$ma_ub[group_zones$ma_ub<min(group_zones$ma_lb)]),3);
			interval_lb[coll_no] <- rebin_collection_with_time_scale(age=ma_lb[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat);
			interval_ub[coll_no] <- rebin_collection_with_time_scale(age=ma_ub[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
#			interval_lb[coll_no] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>=ma_lb[coll_no])];
#			interval_ub[coll_no] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>ma_ub[coll_no])];
			}
#		} else	{
		# in principle, we reach this only if we found 2+ zone taxa in the occcurrences AND if
		#	1+ of those zone taxa has a range within the time entered into the PaleoDB
#		ok_zones <- zone_overlap[zone_overlap$ma_lb!=0,];
#		if (sum(ok_zones$ma_lb>max(ok_zones$ma_ub))>0)	{
#			ma_lb[coll_no] <- round(min(ok_zones$ma_lb[ok_zones$ma_lb>max(ok_zones$ma_ub)]),3);
#			ma_ub[coll_no] <- round(max(ok_zones$ma_ub[ok_zones$ma_ub<min(ok_zones$ma_lb)]),3);
#			interval_lb[coll_no] <- rebin_collection_with_time_scale(age=ma_lb[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat);
#			interval_ub[coll_no] <- rebin_collection_with_time_scale(age=ma_ub[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
#			interval_lb[coll_no] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>=ma_lb[coll_no])];
#			interval_ub[coll_no] <- finest_chronostrat$interval[sum(finest_chronostrat$ma_lb>ma_ub[coll_no])];
#			} else	{
#			badness <- c(badness,paste("partiallyy effed",cz,coll_no));
#			}
		}
	}

if (emend_paleodb)	{
	paleodb_collections$max_ma <- ma_lb;
	paleodb_collections$min_ma <- ma_ub;
	paleodb_collections$early_interval <- interval_lb;
	paleodb_collections$late_interval <- interval_ub;
	} else	{
	paleodb_collections <- cbind(paleodb_collections,ma_lb=as.numeric(ma_lb),ma_ub=as.numeric(ma_ub),interval_lb=as.character(interval_lb),interval_ub=as.character(interval_ub));
	}
return(paleodb_collections);
}

refine_paleodb_collection_dates_with_zone_data_only <- function(paleodb_collections,paleodb_finds,zone_database,time_scale,hierarchical_chronostrat,finest_chronostrat,examine_finds=T,temporal_precision=0.05)	{
# paleodb_collections: collections downloaded from Paleobiology Database
# paleodb_finds: occurrences downloaded from Paleobiology Database	
zone <- unique(paleodb_collections$zone[paleodb_collections$zone!=""]);
zone_ids <- sapply(zone,match_one_collection_zone_to_zone_database,zone_database);
zone_ids <- zone_ids[zone_ids!=""];
zone_ma_lbs <- zone_ma_ubs <- c();
for (zi in 1:length(zone_ids))	{
	zone_nos <- as.numeric(strsplit(zone_ids[zi],";")[[1]]);
	zone_ma_lbs <- c(zone_ma_lbs,zone_database$ma_lb[zone_nos]);
	zone_ma_ubs <- c(zone_ma_ubs,zone_database$ma_ub[zone_nos]);
	}
max_ma <- max(c(paleodb_collections$max_ma),zone_ma_lbs);
min_ma <- min(c(paleodb_collections$min_ma),zone_ma_ubs);

zone_database <- subset(zone_database,zone_database$ma_lb>min_ma);
zone_database <- subset(zone_database,zone_database$ma_ub<max_ma);
zone_database$ma_lb <- temporal_precision*round(zone_database$ma_lb/temporal_precision,0);
zone_database$ma_ub <- temporal_precision*round(zone_database$ma_ub/temporal_precision,0);
ma <- zone_database$ma_lb;
new_interval_lb <- sapply(ma,reassign_intervals_to_uniform_scale,uniform_time_scale=finest_chronostrat,onset=T);
zone_database$interval_lb[new_interval_lb!=""] <- new_interval_lb[new_interval_lb!=""];
ma <- zone_database$ma_ub;
new_interval_ub <- sapply(ma,reassign_intervals_to_uniform_scale,uniform_time_scale=finest_chronostrat,onset=F);
zone_database$interval_ub[new_interval_ub!=""] <- new_interval_lb[new_interval_ub!=""];
if(examine_finds)	{
	zone_info <- match_paleodb_collections_to_possible_zones(paleodb_collections,zone_database,paleodb_finds);
	} else	{
	zone_info <- match_collections_zones_to_zone_database(paleodb_collections,zone_database);
	}

# use zone matches for still more exact dates
zone_matches <- zone_info$zone_matches;
refined_collections <- redate_paleodb_collections_with_zones(paleodb_collections,zone_matches,zone_database,time_scale,emend_paleodb=F);
return(refined_collections);
}

# written 2019-08-15 to simplify life
accio_stratigraphic_ranges_from_sampled_in_bin <- function(finds_per_bin)	{
#present_in_bin <- 1*(finds_per_bin>=0.5);
ttl_bins <- 1:ncol(finds_per_bin);
taxon_ranges <- c();
for (tx in 1:nrow(finds_per_bin))	{
	p_i_b <- 1*(finds_per_bin[tx,]>=0.5);
	if (sum(p_i_b)==0)
		p_i_b[match(max(finds_per_bin[tx,]),finds_per_bin[tx,])] <- 1;
	p_i_b <- ttl_bins*p_i_b;
	p_i_b <- p_i_b[p_i_b>0];
	taxon_ranges <- rbind(taxon_ranges,c(min(p_i_b),max(p_i_b)));
	}
rownames(taxon_ranges) <- rownames(finds_per_bin);
return(taxon_ranges);
}

# written 2019-08-15 to simplify life
accio_synoptic_richness_from_sampled_in_bin <- function(finds_per_bin)	{
taxon_ranges <- accio_stratigraphic_ranges_from_sampled_in_bin(finds_per_bin);
synoptic_richness <- vector(length=max(max(taxon_ranges),ncol(finds_per_bin)));
names(synoptic_richness) <- colnames(finds_per_bin);
for (tx in 1:nrow(taxon_ranges))	{
	synoptic_richness[taxon_ranges[tx,1]:taxon_ranges[tx,2]] <- synoptic_richness[taxon_ranges[tx,1]:taxon_ranges[tx,2]]+1;
	}
return(synoptic_richness);
}

#### ROUTINES TO GET BASIC SAMPLING & DIVERSIFICATION NUMBERS INFORMATION ####
# routine to test out best uniform, exponential, beta & lognormal distributions
accio_sampling_distributions_for_RevBayes <- function(finds_per_bin,sample_units_per_bin,end_FBD="")	{
ok <- 0;
good_bins <- c();
sampled_in_bin <- colSums(round(finds_per_bin+0.001,0)>0);
# get fractional presences for all rocks
all_bin_sampling <- round(sample_units_per_bin,0);
if (end_FBD!="")	{
	last_bin <- match(end_FBD,colnames(finds_per_bin));
	} else	{
	last_bin <- ncol(finds_per_bin);
	}
sampled_in_bin <- colSums(finds_per_bin>=0.5);
interval_richness <- accio_synoptic_richness_from_sampled_in_bin(finds_per_bin);

#colMax(finds_per_bin)
for (bn in 1:last_bin)	{
	interval_name <- names(interval_richness)[bn];
	Sb <- interval_richness[bn];			# bin richness
	bn_coll <- sample_units_per_bin[bn];
	bin_occcurences <- sort(round(finds_per_bin[,bn]+0.01,0),decreasing = T);
	if (sum(bin_occcurences>0)>=4)	{
		if (bn_coll<=max(bin_occcurences))
			bn_coll <- 1+max(bin_occcurences);
		bin_no <- match(interval_name,colnames(finds_per_bin));
		tod <- sort(bin_occcurences,decreasing = T)[1:Sb];
		if (max(tod)>1)	{
			good_bins <- c(good_bins,interval_name);
			if (ok==1)	{
				sampling_uni_bin <- rbind(sampling_uni_bin,data.frame(base::t(optimize_uniform_occupancy(finds=tod,ncoll=bn_coll))));
				sampling_exp_bin <- rbind(sampling_exp_bin,data.frame(base::t(optimize_exponential_occupancy(finds=tod,ncoll=bn_coll))));
				sampling_bta_bin <- rbind(sampling_bta_bin,data.frame(base::t(optimize_beta_occupancy(finds=tod,ncoll=bn_coll))));
				sampling_lgn_bin <- rbind(sampling_lgn_bin,data.frame(base::t(optimize_lognormal_occupancy(finds=tod,ncoll=bn_coll))));
				} else	{
				sampling_uni_bin <- data.frame(base::t(optimize_uniform_occupancy(finds=tod,ncoll=bn_coll)));
				sampling_exp_bin <- data.frame(base::t(optimize_exponential_occupancy(finds=tod,ncoll=bn_coll)));
				sampling_bta_bin <- data.frame(base::t(optimize_beta_occupancy(finds=tod,ncoll=bn_coll)));
				sampling_lgn_bin <- data.frame(base::t(optimize_lognormal_occupancy(finds=tod,ncoll=bn_coll)));
				ok <- 1;
				}
			}
		}
	}
#names(all_bin_sampling) <- finest_chronostrat$interval;
#names(bin_sampling) <- good_bins;
rownames(sampling_uni_bin) <- rownames(sampling_exp_bin) <- rownames(sampling_bta_bin) <- rownames(sampling_lgn_bin) <- good_bins;
#rownames(sampling_uni_bin) <- rownames(sampling_exp_bin) <- rownames(sampling_lgn_bin) <- good_bins;
output <- list(sampling_uni_bin,sampling_exp_bin,sampling_bta_bin,sampling_lgn_bin);
#output <- list(sampling_uni_bin,sampling_exp_bin,sampling_lgn_bin);
names(output) <- c("best_uniform_sampling","best_exponential_sampling","best_beta_sampling","best_lognormal_sampling");
#names(output) <- c("sampling_uni_bin","sampling_exp_bin","sampling_lgn_bin");
return(output);
}

# routine to get the typical per-interval sampling rate given the best fit distributions for individual intervals
accio_median_per_unit_sampling_probability_from_best_distributions <- function(sampling_distributions)	{

diff_dists <-names(sampling_distributions);
diff_dists <- gsub("best_","",diff_dists);
diff_dists <- gsub("_sampling","",diff_dists);
aiccs <- c();
for (dd in 1:length(diff_dists))
	aiccs <- cbind(aiccs,sampling_distributions[[dd]]$AICc);
rownames(aiccs) <- rownames(sampling_distributions[[1]]);
colnames(aiccs) <- diff_dists;
	
bin_spans <- psi_binA <- psi_binB <- psi_binC <- rho_binA <- rho_binB <- c();
per_unit_sampling <- per_unit_sampling_A <- c();
for (bb in 1:nrow(aiccs))	{
#	bn <- match(rownames(aiccs)[bb],chronostrat$interval);
#	bin_span <- chronostrat$ma_lb[bn]-chronostrat$ma_ub[bn];
#	bin_spans <- c(bin_spans,bin_span);
	this_dist <- match(min(aiccs[bb,]),aiccs[bb,]);
	dist_info <- sampling_distributions[[this_dist]][bb,];
	priors <- exp(-aiccs[bb,]/2)/sum(exp(-aiccs[bb,]/2));
	if (diff_dists[this_dist]=="uniform")	{
		per_unit_sampling <- c(per_unit_sampling,dist_info$scale);
		} else if (diff_dists[this_dist]=="exponential")	{
		best_distribution <- scaled_exponential_distribution(sc=dist_info$scale,decay=dist_info$decay);
		per_unit_sampling <- c(per_unit_sampling,median(best_distribution[1:dist_info$richness]));
		} else if (diff_dists[this_dist]=="beta")	{
		best_distribution <- beta_distribution(shape1 = dist_info$alpha,shape2=dist_info$beta,S=dist_info$richness)
		per_unit_sampling <- c(per_unit_sampling,median(best_distribution));
		} else if (diff_dists[this_dist]=="lognormal")	{
		best_distribution <- scaled_lognormal_distribution(sc=dist_info$scale,mag=dist_info$mag_var,S=dist_info$richness);
		per_unit_sampling <- c(per_unit_sampling,median(best_distribution));
		}
	}
names(per_unit_sampling) <- rownames(sampling_distributions[[1]]);
return(per_unit_sampling);
}

# routine to get the typical per-collection or per-formation sampling rate given the best fit distributions for individual intervals
accio_median_prob_find_per_unit_finds_given_multiple_distributions <- function(sampling_distributions)	{
diff_dists <-names(sampling_distributions);
diff_dists <- gsub("best_","",diff_dists);
diff_dists <- gsub("_sampling","",diff_dists);
aiccs <- c();
for (dd in 1:length(diff_dists))
	aiccs <- cbind(aiccs,sampling_distributions[[dd]]$AICc);
aiccs <- data.frame(aiccs,stringsAsFactors=hell_no);
per_unit_sampling <- data.frame(array(0,dim=dim(aiccs)),stringsAsFactors=hell_no);
rownames(per_unit_sampling) <- rownames(aiccs) <- rownames(sampling_distributions[[1]]);
colnames(per_unit_sampling) <- colnames(aiccs) <- diff_dists;

for (bb in 1:nrow(aiccs))	{
#	bn <- match(rownames(aiccs)[bb],chronostrat$interval);
#	bin_span <- chronostrat$ma_lb[bn]-chronostrat$ma_ub[bn];
#	bin_spans <- c(bin_spans,bin_span);
	priors <- exp(-aiccs[bb,]/2)/sum(exp(-aiccs[bb,]/2));
	if (!is.na(match("uniform",diff_dists)))	{
		md <- match("uniform",diff_dists);
		per_unit_sampling$uniform[bb] <- sampling_distributions[[md]]$scale[bb];
		}
	if (!is.na(match("exponential",diff_dists)))	{
		md <- match("exponential",diff_dists);
		best_distribution <- scaled_exponential_distribution(sc=sampling_distributions[[md]]$scale[bb],decay=sampling_distributions[[md]]$decay[bb]);
		per_unit_sampling$exponential[bb] <- median(best_distribution[1:sampling_distributions[[md]]$richness[bb]]);
		}
	if (!is.na(match("beta",diff_dists)))	{
		md <- match("beta",diff_dists);
		best_distribution <- beta_distribution(shape1=sampling_distributions[[md]]$alpha[bb],shape2=sampling_distributions[[md]]$beta[bb],S=sampling_distributions[[md]]$richness[bb])
		per_unit_sampling$beta[bb] <- median(best_distribution);
		}
	if (!is.na(match("lognormal",diff_dists)))	{
		md <- match("lognormal",diff_dists);
		best_distribution <- scaled_lognormal_distribution(sc=sampling_distributions[[md]]$scale[bb],mag=sampling_distributions[[md]]$mag_var[bb],S=sampling_distributions[[md]]$richness[bb]);
		per_unit_sampling$lognormal[bb] <- median(best_distribution);
		}
	}
output <- list(per_unit_sampling,aiccs);
names(output) <- c("per_unit_sampling","AICc")
return(output);
}

weighted_median_sampling_probability <- function(sampling_distributions,sampling_opportunities)	{
per_bin_sampling_basics <- accio_median_prob_find_per_unit_finds_given_multiple_distributions(sampling_distributions);
sampling_priors <- exp(per_bin_sampling_basics$AICc/-2)/rowSums(exp(per_bin_sampling_basics$AICc/-2));
per_unit_sampling <- per_bin_sampling_basics$per_unit_sampling;
per_unit_weighted_prob <- c();
for (bb in 1:nrow(per_unit_sampling))	{
	per_unit_weighted_prob <- c(per_unit_weighted_prob,sum(sampling_priors[bb,]*(1-(1-per_unit_sampling[bb,]))));
	}
names(per_unit_weighted_prob) <- rownames(per_unit_sampling);
return(per_unit_weighted_prob);
}

# expected finds per million years per interval
per_interval_per_ma_psis <- function(sampling_distributions,sample_units_per_bin,bin_spans)	{
per_unit_sampling_probs <- weighted_median_sampling_probability(sampling_distributions,sampling_opportunities=sample_units_per_bin[names(sample_units_per_bin) %in% rownames(sampling_distributions[[1]])]);
bin_psi_pma <- c();
for (bn in 1:length(sample_units_per_bin))	{
	if (!is.na(match(names(sample_units_per_bin)[bn],names(per_unit_sampling_probs))))	{
		bb <- match(names(sample_units_per_bin)[bn],names(per_unit_sampling_probs));
		bin_psi_pma <- c(bin_psi_pma,probability_to_Poisson_rate(1-(1-per_unit_sampling_probs[bb])^sample_units_per_bin[bn])/bin_spans[bn]);
		} else	{
		bin_psi_pma <- c(bin_psi_pma,probability_to_Poisson_rate(1-median((1-per_unit_sampling_probs)^sample_units_per_bin[bn]))/bin_spans[bn]);
		}
	}
names(bin_psi_pma) <- names(sample_units_per_bin);
return(bin_psi_pma)
}

# expected finds per interval
per_interval_psis <- function(sampling_distributions,sample_units_per_bin)	{
# per_unit_sampling_probs gives the average probability of sampling a taxon per collection or rock unit
per_unit_sampling_probs <- weighted_median_sampling_probability(sampling_distributions,sampling_opportunities=sample_units_per_bin[names(sample_units_per_bin) %in% rownames(sampling_distributions[[1]])]);
bin_psi <- c();
for (bn in 1:length(sample_units_per_bin))	{
	if (!is.na(match(names(sample_units_per_bin)[bn],names(per_unit_sampling_probs))))	{
		bb <- match(names(sample_units_per_bin)[bn],names(per_unit_sampling_probs));
		bin_psi <- c(bin_psi,probability_to_Poisson_rate(1-(1-per_unit_sampling_probs[bb])^sample_units_per_bin[bn]));
		} else	{
		bin_psi <- c(bin_psi,probability_to_Poisson_rate(1-median((1-per_unit_sampling_probs)^sample_units_per_bin[bn])));
		}
	}
names(bin_psi) <- names(sample_units_per_bin);
return(bin_psi)
}

# first pass estimate of diversification rates given sampling
accio_initial_diversification_rates <- function(sampled_in_bin,synoptic_richness,psi_bin,chronostrat)	{

gap_taxa <- synoptic_richness-colSums(sampled_in_bin);
if (length(psi_bin)>length(gap_taxa))
	psi_bin <- psi_bin[names(psi_bin) %in% names(gap_taxa)];
chronostrat <- subset(chronostrat,chronostrat$interval %in% names(psi_bin));
bin_spans <- chronostrat$ma_lb-chronostrat$ma_ub;
names(bin_spans) <- chronostrat$interval;
ttl_taxa <- nrow(sampled_in_bin);
per_bin_origination <- per_bin_extinction <- data.frame(best_rate=rep(0,length(psi_bin)),log_likelihood=rep(0,length(psi_bin)));
for (bn in 1:length(psi_bin))	{
	# do origination
	S1 <- sum(sampled_in_bin[,bn]);
	if (bn>2)	{
		pmiss <- dpois(0,psi_bin[bn-1]*bin_spans[bn-1]);
		if (!is.na(pmiss))	{
			two_timer <- sum(sampled_in_bin[,bn]*sampled_in_bin[,bn-1]);
			two_timers <- (1:ttl_taxa)*(sampled_in_bin[,bn]*sampled_in_bin[,bn-1]);
			three_timer <- sum(sampled_in_bin[,bn]*sampled_in_bin[,bn-2]);
			three_timers <- (1:ttl_taxa)*(sampled_in_bin[,bn]*sampled_in_bin[,bn-2]);
			gap_filler <- length(three_timers[!three_timers %in% two_timers]);
			per_bin_origination[bn,] <- accio_best_diversification_given_sampling(pmiss=pmiss,S1=S1,two_timer=two_timer,gap_filler=gap_filler,continuous=T);
			} else	{
			per_bin_origination[bn,] <- c(0,0);
			}
		}  else	{
		per_bin_origination[bn,] <- c(0,0);
		}
	# do extinction
	if (bn<(ncol(sampled_in_bin)-2))	{
		pmiss <- dpois(0,psi_bin[bn+1]*bin_spans[bn+1]);
		if (!is.na(pmiss))	{
			two_timer <- sum(sampled_in_bin[,bn]*sampled_in_bin[,bn+1]);
			two_timers <- (1:ttl_taxa)*(sampled_in_bin[,bn]*sampled_in_bin[,bn+1]);
			three_timer <- sum(sampled_in_bin[,bn]*sampled_in_bin[,bn+2]);
			three_timers <- (1:ttl_taxa)*(sampled_in_bin[,bn]*sampled_in_bin[,bn+2]);
			gap_filler <- length(three_timers[!three_timers %in% two_timers]);
			per_bin_extinction[bn,] <- accio_best_diversification_given_sampling(pmiss=pmiss,S1=S1,two_timer=two_timer,gap_filler=gap_filler,continuous=T);
			} else	{
			per_bin_extinction[bn,] <- c(0,0);
			}
		}	else	{
		per_bin_extinction[bn,] <- c(0,0);
		}
	}
rownames(per_bin_origination) <- rownames(per_bin_extinction) <- names(gap_taxa);
per_bin_origination$best_rate <- per_bin_origination$best_rate/bin_spans;
per_bin_extinction$best_rate <- per_bin_extinction$best_rate/bin_spans;
per_bin_origination <- per_bin_origination[per_bin_origination$best_rate>0,];
per_bin_extinction <- per_bin_extinction[per_bin_extinction$best_rate>0,];

bin_spans <- chronostrat$ma_lb-chronostrat$ma_ub;
names(bin_spans) <- chronostrat$interval;
origination <- sum(bin_spans[match(rownames(per_bin_origination),names(bin_spans))]*per_bin_origination$best_rate)/sum(bin_spans[match(rownames(per_bin_origination),names(bin_spans))]);
extinction <- sum(bin_spans[match(rownames(per_bin_extinction),names(bin_spans))]*per_bin_extinction$best_rate)/sum(bin_spans[match(rownames(per_bin_extinction),names(bin_spans))]);
output <- c(origination,extinction);
names(output) <- c("origination","extinction");
return(output);
}

# altered version of Liow & Starrfelt's TRIPPs estimate
modified_TRIPPs <- function(psi,taxon_intervals_finds,taxon_intervals_spans)	{
dp <- c();
for (i in 1:length(taxon_intervals_finds))	{
	if (is.integer(taxon_intervals_finds[i]))	{
		dp <- c(dp,dpois(x=taxon_intervals_finds[i],lambda=psi*taxon_intervals_spans[i]))
		} else	{
		dp1 <- dpois(x=floor(taxon_intervals_finds[i]),lambda=psi*taxon_intervals_spans[i])
		dp2 <- dpois(x=ceiling(taxon_intervals_finds[i]),lambda=psi*taxon_intervals_spans[i])
		wt2 <- (ceiling(taxon_intervals_finds[i])-taxon_intervals_finds[i]);
		wt1 <- 1-wt2;
		dp <- c(dp,(wt1*dp1)+(wt2*dp2));
		}
	}
dp[dp<MINNO] <- MINNO;
lndp <- log(dp);
return(sum(lndp));
}

#### ROUTINES TO GET BASIC TAXONOMIC INFORMATION ####
accio_updated_taxonomy_for_analyzed_taxa <- function(otu_names,local_directory="",study="")	{
taxon <- otu_names;
print("Getting taxonomic data....")
#taxon_compendium <- data.frame(base::t(sapply(taxon,revelio_taxonomy_for_one_taxon)),stringsAsFactors=hell_no);
#for (nn in 1:length(otu_names))	{
#	this_taxon <- revelio_taxonomy_for_one_taxon(taxon=otu_names[nn]);
#	if (nn==1)	{
#		taxon_compendium <- this_taxon;
#		} else	{
#		taxon_compendium <- rbind(taxon_compendium,this_taxon);
#		}
#	}
initial_compendium <- data.frame(base::t(sapply(taxon,revelio_taxonomy_for_one_taxon)),stringsAsFactors=hell_no);

#for (tx in 1:length(otu_names))	{
#	if (tx==1)	{
#		initial_compendium <- revelio_taxonomy_for_one_taxon(otu_names[tx]);
#		} else	{
#		initial_compendium <- rbind(initial_compendium,revelio_taxonomy_for_one_taxon(otu_names[tx]));
#		}
#	}
initial_compendium <- evanesco_na_from_matrix(initial_compendium,"");
compendium_headers <- colnames(initial_compendium);
taxon_compendium <- data.frame(array("",dim=dim(initial_compendium)),stringsAsFactors=hell_no);
colnames(taxon_compendium) <- compendium_headers;
#for (cn in 1:ncol(initial_compendium))
#	taxon_compendium[,cn] <- simplify2array(initial_compendium[,cn]);
#taxon_compendium[1,] <- initial_compendium[1,];
for (rn in 1:nrow(initial_compendium))
	for (cn in 1:ncol(initial_compendium))
		taxon_compendium[rn,cn] <- as.character(initial_compendium[rn,cn]);
taxon_compendium <- evanesco_na_from_matrix(taxon_compendium,replacement="");

missing_taxa_info <- subset(taxon_compendium,taxon_compendium$accepted_name=="?");
found_taxa_info <- subset(taxon_compendium,taxon_compendium$n_occs>0);

if (nrow(missing_taxa_info)>0)	{
	print("Not all of your taxa are entered into the PaleoDB: you should fix that!");
	if (study=="")
		study <- paste(otu_names[1],"Clade");
	file_name <- paste(local_directory,study,"_Unentered_Taxa.csv",sep="");
	write.csv(missing_taxa_info,file=file_name,row.names = FALSE);
	}
#taxon_list <- simplify2array(found_taxa_info$accepted_name);
#output <- list(found_taxa_info,missing_taxa_info);
return(taxon_compendium);
}

# given a list of SPECIES (or subspecies), get the taxonomic data that the PaleoDB has (if any)
accio_taxonomic_data_for_one_taxon <- function(taxon)	{
orig_taxon <- taxon;
taxon <- gsub(" ","%20",taxon);
httpT <- paste("http://paleobiodb.org/data1.2/taxa/list.txt?match_name=",taxon,"&show=attr,parent,class,refattr,crmod",sep="");
#             http://paleobiodb.org/data1.2/taxa/opinions.txt?base_name=Cypraeidae&rank=species,subspecies&op_type=all
fetch <- RCurl::getURL(httpT);
taxon_information <- utils::read.csv(text = fetch, header = TRUE, stringsAsFactors=hell_no);
if (!is.na(match(orig_taxon,taxon_information$taxon_name)))	{
	taxon_information <- subset(taxon_information,taxon_information$taxon_name==orig_taxon);
	} else if (!is.na(match(orig_taxon,taxon_information)))
taxon_information <- subset(taxon_information,taxon_information$flags=="B");
return(taxon_information);
}

# get the taxonomic opinions that the PaleoDB has (if any) for a taxon given its name
accio_taxonomic_opinions_for_one_taxon <- function(taxon,exact_match=T) {
# taxon: name of taxon
# exact_match: if true, then get only this taxon; if false, then get all constituent taxa, too.
taxon <- gsub(" ","%20",taxon);
taxon_opinions <- c();
if (exact_match)	{
	httpTO <- paste("http://www.paleobiodb.org/data1.2/taxa/opinions.csv?match_name=",taxon,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	} else	{
	httpTO <- paste("http://www.paleobiodb.org/data1.2/taxa/opinions.csv?base_name=",taxon,"&op_type=all&private&show=basis,ref,refattr,ent,entname,crmod",sep="");
	}
accio <- RCurl::getURL(httpTO);
taxon_opinions <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=hell_no));
return(taxon_opinions);
}

#taxon_name <- "Redlichia chinensis - Kootenia gimmelfarbi" taxon_name <- "lower Fungochitina spinfera"
# routine to separate genus name from whole name
diffindo_genus_names_from_species_names <- function(taxon_name)	{
#print(taxon_name);
j <- stringr::str_split(taxon_name,pattern=" ")[[1]];
zone_detritus <- c("basal","lowermost","lower","middle","upper","uppermost","top");
j <- j[!j %in% zone_detritus];
j <- j[j!=""];
j <- j[!j %in% c("aff.","cf.","informal")]
genus_name <- "";
if (length(j)==2)	{
	jj <- strsplit(j[2],split="",fixed=TRUE)[[1]];
	if (jj[1]=="(")	{
		genus_name <- paste(j[1:2],collapse=" ");
		} else	{
		genus_name <- j[1];
		}
	} else if (length(j)==4)	{
	genus_name <- paste(j[1:2],collapse=" ");
	} else if (length(j)==3)	{
	jj <- strsplit(j[2],split="",fixed=TRUE)[[1]];
	if (jj[1]=="(")	{
		genus_name <- paste(j[1:2],collapse=" ");
		} else	{
		genus_name <- j[1];
		}
	} else if (length(j)>1)	{
	jj <- strsplit(j[2],split="",fixed=TRUE)[[1]];
	if (jj[1]=="(")	{
		genus_name <- paste(j[1:2],collapse=" ");
		} else	{
		genus_name <- j[1];
		}
	} else if (length(j)==1)	{
		jj <- strsplit(j,split="",fixed=TRUE)[[1]];
		if (jj[1]==toupper(jj[1]))
			genus_name <- j;
	}
return(genus_name)
}

# routine to separate subgenus & genus names from genus (subgenus) name
diffindo_subgenus_names_from_genus_names <- function(genus_name)	{
#j <- stringr::str_split(genus_name,pattern=" ")[[1]];
j <- strsplit(genus_name,split=" ")[[1]];
name_1 <- j[1];
if (length(j)==1)	{
	name_2 <- "";
	} else {
	jj <- strsplit(genus_name,split="")[[1]];
	if (sum(jj=="(")==1)	{
		name_2 <- paste(strsplit(j[2],"")[[1]][2:(length(strsplit(j[2],"")[[1]])-1)],collapse="");
		} else	{
		name_2 <- "";
		}
	}
return(c(name_1,name_2));
}

# routine to separate species name from whole name
diffindo_species_epithets <- function(taxon_name)	{
#split species name (or subspecies names) off of genus species combo
#print(taxon_name);
j <- simplify2array(strsplit(taxon_name," ")[[1]]);
j <- j[j!=""];
species_name <- "";
wordy <- length(j);
if (wordy==2)	{
	if (strsplit(j[2],split="")[[1]][1]!="(")
		species_name <- j[2];
	} else if (wordy > 1) {
		if (simplify2array(strsplit(j[2],"")[[1]])[1]=="(")	{
		# subgenus
		species_name <- paste(j[3:wordy],collapse=" ");
		} else	{
		species_name <- paste(j[2:wordy],collapse=" ");
		}
	} else if (wordy==1)	{
		if(strsplit(j,split="")[[1]][1] == tolower(strsplit(j,split="")[[1]][1]))
			species_name <- j;
	}
return(species_name);
}

accio_genus_subgenus_combinations <- function(genus_species_combo)	{
print(paste(genus_species_combo,sep=" "));
genus_name <- genus_species_combo[1];
species_name <- genus_species_combo[2];
if (species_name!="")	{
	gsg <- strsplit(genus_name,split=" ")[[1]];
	if (length(gsg)==2)	{
		g_s_g <- strsplit(gsg[2],split="")[[1]];
		if (g_s_g[1]=="(")	{
			gsg[2] <- paste(g_s_g[2:(length(g_s_g)-1)],collapse="");
			} else	{
			gsg[2] <- gsg[1];
			}
		} else {
		gsg <- c(gsg,gsg);
		}
	} else	{
	gsg <- c("","");
	}
#output <- matrix("",nrow=1,ncol=2);
#output[1,1] <- paste(gsg[1],species_name);
#output[1,2] <- paste(gsg[2],species_name);
genus_species <- paste(gsg[1],species_name);
subgenus_species <- paste(gsg[2],species_name);
#return(output)
return(data.frame(genus_species=as.character(genus_species),subgenus_species=as.character(subgenus_species)));
}

revelio_taxonomy_for_one_taxon <- function(taxon,settle=F)	{
dud <- F;
orig_name_no <- 0;
orig_taxon <- taxon;
taxon <- gsub(" ","%20",taxon);
httpT <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=full,immparent,classext,attr,app",sep="");
## If the taxon is absent, then something weird will happen: kill it with fire.
accio <- RCurl::getURL(httpT);
taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=hell_no));
if (ncol(taxon_info)<=2 || (taxon_info[1,1]=="THIS REQUEST RETURNED NO RECORDS" || taxon_info=="THIS.REQUEST.RETURNED.NO.RECORDS" || nrow(which(taxon_info=="THIS REQUEST RETURNED NO RECORDS",arr.ind = T))>0))	{
#	ttl_finds <- rbind(ttl_finds,cbind(taxon_list[tx],0));
	dud <- T;
	# if unknown species, then get information about the genus
	if (settle && length(strsplit(orig_taxon,split=" ")[[1]])>1)	{
		taxon <- strsplit(orig_taxon,split=" ")[[1]][1];
		httpT <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=",taxon,"&show=full,immparent,classext,attr,app",sep="");
		accio <- RCurl::getURL(httpT);
		taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=hell_no));
		taxon_info <- taxon_info[match(taxon,taxon_info$taxon_name),];	# get rid of species & subgenera
		if (nrow(which(taxon_info=="THIS REQUEST RETURNED NO RECORDS",arr.ind = T))==0)
			dud <- F;
		}
	if (dud)	{
		httpTT <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?base_name=Lophospira&show=full,immparent,classext,attr,app",sep="");
		accioTT <- RCurl::getURL(httpTT);
		dummy_info <- data.frame(utils::read.csv(text = accioTT, header = TRUE, stringsAsFactors=hell_no));
#	output_labels <- c("orig_no","taxon_no","record_type","flags","taxon_rank","taxon_name","taxon_attr","difference","accepted_no","accepted_rank","accepted_name","parent_no","reference_no","is_extant","n_occs","firstapp_max_ma","firstapp_min_ma","lastapp_max_ma","lastapp_min_ma","early_interval","late_interval");
		output_labels <- colnames(dummy_info);
		taxon_info <- data.frame(matrix("",1,length(output_labels)),stringsAsFactors=hell_no);
		colnames(taxon_info) <- output_labels;
		taxon_info$taxon_name <- orig_taxon;
		taxon_info$taxon_attr <- taxon_info$accepted_name <- "?";
		taxon_info$n_occs <- 0;
		}
	}
if (nrow(taxon_info)>1)	{
	orig_info <- taxon_info;
	taxon_info <- subset(taxon_info, taxon_info$taxon_name==orig_taxon);
	if (nrow(taxon_info)==0)	{
		taxon_info <- orig_info;
#		opinions <- accio_taxonomic_opinions_for_one_species(species_name=orig_taxon);
#		opinions <- accio_taxonomic_opinions_for_one_species(species_name=orig_taxon);
		opinions <- accio_taxonomic_opinions_for_one_taxon(taxon=orig_taxon);
		ttl_opinions <- nrow(opinions);
		kluge_mc_kluge <- (1:ttl_opinions)[is.na(opinions$child_name)];
		if (length(kluge_mc_kluge)>0)
			opinions$child_name[kluge_mc_kluge[opinions$child_spelling_no==opinions$orig_no]] <- opinions$taxon_name[kluge_mc_kluge[opinions$child_spelling_no==opinions$orig_no]];
		orig_name <- (1:ttl_opinions)[opinions$child_name %in% orig_taxon];
		if (length(orig_name)>0)	{
			orig_name_no <- unique(opinions$child_spelling_no[orig_name]);
			taxon_info <- unique(taxon_info[!is.na(match(opinions$taxon_name[orig_name],taxon_info$taxon_name)),]);
			}
		if (nrow(taxon_info)!=1 && length(orig_name_no)==1)	{
			httpTn <- paste("http://paleobiodb.org/data1.2/taxa/list.csv?id=",orig_name_no,"&show=full,immparent,classext,attr,app",sep="");
			accio <- RCurl::getURL(httpTn);
			taxon_info <- data.frame(utils::read.csv(text = accio, header = TRUE, stringsAsFactors=hell_no));
#			taxon_info <- orig_info;
#			orig_name <- (1:ttl_opinions)[opinions$taxon_name %in% orig_taxon];
#			if (length(orig_name)==1)	{
#				taxon_info <- taxon_info[match(opinions$taxon_name[orig_name],taxon_info$taxon_name),];
#				}
			}
		### make sure that it finds some match!!!!! 
		### use child_name && child_spelling_no
		### this can be screwed up if there are two versions of the same taxon in the PaleoDB
#		if (length(orig_name)==0 || sum(is.na(match(opinions$taxon_name[orig_name],taxon_info$taxon_name)))==length(orig_name))	{
#			which(opinions==orig_taxon,arr.ind = T)[,1]
#			orig_name <- (1:ttl_opinions)[opinions$child_name %in% taxon];
#			opinions$child_spelling_no[orig_name]
#			taxon_info[unique(which(taxon_info==unique(opinions$child_spelling_no[orig_name]),arr.ind = T)[,1]),]
#			orig_name <- (1:ttl_opinions)[opinions$taxon_name %in% orig_info$taxon_name];
#			}
#		taxon_info$n_occs <- 0;
		}
	}
return(taxon_info);
}

### routines to standardize taxon names
scourgify_taxon_names <- function(taxon_name,keep_uncertainty=F)	{
# taxon_name: string giving species name
taxon_name <- gsub(" n\\. sp\\.\\?","",taxon_name);
taxon_name <- gsub(" n\\. sp\\.","",taxon_name);
taxon_name <- gsub("n\\. gen\\. ","",taxon_name);
taxon_name <- gsub(" n\\. subgen\\.","",taxon_name);
if (!keep_uncertainty)	{
	taxon_name <- gsub(" cf\\.","",taxon_name);
	taxon_name <- gsub(" aff\\.","",taxon_name);
	taxon_name <- gsub(" ex gr\\.","",taxon_name);
	taxon_name <- gsub("cf\\. ","",taxon_name);
	taxon_name <- gsub("aff\\. ","",taxon_name);
	}
taxon_name <- gsub(" informal","",taxon_name);
taxon_name <- gsub("sensu lato","",taxon_name);
taxon_name <- gsub(" sensu lato","",taxon_name);
taxon_name <- gsub("\"", "",taxon_name);
taxon_name <- gsub(" \\?" ,"",taxon_name);
taxon_name <- gsub("\\? " ,"",taxon_name);
taxon_name <- gsub("\\?" ,"",taxon_name);
taxon_name <- gsub("  " ," ",taxon_name);
taxon_name <- paste(strsplit(taxon_name,split=" ")[[1]][!strsplit(taxon_name,split=" ")[[1]] %in% paste(LETTERS,".",sep="")],collapse=" ");
#paste(LETTERS,".",sep="");
return(taxon_name)
}

revelio_uncertain_species_assignments <- function(taxon_name)	{
# taxon_name: string giving species name
flags <- "";
cleaned_name <- scourgify_taxon_names(taxon_name);
taxon_name <- gsub("n. gen. ","",taxon_name);
taxon_name <- gsub("n. sp. ","",taxon_name);
taxon_name <- gsub("n. subgen. ","",taxon_name);
species_epithet <- scourgify_taxon_names(taxon_name=diffindo_species_epithets(cleaned_name));
genus_name <- diffindo_genus_names_from_species_names(taxon_name);
subgenus_name <- diffindo_subgenus_names_from_genus_names(genus_name)[2];
if (subgenus_name!="")
	subgenus_name <- paste("(",subgenus_name,")",sep="");
modified_taxon_name <- gsub(" ex gr\\."," ex_gr\\.",taxon_name);
taxon_components <- simplify2array(strsplit(modified_taxon_name," ")[[1]]);
t_c <- 1:length(taxon_components);
if (!is.na(match("?",taxon_components)) && max(t_c[taxon_components %in% "?"])==length(taxon_components))	{
	flags <- "uncertain species";
	} else if (sum(taxon_components %in% uncertains) > 0)	{
	btc <- t_c[taxon_components %in% uncertains];
	if (!is.na(match(subgenus_name,taxon_components)) && (match(subgenus_name,taxon_components)+1) %in% btc)	{
		flags <- "uncertain species";
		} else if (subgenus_name=="" && (match(genus_name,taxon_components)+1) %in% btc)	{
		flags <- "uncertain species";
		}
	}
return(flags);
}

revelio_uncertain_genus_assignments <- function(taxon_name)	{
# taxon_name: string giving species name
flags <- "";
cleaned_name <- scourgify_taxon_names(taxon_name);
taxon_name <- gsub("n. gen. ","",taxon_name);
taxon_name <- gsub("n. sp. ","",taxon_name);
taxon_name <- gsub("n. subgen. ","",taxon_name);
#species_epithet <- scourgify_taxon_names(taxon_name=diffindo_species_epithets(cleaned_name));
whole_genus_name <- diffindo_genus_names_from_species_names(taxon_name=cleaned_name);
gen_subgen <- diffindo_subgenus_names_from_genus_names(genus_name=whole_genus_name);
genus_name <- scourgify_taxon_names(taxon_name=gen_subgen[1]);
genus_name_check <- diffindo_subgenus_names_from_genus_names(genus_name=diffindo_genus_names_from_species_names(taxon_name))[1];
subgenus_name <- scourgify_taxon_names(taxon_name=gen_subgen[2]);
if (subgenus_name!="")
	subgenus_name <- paste("(",subgenus_name,")",sep="");
modified_taxon_name <- gsub(" ex gr\\."," ex_gr\\.",taxon_name);
taxon_components <- simplify2array(strsplit(modified_taxon_name," ")[[1]]);
if (genus_name != genus_name_check || taxon_components[2]=="?" || sum(taxon_components[1] %in% uncertains)==1)	{
	genus_name <- genus_name_check
	flags <- "uncertain genus";
	} else if (sum(taxon_components %in% "?)")==1)	{
	flags <- "uncertain genus";
	}
return(flags);
}

### routine to identify uncertain taxon assignment type
identify_taxonomic_uncertainty <- function(taxon_name)	{
# taxon_name: string giving species name
flags <- revelio_uncertain_genus_assignments(taxon_name);
flags <- c(flags,revelio_uncertain_species_assignments(taxon_name));
#paste(flags,collapse=", ");
if (sum(flags!="")==2) {
	return(paste(flags,collapse=", "));
	} else if (sum(flags=="")==2)	{
	return("");
	} else	{
	return(flags[flags!=""]);
	}
}

##### TIME SCALE ROUTINES #####
# use ages to reassign collection to an interval on a particular times scale
reassign_intervals_to_uniform_scale <- function(ma,uniform_time_scale,onset,round_extremes=T)	{
# age (in millions of years)
# uniform_time_scale: time scale that is one time line, so that no interval belongs to another interval
# onset: if T, then this ma is the lower bound; if F, then this ma is the upper bound
# round_extremes: if T, then lower bound predating the first interval goes in interval 1;
#	an upper bound postdating the last interval goes in the last
if (onset)	{
	if (ma <= max(uniform_time_scale$ma_lb) && ma > min(uniform_time_scale$ma_ub))	{
		return(uniform_time_scale$interval[sum(ma<=uniform_time_scale$ma_lb)])
		} else	{
		if (round_extremes)	{
			return(uniform_time_scale$interval[1]);
			} else	{
			return("");
			}
		}
	} else	{
	if (ma >= min(uniform_time_scale$ma_ub) && ma < max(uniform_time_scale$ma_lb))	{
		return(uniform_time_scale$interval[sum(ma<=uniform_time_scale$ma_lb)]);
	} else	{
		if (round_extremes)	{
			return(uniform_time_scale$interval[nrow(uniform_time_scale)]);
			} else	{
			return("");
			}
		}
	}
}

# routine to convert intervals to some standardized chronostratigraphic scale (e.g., international unts)
accio_hierarchical_timescale <- function(chronostrat_units,time_scale,regional_scale="International",ma_fuzz=0)	{
chronostrat <- time_scale[as.character(time_scale$interval) %in% chronostrat_units,];
onset <- max(time_scale$ma_lb[as.character(time_scale$interval) %in% chronostrat_units])+ma_fuzz;
end <- min(time_scale$ma_ub[as.character(time_scale$interval) %in% chronostrat_units])-ma_fuzz;
relv_scales <- unique(c(chronostrat$scale,regional_scale));

xxx <- subset(time_scale,time_scale$ma_lb<=onset);
chronostrat <- subset(xxx,xxx$ma_ub>=end);
chronostrat <- subset(chronostrat,chronostrat$scale==regional_scale);
chronostrat <- chronostrat[chronostrat$scale %in% relv_scales,];
chronostrat <- chronostrat[order(chronostrat$ma_lb,chronostrat$ma_lb,decreasing = T),];
chronostrat <- chronostrat[order(-chronostrat$ma_lb,chronostrat$ma_lb),];

bin_onsets <- sort(unique(chronostrat$ma_lb),decreasing=T);
bin_ends <- sort(unique(chronostrat$ma_ub),decreasing=T);

n_intr <- nrow(chronostrat);	# total intervals;
if (length(bin_onsets)==length(bin_ends))	{
	parent_interval <- array("",dim=c(n_intr));
	bin_first <- match(chronostrat$ma_lb,bin_onsets);
	bin_last <- match(chronostrat$ma_ub,bin_ends);
	bin_spans <- 1+bin_last-bin_first;
	nbins <- max(bin_last);
	chronostrat <- tibble::add_column(chronostrat, bin_last=as.numeric(bin_last), .after = match("ma_ub",colnames(chronostrat)));
	chronostrat <- tibble::add_column(chronostrat, bin_first=as.numeric(bin_first), .after = match("ma_ub",colnames(chronostrat)));
	unique_bin_spans <- sort(unique(bin_spans[bin_spans>1]),decreasing=T);
	bs <- 0;
	while (bs < length(unique_bin_spans))	{
		bs <- bs+1;
		broad_bins <- (1:n_intr)[bin_spans==unique_bin_spans[bs]];
		for (bb in 1:length(broad_bins))	{
			lng_bin <- broad_bins[bb];
			if (parent_interval[lng_bin]=="")
				parent_interval[lng_bin] <- chronostrat$interval[lng_bin];
			daughter_bins <- (1:n_intr)[bin_first >= bin_first[lng_bin]][(1:n_intr)[bin_first >= bin_first[lng_bin]] %in% (1:n_intr)[bin_last <= bin_last[lng_bin]]]
			daughter_bins <- daughter_bins[!daughter_bins %in% lng_bin];
			parent_interval[daughter_bins] <- chronostrat$interval[lng_bin];
			}
		}
	if (bs==0 && (regional_scale=="Stage Slice") || (regional_scale=="Time Slice"))	{
		slice <- chronostrat$interval;
		parent_interval <- sapply(slice,accio_parent_intervals_for_stage_slice);
		} else	{
		parent_interval[(1:n_intr)[parent_interval==""]] <- chronostrat$interval[(1:n_intr)[parent_interval==""]];
		}
	chronostrat <- tibble::add_column(chronostrat, parent_interval=as.character(parent_interval), .after = match("interval",colnames(chronostrat)));
	} else	{
	bin_boundaries <- sort(unique(c(bin_onsets,bin_ends)),decreasing = T);
	bin_first <- match(chronostrat$ma_lb,bin_boundaries);
	bin_last <- match(chronostrat$ma_ub,bin_boundaries)-1;
	bin_spans <- 1 + bin_last - bin_first;
	parent_interval <- array("",dim=c(n_intr));
	chronostrat <- tibble::add_column(chronostrat, bin_last=as.numeric(bin_last), .after = match("ma_ub",colnames(chronostrat)));
	chronostrat <- tibble::add_column(chronostrat, bin_first=as.numeric(bin_first), .after = match("ma_ub",colnames(chronostrat)));
#	chronostrat <- tibble::add_column(chronostrat, parent_interval=as.character(chronostrat$interval[chronostrat_subintervals]), .after = match("interval",colnames(chronostrat)));
	unique_bin_spans <- sort(unique(bin_spans[bin_spans>1]),decreasing=T);
	bs <- 0;
	while (bs < length(unique_bin_spans))	{
		bs <- bs+1;
		broad_bins <- (1:n_intr)[bin_spans==unique_bin_spans[bs]];
		for (bb in 1:length(broad_bins))	{
			lng_bin <- broad_bins[bb];
			if (parent_interval[lng_bin]=="")
				chronostrat$parent_interval[lng_bin] <- chronostrat$interval[lng_bin];
			daughter_bins <- (1:n_intr)[bin_first >= bin_first[lng_bin]][(1:n_intr)[bin_first >= bin_first[lng_bin]] %in% (1:n_intr)[bin_last <= bin_last[lng_bin]]]
			daughter_bins <- daughter_bins[!daughter_bins %in% lng_bin];
			chronostrat$parent_interval[daughter_bins] <- chronostrat$interval[lng_bin];
			}
		}
	}
return(chronostrat);
}

# get chronostratigraphic unit to which another chronostratigraphic unit belongs (e.g., the Katian belongs to the Ordovician belongs to the Paleozoic belongs to the Phanerozoic)
accio_parent_intervals_for_stage_slice <- function(slice)	{
j <- strsplit(slice,split="",fixed=TRUE)[[1]];
jl <- length(j);
if (sum(c("a","b","c","d") %in% j[jl])==1)	jl <- jl-1;
jl <- jl-1;
slice <- paste(j[1:jl],collapse="");
}

reset_paleodb_intervals_to_desired_time_scale <- function(collections,finest_chronostrat,time_scale)	{
ncolls <- nrow(collections);
collections$late_interval[collections$late_interval==""] <- collections$early_interval[collections$late_interval==""];
problem_early_intervals <- (1:ncolls)[is.na(match(collections$early_interval,finest_chronostrat$interval))];
problem_late_intervals <- (1:ncolls)[is.na(match(collections$late_interval,finest_chronostrat$interval))];
pei <- 0;
while (pei < length(problem_early_intervals))	{
	pei <- pei+1;
	coll_no <- problem_early_intervals[pei];
	int_no <- match(collections$early_interval[coll_no],time_scale$interval);
	collections$early_interval[coll_no] <- finest_chronostrat$interval[max(1,sum(time_scale$ma_lb[int_no]<=finest_chronostrat$ma_lb))];
	if (is.na(match(collections$early_interval[coll_no],finest_chronostrat$interval)))
		collections$early_interval[coll_no] <- rebin_collection_with_time_scale(age=collections$max_ma[coll_no],onset_or_end = "onset",fine_time_scale = finest_chronostrat)
	}

pli <- 0;
while (pli < length(problem_late_intervals))	{
	pli <- pli+1;
	coll_no <- problem_late_intervals[pli];
	int_no <- match(collections$late_interval[coll_no],time_scale$interval);

	collections$late_interval[coll_no] <- finest_chronostrat$interval[max(1,sum(time_scale$ma_ub[int_no]<finest_chronostrat$ma_lb))];
	if (is.na(collections$late_interval[coll_no]))
		collections$late_interval[coll_no] <- rebin_collection_with_time_scale(age=collections$min_ma[coll_no],onset_or_end = "end",fine_time_scale = finest_chronostrat);
	}
return(collections);
}

# redo collection dates based on stratigraphic intervals
rebin_collection_with_time_scale <- function(age,onset_or_end,fine_time_scale)	{
# age: age in millions of years
# onset_or_end: "onset" for lower bound, "end" for upper bound		
# fine_time_scale: dataframe where:
#	fine_time_scale$interval gives interval name
#	fine_time_scale$ma_lb gives interval onset (given Gradstein et al. 2012)
#	fine_time_scale$ma_ub gives interval end (given Gradstein et al. 2012)
#	NOTE: fine_time_scale cannot include subintervals of other intervals in that time scale:
#	  e.g., just Cambrian, Ordovician, Silurian or Sandbian, Katian, Rhuddanian, etc.
age <- round(abs(age),3);
if (onset_or_end=="onset" || onset_or_end=="Onset")	{
	return(as.character(fine_time_scale$interval[sum(age<=round(fine_time_scale$ma_lb,3))]));
	} else	{
	return(as.character(fine_time_scale$interval[sum(age<round(fine_time_scale$ma_lb,3))]));
	}
}

# replace reported stratigraphic intervals with appropriate intervals from one time scale (usually internationally)
completely_rebin_collections_with_uniform_time_scale <- function(collections,uniform_time_scale)	{
uniform_time_scale$ma_lb <- 0.001*round(uniform_time_scale$ma_lb/0.001,0);
uniform_time_scale$ma_ub <- 0.001*round(uniform_time_scale$ma_ub/0.001,0);
if (!is.null(collections$ma_lb))	{
	ma_lb_col <- match("ma_lb",colnames(collections));
	ma_ub_col <- match("ma_ub",colnames(collections));
	} else	{
	ma_lb_col <- match("max_ma",colnames(collections));
	ma_ub_col <- match("min_ma",colnames(collections));
	}
collections[,ma_lb_col] <- 0.001*round(collections[,ma_lb_col]/0.001,0);
collections[,ma_ub_col] <- 0.001*round(collections[,ma_ub_col]/0.001,0);
if (!is.null(collections$interval_lb))	{
	int_lb_col <- match("interval_lb",colnames(collections));
	int_ub_col <- match("interval_ub",colnames(collections));
	} else	{
	int_lb_col <- match("early_interval",colnames(collections));
	int_ub_col <- match("late_interval",colnames(collections));
	}
age <- collections[,ma_lb_col];
collections[,int_lb_col] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="onset",fine_time_scale=uniform_time_scale);
age <- collections[,ma_ub_col];
collections[,int_ub_col] <- sapply(age,rebin_collection_with_time_scale,onset_or_end="end",fine_time_scale=uniform_time_scale);
if (!is.null(collections$bin_lb))	{
	collections$bin_lb <- match(collections[,int_lb_col],uniform_time_scale$interval);
	collections$bin_ub <- match(collections[,int_ub_col],uniform_time_scale$interval);
	}
return(collections);
}

# standardize "upper / late and "lower / early" -> "upper" and "early"
turgio_stage <- function(stage,dbug=FALSE)	{
if (dbug)	print(zone);
stage <- gsub(" / ","/",stage);
stage <- gsub("Early/Lower","Early",stage);
stage <- gsub("Lower","Early",stage);
stage <- gsub("Late/Upper","Late",stage);
stage <- gsub("Upper","Late",stage);
stage <- transmogrify_diacritics(funky_text=stage);
return(stage);
}

##### ROUTINES TO ORGANIZE PALEODB DATA WITH EXTERNAL DATABASE ######
# routine to download PaleoDB data and use external stratigraphic and zone databases to put
#	more exact dates on collections than PaleoDB provides
# study: name of study (anything you like)
# collections: downloaded PaleoDB data.  MUST include:
#	collection_no
#	formation, member and stratgroup
#	zone
# external_strat_database: name of an external database giving detailed information on rocks.
# MUST give:
##	rock_no: number of rock unit (specific formation + member combination)
##	rock_no_sr: number of the senior synonym (which usually is the same as rock_no)
##	formation_no: the number of the formation (to link different members of a formation)
##	formation: name of formation (might be blank if group is the only information)
##	member: name of member (usually blank if only formation or group information given)
##	full_name: formation_Name (member_Name) or formation_Name if no member.  (Blank if only group)
##	group: Stratigraphic group (often blank)
##	interval_lb: name of chronostratigraphc unit in which rock-unit starts (I prefer stage-slices)
##	interval_ub: name of chronostratigraphc unit in which rock-unit ends (I prefer stage-slices)
##  ma_lb: oldest possible age of rock unit given some chronostratigraphic scale (use 500 for 500 Ma);
##	ma_ub; oldest possible age of rock unit given some chronostratigraphic scale (use 495 for 495 Ma);
# external_zone_database: name of external database giving zones associated with rock units.
# collections <- strophomenoid_collections;
refine_collection_dates_with_external_database <- function(study="",collections,rock_database,zone_database,rock_to_zone_database,time_scale,directory="",save_files=TRUE,output_type=".csv")	{
n_rocks <- nrow(rock_database);

redone_collections <- match_paleodb_collections_to_external_stratigraphic_database(collections,wagner_rocks=rock_database);
ncoll <- nrow(redone_collections);
#lost_time <- (1:ncoll)[!redone_collections$early_interval %in% timey_wimey$interval]
#redone_collections$early_interval[lost_time]
#(1:ncolls)[redone_collections$ma_lb<=redone_collections$ma_ub]
#cbind(redone_collections$ma_lb[fff],redone_collections$ma_ub[fff])
#### But wait, there's more!  The PaleoDB has zone data for some localities.  This can further restrict dates
zone_thesaurus <- accio_zone_thesaurus(zone_database);
zd <- nrow(zone_database);
zt <- nrow(zone_thesaurus);
rzd <- nrow(rock_to_zone_database);

collections_w_zones <- (1:ncoll)[redone_collections$zone!=""];
cwz <- length(collections_w_zones);
redone_collections$zone <- as.character(redone_collections$zone)
zone <- as.character(redone_collections$zone[collections_w_zones]);
redone_collections$zone[collections_w_zones] <- sapply(zone,turgio_zone);
zone_species <- as.character(redone_collections$zone);
non_taxon_zone <- non_taxon_zone_label <- array("",dim=nrow(redone_collections));
zone_species[collections_w_zones] <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
non_taxon_zone_info <- raster::t(sapply(zone,aparecium_nontaxon_zone));
rownames(non_taxon_zone_info) <- NULL;
non_taxon_zone[collections_w_zones] <- non_taxon_zone_info[,1];
non_taxon_zone_label[collections_w_zones] <- non_taxon_zone_info[,2];
redone_collections <- cbind(redone_collections,zone_species,non_taxon_zone,non_taxon_zone_label);

collections_w_zones <- (1:ncoll)[redone_collections$zone!=""];
cwz <- length(collections_w_zones);
zones_matched <- c();

cz <- 0;
while (cz < cwz)	{
	### blows up  when a stage is  entered as the zone!!!
	cz <- cz+1;
	paldbc <- collections_w_zones[cz];	#redone_collections$collection_no[paldbc]
	zone_paleodb <- turgio_zone(zone=as.character(redone_collections$zone[paldbc]));
	if (!is.na(match(zone_paleodb,zone_database$zone)))	{
		zone_paleodb_sr <- zone_database$zone_sr[match(zone_paleodb,zone_database$zone)];
		}	else	{
		zone_paleodb_sr <- "flibberty_gibberty";
		}
	zone_spc_paleodb <- as.character(redone_collections$zone_species[paldbc]);
	zone_nontax_paleodb <- as.character(redone_collections$non_taxon_zone[paldbc]);
	zone_nontax_lab_paleodb <- as.character(redone_collections$non_taxon_zone_label[paldbc]);
	poss_matches <- (1:rzd)[rock_to_zone_database$rock_no_sr %in% redone_collections$rock_no_sr[paldbc]];
	wagner_rock_unit <- match(redone_collections$rock_no_sr[paldbc],rock_database$rock_no);
	
	### if we cannot find formation (member) combination, then just look for formation
	if (length(poss_matches)==0 || (redone_collections$rock_no_sr[paldbc]>0 && rock_database$member[wagner_rock_unit]=="")) 
		poss_matches <- (1:rzd)[rock_to_zone_database$formation_no %in% redone_collections$formation_no[paldbc]];

	# now, try to find this zone in zones associated with rock-unit
	matches <- c();
	relevant_zones <- array("",dim=c(0,5));
	relevant_zones <- data.frame(relevant_zones);
	if (length(poss_matches)>0)	{
		# match listed zones to PaleoDB entered zones
		matches1 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% zone_paleodb];
		# match senior names of zones to PaleoDB entered zones
		matches2 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% zone_paleodb];
		matches3 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% zone_paleodb_sr];
		# match PaleoDB zone to species name alone 
		matches4 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_paleodb];
		# match PaleoDB zone species name only to species name alone
		# match PaleoDB zone species name only to senior species name alone
		if (zone_spc_paleodb!="")	{
			matches5 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_spc_paleodb];
			matches6 <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% zone_spc_paleodb]
			} else	{
			matches5 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_paleodb];
			matches6 <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% zone_paleodb];
			}
		# match PaleoDB zone species name only to senior species name alone
		matches7 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% zone_nontax_paleodb];
		# match PaleoDB zone species name only to senior species name alone
		if (zone_nontax_paleodb != "")	{
			matches8 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% zone_nontax_paleodb];
			matches9 <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% zone_nontax_paleodb];
			} else	{
			matches8 <- poss_matches[rock_to_zone_database$non_taxon_zone[poss_matches] %in% zone_paleodb];
			matches9 <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% zone_paleodb];
			}
		# match PaleoDB zone name only to a zone label
		matches <- sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7,matches8,matches9)));
		if (length(matches)>0)
			relevant_zones <- rbind(relevant_zones,cbind(as.character(rock_to_zone_database$zone_sr[matches]),as.numeric(rock_to_zone_database$ma_lb[matches]),as.numeric(rock_to_zone_database$ma_ub[matches]),as.character(rock_to_zone_database$interval_lb[matches]),as.character(rock_to_zone_database$interval_ub[matches])));
		}

	## Look out for multiple zones: if that is the caes, then get ages from all of them.
	if (length(poss_matches)>0 && length(matches)==0)	{
#		print(paste("couldn't find",as.character(zone_spc_paleodb)));
		multizones <- diffindo_zone(zone=zone_paleodb);
		mz <- length(multizones);
		if (mz>1)	{
#			multimatches <- c();
			for (zz in 1:mz)	{
				# match whole name of zones to PaleoDB entered zone
				multimatches <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% multizones[zz]];
				# match senior names of zones to PaleoDB entered zones
				if (length(multimatches)==0)	{
					multimatches <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% multizones[zz]];
					# match species name only to PaleoDB entered zone
					if (length(multimatches)==0)	{
						multimatches <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% multizones[zz]];
						# match senior synonym species name only to PaleoDB entered zone
						if (length(multimatches)==0)	{
							multimatches <- poss_matches[rock_to_zone_database$zone_species_sr[poss_matches] %in% multizones[zz]];
							if (length(multimatches)==0)	{
								multimatches <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% multizones[zz]];
								if (length(multimatches)==0)	{
									multimatches <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% multizones[zz]];
									if (length(multimatches)==0)	{
										multimatches <- poss_matches[rock_to_zone_database$non_taxon_zone_label[poss_matches] %in% multizones[zz]];
										}
									}
								}
							}
						}
					}
				if (length(multimatches)>0)	{
					relevant_zones <- rbind(relevant_zones,cbind(as.character(rock_to_zone_database$zone_sr[multimatches]),as.numeric(rock_to_zone_database$ma_lb[multimatches]),as.numeric(rock_to_zone_database$ma_ub[multimatches]),as.character(rock_to_zone_database$interval_lb[multimatches]),as.character(rock_to_zone_database$interval_ub[multimatches])));
					} else if (length(multimatches)==0)	{
					mm1 <- (1:zt)[zone_thesaurus$zone %in% multizones[zz]];
					mm2 <- (1:zt)[zone_thesaurus$zone_species %in% multizones[zz]];
					mm3 <- (1:zt)[zone_thesaurus$zone_species_sr %in% multizones[zz]];
					mm4 <- (1:zt)[zone_thesaurus$non_taxon_zone %in% multizones[zz]];
					mm5 <- (1:zt)[zone_thesaurus$non_taxon_zone_label %in% multizones[zz]];
					mm6 <- (1:zt)[zone_thesaurus$non_taxon_zone_sr %in% multizones[zz]];
					mm7 <- (1:zt)[zone_thesaurus$non_taxon_zone_label_sr %in% multizones[zz]];
#					multimatches <- (1:zt)[zone_thesaurus$zone_sr[zone_thesaurus$zone_species %in% multizones[zz]] %in% rock_to_zone_data$zone[poss_matches]];
					multimatches <- sort(unique(c(mm1,mm2,mm3,mm4,mm5,mm6,mm7)));
					if (length(multimatches)>0)	{
						zone_taxa <- sort(unique(zone_thesaurus$zone_sr[multimatches]));

						poss_zones_all <- (1:zd)[zone_database$zone_sr %in% zone_taxa];
						rel_z <- cbind(as.character(zone_database$zone_sr[poss_zones_all]),as.numeric(zone_database$ma_lb[poss_zones_all]),as.numeric(zone_database$ma_ub[poss_zones_all]),as.character(zone_database$interval_lb[poss_zones_all]),as.character(zone_database$interval_ub[poss_zones_all]));
						for (ztx  in 1:length(zone_taxa))	{
							poss_zones <- subset(rel_z,rel_z[,1]==zone_taxa[ztx]);
							zlb <- max(as.numeric(poss_zones[,2]));
							zub <- min(as.numeric(poss_zones[,3]));
							if (zlb >= redone_collections$min_ma[paldbc] && zub <= redone_collections$max_ma[paldbc])
								relevant_zones <- rbind(relevant_zones,poss_zones);
							}
						}
					}
				if (length(multimatches)>0)	matches <- c(matches,multimatches);
				}
			}
	 	}
		### routine to find if a senior synonym of the reported zone is in the database
#		if (length(matches)==0 && length(strsplit(zone_paleodb," ")[[1]])==1)	{
	## Look out for zones that are not "Amorphognathus ordovicicus" but "Zone D"
	if (length(poss_matches)>0 && length(matches)==0)	{
		aa <- (1:zd)[zone_database$zone %in% zone_paleodb];
		if (length(aa)==0)	{
			bb <- (1:zd)[zone_database$zone_species %in% zone_paleodb];
			} else	{
			bb <- c();
			}
		cc <- (1:zd)[zone_database$non_taxon_zone_label %in% zone_paleodb];
		dd <- (1:zd)[zone_database$non_taxon_zone_label_sr %in% zone_paleodb];
		ee <- sort(unique(c(aa,bb,cc,dd)));
		if (length(ee)>0)	{
			try_this_zone <- unique(zone_database$zone_sr[ee]);
			try_this_zone_spc <- transmogrify_full_zone_names_to_species_names_only(try_this_zone);
			## now, repeat the matching game!
			matches1 <- poss_matches[rock_to_zone_database$zone[poss_matches] %in% try_this_zone];
			# match senior names of zones to PaleoDB entered zones
			matches2 <- poss_matches[rock_to_zone_database$zone_sr[poss_matches] %in% try_this_zone];
			matches3 <- poss_matches[rock_to_zone_database$zone_species[poss_matches] %in% try_this_zone_spc];
			matches <- c(matches,sort(unique(c(matches1,matches2,matches3))));
			if (length(matches)==0)	{
				relevant_zones <- rbind(relevant_zones,cbind(zone_database$zone_sr[ee],zone_database$ma_lb[ee],zone_database$ma_ub[ee],as.character(zone_database$interval_lb[ee]),as.character(zone_database$interval_ub[ee])));
				matches <- length(zone_database$zone_sr[ee])
				}
			}
		}

	if (length(poss_matches)==0)	{
		matches1 <- (1:zd)[zone_database$zone %in% zone_paleodb];
		if (length(matches1)==0)	{
			matches2 <- (1:zd)[zone_database$zone_species %in% zone_spc_paleodb];
			matches2 <- matches2[!zone_database$zone_species[matches2]==""];
			} else	{
			matches2 <- c();
			}
		if (zone_nontax_paleodb!="")	{
			matches3 <- (1:zd)[zone_database$non_taxon_zone %in% zone_nontax_paleodb];
			} else	{
			matches3 <- c();
			}
		if (zone_nontax_lab_paleodb!="")	{
			matches4 <- (1:zd)[zone_database$non_taxon_zone_label %in% zone_nontax_lab_paleodb];
			} else	{
			matches4 <- c();
			}
		matches <- sort(unique(c(matches1,matches2,matches3,matches4)));
		if (length(matches)==0)	{
			multizones <- diffindo_zone(zone_paleodb);
			mz <- length(multizones);
			if (mz>1)	{
				rel_z <- c();
	#			multimatches <- c();
				for (zz in 1:mz)	{
					# match zone names of zones to PaleoDB entered zones
					matches1 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					# match senior names of zones to PaleoDB entered zones
					matches2 <- (1:zd)[zone_database$zone_sr %in% multizones[zz]];
					# match species names only of zones to PaleoDB entered zones
					matches3 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
					matches4 <- (1:zd)[zone_database$zone_species_sr %in% multizones[zz]];
					matches5 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					matches6 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
					matches7 <- (1:zd)[zone_database$non_taxon_zone_label %in% multizones[zz]];
					multimatches <- c(sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7))));
					if (length(multimatches)>0)	{
						rel_z <- rbind(rel_z,cbind(as.character(zone_database$zone_sr[multimatches]),as.numeric(zone_database$ma_lb[multimatches]),as.numeric(zone_database$ma_ub[multimatches]),as.character(zone_database$interval_lb[multimatches]),as.character(zone_database$interval_ub[multimatches])));
						matches <- c(matches,multimatches);
						}
					}
			} else	{
				rel_z <- c();	
				}
			} else	{
			rel_z <- cbind(as.character(zone_database$zone_sr[matches]),as.numeric(zone_database$ma_lb[matches]),as.numeric(zone_database$ma_ub[matches]),as.character(zone_database$interval_lb[matches]),as.character(zone_database$interval_ub[matches]));
			}
		if (length(matches)>0)	{
			zone_taxa <- sort(unique(rel_z[,1]));
		# there should be only one taxon here.  If 2+ were found (e.g., "Ctenognathodus murchisoni",
		#	"Cyrtograptus murchisoni" and "Didymograptus murchisoni" for "murchisoni"), then find the appropriate one
			if (length(zone_taxa)>=1)	{
#			possible_zones <- c();
				for (ztx  in 1:length(zone_taxa))	{
					poss_zones <- subset(rel_z,rel_z[,1]==zone_taxa[ztx]);
					zlb <- max(as.numeric(poss_zones[,2]));
					zub <- min(as.numeric(poss_zones[,3]));
					if (zlb >= redone_collections$min_ma[paldbc] && zub <= redone_collections$max_ma[paldbc])
						relevant_zones <- rbind(relevant_zones,poss_zones);
					}
				}
			}
		}	# end case where no rock unit is reported

	# now redate the collection based on zone or zones
	if (length(relevant_zones[,2])==0)	{
		match1 <- (1:zd)[zone_database$zone %in% redone_collections$zone[paldbc]];
		match2 <- (1:zd)[zone_database$zone_sr %in% redone_collections$zone[paldbc]];
		if (redone_collections$zone_species[paldbc]!="")	{
			match3 <- (1:zd)[zone_database$zone_species %in% redone_collections$zone_species[paldbc]];
			match4 <- (1:zd)[zone_database$zone_species_sr %in% redone_collections$zone_species[paldbc]];
			}	else	{
			match3 <- c();
			match4 <- c();
#			match3 <- match3[!zone_data$zone_species[match3] %in% ""];
			}
		matches <- sort(unique(c(match1,match2,match3,match4)));
		if (length(matches)==0)	{
			multizones <- diffindo_zone(zone_paleodb);
			mz <- length(multizones);
			if (mz>1)	{
				rel_z <- c();
	#			multimatches <- c();
				for (zz in 1:mz)	{
					# match zone names of zones to PaleoDB entered zones
					matches1 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					# match senior names of zones to PaleoDB entered zones
					matches2 <- (1:zd)[zone_database$zone_sr %in% multizones[zz]];
					# match species names only of zones to PaleoDB entered zones
					matches3 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
					matches4 <- (1:zd)[zone_database$zone_species_sr %in% multizones[zz]];
					matches5 <- (1:zd)[zone_database$zone %in% multizones[zz]];
					matches6 <- (1:zd)[zone_database$zone_species %in% multizones[zz]];
					matches7 <- (1:zd)[zone_database$non_taxon_zone_label %in% multizones[zz]];
					multimatches <- c(sort(unique(c(matches1,matches2,matches3,matches4,matches5,matches6,matches7))));
					if (length(multimatches)>0)	{
						rel_z <- rbind(rel_z,cbind(as.character(zone_database$zone_sr[multimatches]),as.numeric(zone_database$ma_lb[multimatches]),as.numeric(zone_database$ma_ub[multimatches]),as.character(zone_database$interval_lb[multimatches]),as.character(zone_database$interval_ub[multimatches])));
						matches <- c(matches,multimatches);
						}
					}
				} else	{
				rel_z <- c();
				}
			} else	{
			rel_z <- cbind(as.character(zone_database$zone_sr[matches]),as.numeric(zone_database$ma_lb[matches]),as.numeric(zone_database$ma_ub[matches]),as.character(zone_database$interval_lb[matches]),as.character(zone_database$interval_ub[matches]));
			}
		if (length(matches)>0)	{
			zone_taxa <- sort(unique(rel_z[,1]));
		# there should be only one taxon here.  If 2+ were found (e.g., "Ctenognathodus murchisoni",
		#	"Cyrtograptus murchisoni" and "Didymograptus murchisoni" for "murchisoni"), then find the appropriate one
			if (length(zone_taxa)>=1)	{
#			possible_zones <- c();
				for (ztx  in 1:length(zone_taxa))	{
					poss_zones <- subset(rel_z,rel_z[,1]==zone_taxa[ztx]);
					zlb <- max(as.numeric(poss_zones[,2]));
					zub <- min(as.numeric(poss_zones[,3]));
					if (zlb >= redone_collections$min_ma[paldbc] && zub <= redone_collections$max_ma[paldbc])
						relevant_zones <- rbind(relevant_zones,poss_zones);
					}
				}
			}
		}
	
	# now, make adjustments
	colnames(relevant_zones) <- c("zone_sr","ma_lb","ma_ub","interval_lb","interval_ub");
	if (length(relevant_zones$ma_lb)>0)	{
		zones_matched <- c(zones_matched,paldbc);
		# if zone starts after earliest possible age for formation, adjust oldest age & stage
		ma_lb <- as.numeric(as.character(relevant_zones$ma_lb));
		ma_mx <- max(ma_lb);
		lbz <- match(ma_mx,ma_lb);
		if (redone_collections$ma_lb[paldbc] > ma_mx || redone_collections$ma_lb[paldbc]==0)	{
			redone_collections$ma_lb[paldbc] <- ma_mx;
			redone_collections$interval_lb[paldbc] <- as.character(relevant_zones$interval_lb[lbz]);
			}
		# if zone ends before latest possible age for formation, adjust youngest age & stage
		ma_ub <- as.numeric(as.character(relevant_zones$ma_ub));
		ma_mn <- min(ma_ub);
		ubz <- match(ma_mn,ma_ub);
		if (redone_collections$ma_ub[paldbc] < ma_mn)	{
			redone_collections$ma_ub[paldbc] <- ma_mn;
			redone_collections$interval_ub[paldbc] <- as.character(relevant_zones$interval_ub[ubz]);
			}
		}
#	print(zones_matched);
	}

# check for poorly restricted rocks that had better ages from Paleodb
ma_max <- max(collections$max_ma)+25;
ma_min <- min(collections$min_ma)-25;
timey_wimey <- subset(time_scale,time_scale$ma_lb<=ma_max);
timey_wimey <- subset(timey_wimey,timey_wimey$ma_ub>=ma_min);
stage <- timey_wimey$interval;
timey_wimey$interval <- sapply(stage,turgio_stage);

disjunct <- c();
for (cn in 1:ncoll)	{
	if (redone_collections$max_ma[cn]>redone_collections$ma_ub[cn] && redone_collections$ma_lb[cn]>redone_collections$min_ma[cn])	{
		wibbly_wobbly <- accio_temporal_overlap(lb1=redone_collections$max_ma[cn],ub1=redone_collections$min_ma[cn],lb2=redone_collections$ma_lb[cn],ub2=redone_collections$ma_ub[cn]);
		if (wibbly_wobbly$ma_lb==redone_collections$max_ma[cn])
			redone_collections$interval_lb[cn] <- redone_collections$early_interval[cn];
		if (wibbly_wobbly$ma_ub==redone_collections$min_ma[cn])
			redone_collections$interval_ub[cn] <- redone_collections$late_interval[cn];
		redone_collections$ma_lb[cn] <- wibbly_wobbly$ma_lb;
		redone_collections$ma_ub[cn] <- wibbly_wobbly$ma_ub;
		} else	{
		disjunct <- c(disjunct,cn);
		}
	}

# cull out any collections using time units not in the external time scale
stage_numb_1 <- match(collections$early_interval,timey_wimey$interval);
stage_numb_2 <- match(collections$late_interval,timey_wimey$interval);
cut_these <- (1:ncoll)[is.na(stage_numb_1)];
if (length(cut_these) > 0)	{
	pass <- (1:rd_coll)[!is.na(stage_numb_1)];
	rdc <- rdc[pass,];
	stage_numb_1 <- stage_numb_1[pass];
	stage_numb_2 <- stage_numb_2[pass];
	early_interval <- early_interval[pass];
	late_interval <- late_interval[pass];
	double_check_extdb <- double_check_extdb[pass];
	rd_coll <- nrow(rdc);
	} else	{
	rd_coll <- ncoll;
	}
cut_these <- (1:rd_coll)[is.na(stage_numb_2)];
if (length(cut_these) > 0)	{
	pass <- (1:rd_coll)[!is.na(stage_numb_2)];
	rdc <- rdc[pass,];
	stage_numb_1 <- stage_numb_1[pass];
	stage_numb_2 <- stage_numb_2[pass];
	early_interval <- early_interval[pass];
	late_interval <- late_interval[pass];
	double_check_extdb <- double_check_extdb[pass];
	rd_coll <- nrow(rdc);
	}

missing <- collections_w_zones[!collections_w_zones %in% zones_matched];
untraceable_zone_information <- cbind(redone_collections$collection_no[missing],redone_collections$formation[missing],redone_collections$member[missing],redone_collections$zone[missing]);
colnames(untraceable_zone_information) <- c("collection_no","formation","member","zone");
untraceable_zone_information <- data.frame(untraceable_zone_information);
# 2019-03-07: routine to get unidentified rock units limited....
#ranges_paledb <- redone_collections$max_ma-redone_collections$min_ma;
#ranges_redone <- redone_collections$ma_lb-redone_collections$ma_ub;
#max_ma_downloaded <- max(redone_collections$max_ma);
#min_ma_downloaded <- min(redone_collections$min_ma);
#too_old <- (1:ncoll)[redone_collections$ma_lb > max_ma_downloaded];
#redone_collections$ma_lb[too_old] <- max_ma_downloaded;
#redone_collections$interval_lb[match(max_ma_downloaded,redone_collections$max_ma)]
#redone_collections$interval_lb
#too_young <- (1:ncoll)[redone_collections$ma_ub < min_ma_downloaded];

if (save_files)	{
	output1 <- paste(study,"_Recalibrated_Collections",output_type,sep="");
	output2 <- paste(study,"_Collection_Zone_Data_Untraced",output_type,sep="");
	if (directory!="")	{
		output1 <- paste(directory,output1,sep="");
		output2 <- paste(directory,output2,sep="");
		}
	if (output_type==".csv")	{
		write.csv(redone_collections,output1,row.names=FALSE);
		write.csv(untraceable_zone_information,output2,row.names=FALSE);
		} else	{
		write.table(redone_collections,output1,row.names=FALSE,col.names=TRUE,sep="\t");
		write.table(untraceable_zone_information,output2,row.names=FALSE,col.names=TRUE,sep="\t");
		}
	}

output <- list(redone_collections,untraceable_zone_information);
names(output) <- c("Recalibrated_Collections","Untraced_Collection_Zones");

return(output);
}

# link PaleoDB collections to external database based on rock units, zones, etc.
match_paleodb_collections_to_external_stratigraphic_database <- function(collections,wagner_rocks)	{

n_rocks <- nrow(wagner_rocks);
### quick kluge for John....
#entered_rock_and_onsets <-cbind(collections$formation,collections$member,collections$max_ma);
#turgio_rocks <- paleodb_editing_kluges(entered_rock_and_onsets);
#collections$formation <- turgio_rocks[,1];
#collections$member <- turgio_rocks[,2];
### some day I won't need this..... :-(
### today is that day! :-)

### Now, take PaleoDB collections and match them to external database
# separate collections with different types of stratigraphic information
n_coll <- nrow(collections);
paleodb_coll_w_formation <- (1:n_coll)[collections$formation!=""];
paleodb_coll_w_member <- (1:n_coll)[collections$member!=""];
paleodb_coll_w_group <- (1:n_coll)[collections$stratgroup!=""];
paleodb_coll_w_formation_and_member <- paleodb_coll_w_member[paleodb_coll_w_member %in% paleodb_coll_w_formation];
paleodb_coll_w_member_only <-  paleodb_coll_w_member[!paleodb_coll_w_member %in% paleodb_coll_w_formation];
paleodb_coll_w_formation_only <- paleodb_coll_w_formation[!paleodb_coll_w_formation %in% paleodb_coll_w_member];
paleodb_coll_w_group_only <- paleodb_coll_w_group[!paleodb_coll_w_group %in% c(paleodb_coll_w_formation,paleodb_coll_w_member)];
collections_w_rock_names <- sort(unique(c(paleodb_coll_w_formation,paleodb_coll_w_member,paleodb_coll_w_group)))
paleodb_coll_w_rock_names <- length(collections_w_rock_names);

paleodb_coll_entered_rock_unit <- array("",dim=n_coll);
paleodb_coll_entered_rock_unit[paleodb_coll_w_formation_only] <- collections$formation[paleodb_coll_w_formation_only];
paleodb_coll_entered_rock_unit[paleodb_coll_w_member_only] <- collections$member[paleodb_coll_w_member_only];
paleodb_coll_entered_rock_unit[paleodb_coll_w_formation_and_member] <- paste(as.character(collections$formation[paleodb_coll_w_formation_and_member])," (",as.character(collections$member[paleodb_coll_w_formation_and_member]),")",sep="");

### first pass: just check entered names (clean rock units basic)
paleodb_clean_formation_basic <- paleodb_clean_formation_no_rock <- paleodb_clean_formation_no_rock_formal <- collections$formation;
named_rock_unit <- paleodb_clean_formation_basic[paleodb_coll_w_formation];
paleodb_clean_formation_basic[paleodb_coll_w_formation] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal=FALSE);
paleodb_clean_formation_no_rock[paleodb_coll_w_formation] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal=FALSE);
paleodb_clean_formation_no_rock_formal[paleodb_coll_w_formation] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal=TRUE);

paleodb_clean_member_basic <- paleodb_clean_member_no_rock <- paleodb_clean_member_no_rock_formal <- collections$member;
named_rock_unit <- paleodb_clean_member_basic[paleodb_coll_w_member];
paleodb_clean_member_basic[paleodb_coll_w_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=FALSE,delete_informal=FALSE);
paleodb_clean_member_no_rock[paleodb_coll_w_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=FALSE);
paleodb_clean_member_no_rock_formal[paleodb_coll_w_member] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=TRUE);
# get lists of collections with members, members only, formations+members and formations only
#		after different cleaning
paleodb_coll_w_member_nr <- paleodb_coll_w_member[paleodb_clean_member_no_rock[paleodb_coll_w_member]!=""];
paleodb_coll_w_member_nr_f <- paleodb_coll_w_member[paleodb_clean_member_no_rock_formal[paleodb_coll_w_member]!=""];
paleodb_coll_w_member_only_nr <- paleodb_coll_w_member_nr[paleodb_clean_formation_no_rock[paleodb_coll_w_member_nr]==""];
paleodb_coll_w_member_only_nr_f <- paleodb_coll_w_member_nr_f[paleodb_clean_formation_no_rock_formal[paleodb_coll_w_member_nr_f]==""];
paleodb_coll_w_formation_and_member_nr <- paleodb_coll_w_member_nr[!paleodb_coll_w_member_nr %in% paleodb_coll_w_member_only_nr];
paleodb_coll_w_formation_and_member_nr_f <- paleodb_coll_w_member_nr_f[!paleodb_coll_w_member_nr_f %in% paleodb_coll_w_member_only_nr_f];
paleodb_coll_w_formation_only_nr <- paleodb_coll_w_formation[!paleodb_coll_w_formation %in% paleodb_coll_w_formation_and_member_nr];
paleodb_coll_w_formation_only_nr_f <- paleodb_coll_w_formation[!paleodb_coll_w_formation %in% paleodb_coll_w_formation_and_member_nr_f];

paleodb_coll_w_formation_or_member <- sort(unique(c(paleodb_coll_w_formation,paleodb_coll_w_member)));

paleodb_clean_rock_unit_basic <- paleodb_clean_formation_basic;
paleodb_clean_rock_unit_basic[paleodb_coll_w_member_only] <- paleodb_clean_member_basic[paleodb_coll_w_member_only];
paleodb_clean_rock_unit_basic[paleodb_coll_w_formation_and_member] <- paste(paleodb_clean_formation_basic[paleodb_coll_w_formation_and_member]," (",paleodb_clean_member_basic[paleodb_coll_w_formation_and_member],")",sep="");
paleodb_clean_rock_unit_no_rock <- paleodb_clean_formation_no_rock;
paleodb_clean_rock_unit_no_rock[paleodb_coll_w_member_only_nr] <- paleodb_clean_member_no_rock[paleodb_coll_w_member_only_nr];
paleodb_clean_rock_unit_no_rock[paleodb_coll_w_formation_and_member_nr] <- paste(paleodb_clean_formation_no_rock[paleodb_coll_w_formation_and_member_nr]," (",paleodb_clean_member_no_rock[paleodb_coll_w_formation_and_member_nr],")",sep="");
paleodb_clean_rock_unit_no_rock_formal <- paleodb_clean_formation_no_rock_formal;
paleodb_clean_rock_unit_no_rock_formal[paleodb_coll_w_member_only_nr_f] <- paleodb_clean_member_no_rock_formal[paleodb_coll_w_member_only_nr_f];
paleodb_clean_rock_unit_no_rock_formal[paleodb_coll_w_formation_and_member_nr_f] <- paste(paleodb_clean_formation_no_rock_formal[paleodb_coll_w_formation_and_member_nr_f]," (",paleodb_clean_member_no_rock_formal[paleodb_coll_w_formation_and_member_nr_f],")",sep="");

paleodb_clean_group_basic <- paleodb_clean_group_no_rock <- paleodb_clean_group_no_rock_formal <- collections$stratgroup;
named_rock_unit <- paleodb_clean_group_basic[paleodb_coll_w_group];
paleodb_clean_group_basic[paleodb_coll_w_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = FALSE,delete_informal=FALSE);
paleodb_clean_group_no_rock[paleodb_coll_w_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal=FALSE);
paleodb_clean_group_no_rock_formal[paleodb_coll_w_group] <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type = TRUE,delete_informal=TRUE);

# if exact same name entered into formation & member, then make it just formation name.
same_name <- paleodb_coll_w_formation_and_member[paleodb_clean_formation_basic[paleodb_coll_w_formation_and_member]==paleodb_clean_member_basic[paleodb_coll_w_formation_and_member]];
paleodb_clean_member_no_rock_formal[same_name] <- paleodb_clean_member_no_rock[same_name] <- paleodb_clean_member_basic[same_name]<- "";
paleodb_clean_rock_unit_no_rock_formal[same_name] <- paleodb_clean_formation_no_rock_formal[same_name];
paleodb_clean_rock_unit_no_rock[same_name] <- paleodb_clean_formation_no_rock[same_name];
paleodb_clean_rock_unit_basic[same_name] <- paleodb_clean_formation_basic[same_name];

# set up new fields to be added
rock_unit_senior2 <- rock_unit_senior <- interval_lb <- interval_ub <- array("?",dim=n_coll);
ma_lb <- ma_ub <- rock_no <- rock_no_sr <- formation_no <- rock2_no <- rock2_no_sr <- formation2_no <- clean_match <- array(0,dim=n_coll);

# test 1: does name match any of the basic full names?
test_cols <- match(c("full_name","rock_unit_clean_basic","rock_unit_clean_no_rock"),colnames(wagner_rocks));
examined_matrix <- wagner_rocks[,test_cols];
for (i in 1:ncol(examined_matrix))	examined_matrix[,i] <- tolower(examined_matrix[,i]);

#find_me <- "Haverford";
find_me <- paleodb_clean_rock_unit_basic[paleodb_coll_w_formation_or_member];
external_rock_rows <- sapply(tolower(find_me),matrix_row_match,examined_matrix);
find_me <- paleodb_clean_rock_unit_no_rock[paleodb_coll_w_formation_or_member];
e_r_r_2 <- sapply(tolower(find_me),matrix_row_match,examined_matrix);
find_me <- paleodb_clean_rock_unit_no_rock_formal[paleodb_coll_w_formation_or_member];
e_r_r_3 <- sapply(tolower(find_me),matrix_row_match,examined_matrix);
external_rock_rows[external_rock_rows < 1] <- e_r_r_2[external_rock_rows < 1];
external_rock_rows[external_rock_rows < 1] <- e_r_r_3[external_rock_rows < 1];
# positive integers in testing are row numbers from stratigraphic database for matching rocks
# -1 means that no matches were made.
# -2 means that 2+ matches were made.
# paleodb_coll_w_formation_or_member[testing>0] gives row numbers of paleodb collections for which we've matched rocks
emended_rocks <- external_rock_rows[external_rock_rows>0];
updated_collections <- paleodb_coll_w_formation_or_member[external_rock_rows>0];
ma_lb[updated_collections] <- wagner_rocks$ma_lb[external_rock_rows[external_rock_rows>0]];
ma_ub[updated_collections] <- wagner_rocks$ma_ub[external_rock_rows[external_rock_rows>0]];
rock_no[updated_collections] <- wagner_rocks$rock_no[external_rock_rows[external_rock_rows>0]];
rock_no_sr[updated_collections] <- wagner_rocks$rock_no_sr[external_rock_rows[external_rock_rows>0]];
formation_no[updated_collections] <- wagner_rocks$formation_no[external_rock_rows[external_rock_rows>0]];
rock_unit_senior[updated_collections] <- wagner_rocks$rock_unit_senior[external_rock_rows[external_rock_rows>0]];
interval_lb[updated_collections] <- wagner_rocks$interval_lb[external_rock_rows[external_rock_rows>0]];
interval_ub[updated_collections] <- wagner_rocks$interval_ub[external_rock_rows[external_rock_rows>0]];
clean_match[updated_collections] <- 1;

# now, deal with stubborn members & formations....
peles <- sort(c(paleodb_coll_w_formation_only,paleodb_coll_w_member_only));	# one name units
need_info_still <- peles[ma_lb[peles]==0];
n_i_p <- length(need_info_still);
test_cols_p <- match(c("group","formation","member","group_clean_basic","formation_clean_basic","member_clean_basic","group_clean_no_rock","formation_clean_no_rock","member_clean_no_rock","group_clean_no_rock_formal","formation_clean_no_rock_formal","member_clean_no_rock_formal"),colnames(wagner_rocks));
examined_matrix_p <- wagner_rocks[,test_cols_p];
for (i in 1:length(test_cols_p))	examined_matrix_p[,i] <- tolower(examined_matrix_p[,i]);

test_cols_g <- match(c("group","group_clean_basic","group_clean_no_rock","group_clean_no_rock_formal"),colnames(wagner_rocks));
examined_matrix_g <- wagner_rocks[,test_cols_g];
for (i in 1:length(test_cols_g))	examined_matrix_g[,i] <- tolower(examined_matrix_g[,i]);
#matrix_row_match("el paso",examined_matrix_p)

find_me <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
external_rock_rows_p <- sapply(find_me,matrix_row_match,examined_matrix_p);
find_me <- tolower(paleodb_clean_rock_unit_no_rock[need_info_still]);
e_r_r_p_2 <- sapply(find_me,matrix_row_match,examined_matrix_p);
e_r_r_p_2[find_me==""] <- -1;
external_rock_rows_p[external_rock_rows_p<1] <- e_r_r_p_2[external_rock_rows_p<1];
find_me <- tolower(paleodb_clean_rock_unit_no_rock_formal[need_info_still]);
e_r_r_p_3 <- sapply(find_me,matrix_row_match,examined_matrix_p);
e_r_r_p_3[find_me==""] <- -1;
external_rock_rows_p[external_rock_rows_p<1] <- e_r_r_p_3[external_rock_rows_p<1];
#cbind(need_info_still,external_rock_rows_p)
# collections with members with 2+ ids in external database.
#still_need <- (1:length(need_info_member))[external_rock_rows_m==-2];
for_homer <- (1:n_i_p)[external_rock_rows_p==-2];
two_plus <- need_info_still[external_rock_rows_p==-2];
find_me_2 <- tolower(paleodb_clean_rock_unit_no_rock_formal[two_plus]);
two_plus_names <- tolower(unique(find_me_2));
# use try_again again!
tp <- 0;
while (tp < length(two_plus_names))	{
#for (tp in 1:length(two_plus_names))	{
	tp <- tp+1;
	xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
	if (nrow(xxx)>0)	{
		yyy <- unique(xxx[,1]);
		if (length(unique(wagner_rocks$rock_no_sr[yyy]))==1)	{
			# gets cases where rock unit is a member in external data base
			# get needy collections with this as member
			tried_m <- (1:n_i_p)[tolower(paleodb_clean_member_no_rock_formal[need_info_still])==two_plus_names[tp]];
			tried_f <- (1:n_i_p)[tolower(paleodb_clean_formation_no_rock_formal[need_info_still])==two_plus_names[tp]];
			tried <- unique(c(tried_m,tried_f));
			external_rock_rows_p[tried] <- match(wagner_rocks$rock_no_sr[yyy[1]],wagner_rocks$rock_no)
			} else if (length(unique(wagner_rocks$formation_no[yyy]))==1)	{
			# get needy collections with this as formation
			# gets cases where rock unit is a formation in external data base
			tried <- (1:n_i_p)[tolower(paleodb_clean_formation_no_rock_formal[need_info_still])==two_plus_names[tp]];
			external_rock_rows_p[tried] <- match(wagner_rocks$formation_no[yyy[1]],wagner_rocks$rock_no);
			} else	{
			vvv <- which(examined_matrix_g==two_plus_names[tp],arr.ind=TRUE);
			uuu <- unique(vvv[,1]);
#			if (length(uuu[wagner_rocks$formation[uuu]=="•"])==1)	{
			if (length(uuu[wagner_rocks$formation[uuu]==""])==1)	{
				# gets cases where rock unit is a group in external data base
				tried <- (1:n_i_p)[find_me %in% two_plus_names[tp]];
#				tried <- (1:n_i_p)[tolower(paleodb_clean_group_no_rock_formal[need_info_still])==two_plus_names[tp]];
				external_rock_rows_p[tried] <- match(uuu[wagner_rocks$formation[uuu]==""],wagner_rocks$rock_no);
#				external_rock_rows_p[tried] <- match(uuu[wagner_rocks$formation[uuu]=="•"],wagner_rocks$rock_no);
				}
			}
		}
#	print(cbind(paleodb_clean_rock_unit_basic[two_plus],external_rock_rows_p[for_homer]))
	}

emended_rocks <- (1:n_i_p)[external_rock_rows_p>0];
emended_rocks <- emended_rocks[!is.na(emended_rocks)];	# KLUGE!!!!
updated_collections <- need_info_still[emended_rocks];
if (length(updated_collections) > 0)	{
	ma_lb[updated_collections] <- wagner_rocks$ma_lb[external_rock_rows_p[emended_rocks]];
	ma_ub[updated_collections] <- wagner_rocks$ma_ub[external_rock_rows_p[emended_rocks]];
	rock_no[updated_collections] <- wagner_rocks$rock_no[external_rock_rows_p[emended_rocks]];
	rock_no_sr[updated_collections] <- wagner_rocks$rock_no_sr[external_rock_rows_p[emended_rocks]];
	formation_no[updated_collections] <- wagner_rocks$formation_no[external_rock_rows_p[emended_rocks]];
	#rock_unit_senior[updated_collections] <- as.character(wagner_rocks$rock_unit_senior[emended_rocks]);
	rock_unit_senior[updated_collections] <- wagner_rocks$rock_unit_senior[external_rock_rows_p[emended_rocks]];
	interval_lb[updated_collections] <- wagner_rocks$interval_lb[external_rock_rows_p[emended_rocks]];
	interval_ub[updated_collections] <- wagner_rocks$interval_ub[external_rock_rows_p[emended_rocks]];
	clean_match[updated_collections] <- 1;
	}
### done with stubborn rocks

# look at things with just group ids
need_info_still <- collections_w_rock_names[ma_lb[collections_w_rock_names]==0];
need_info_still_f <- need_info_still[paleodb_clean_formation_basic[need_info_still]!=""];
need_info_still <- need_info_still[!need_info_still %in% need_info_still_f];
n_i_s <- length(need_info_still);
find_me <- tolower(paleodb_clean_group_no_rock_formal[need_info_still]);
external_rock_rows <- sapply(find_me,matrix_row_match,examined_matrix_p);
#cbind(paleodb_clean_group_no_rock_formal[need_info_still_g],external_rock_rows_g);
for_bart <- (1:n_i_s)[external_rock_rows==-2];
two_plus <- need_info_still[external_rock_rows==-2];
find_me <- tolower(paleodb_clean_group_no_rock_formal[two_plus]);
two_plus_names <- tolower(unique(find_me));
# use try_again again!
tp <- 0;
while (tp < length(two_plus_names))	{
#for (tp in 1:length(two_plus_names))	{
	tp <- tp+1;
	xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
	yyy <- unique(xxx[,1]);
	vvv <- c();
#	if (length(yyy[tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])=="• (¢)"])==1)	{
#		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])=="• (¢)")];
	if (length(yyy[tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==""])==1)	{
		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])=="")];
		}	else if (length(yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])])==1)	{
		# failing that, look to see if the group name is a stand-alone formation name
		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])];
		} else if (length(yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])])>1)	{
		vvv <- yyy[wagner_rocks$member[yyy]==""];
		}
	# get all of the collections with this rock unit
	if (length(vvv)==1)	{
		tried <- (1:n_i_s)[tolower(paleodb_clean_group_no_rock_formal[need_info_still])==two_plus_names[tp]];
		external_rock_rows[tried] <- vvv[1];
		}
#	print(cbind(paleodb_clean_group_basic[two_plus],external_rock_rows_g[for_bart]));
	}

emended_rocks <- (1:n_i_s)[external_rock_rows>0];
emended_rocks <- emended_rocks[!is.na(emended_rocks)];		# KLUGE!!!!!
updated_collections <- need_info_still[emended_rocks];
if (length(updated_collections)>0)	{
	ma_lb[updated_collections] <- wagner_rocks$ma_lb[external_rock_rows[emended_rocks]];
	ma_ub[updated_collections] <- wagner_rocks$ma_ub[external_rock_rows[emended_rocks]];
	rock_no[updated_collections] <- wagner_rocks$rock_no[external_rock_rows[emended_rocks]];
	rock_no_sr[updated_collections] <- wagner_rocks$rock_no_sr[external_rock_rows[emended_rocks]];
	formation_no[updated_collections] <- wagner_rocks$formation_no[external_rock_rows[emended_rocks]];
	#rock_unit_senior[updated_collections] <- as.character(wagner_rocks$rock_unit_senior[emended_rocks]);
	rock_unit_senior[updated_collections] <- wagner_rocks$rock_unit_senior[external_rock_rows[emended_rocks]];
	interval_lb[updated_collections] <- wagner_rocks$interval_lb[external_rock_rows[emended_rocks]];
	interval_ub[updated_collections] <- wagner_rocks$interval_ub[external_rock_rows[emended_rocks]];
	clean_match[updated_collections] <- 1;
	}
# we are getting near the point of just burning it all and making it look like an electrical thing
need_info_still <- paleodb_coll_w_formation_or_member[ma_lb[paleodb_coll_w_formation_or_member]==0];
n_i_s <- length(need_info_still);
find_me <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
external_rock_rows <- sapply(find_me,matrix_row_match,examined_matrix);
for_bender <- (1:n_i_s)[external_rock_rows==-2];
two_plus <- need_info_still[external_rock_rows==-2];
find_me <- tolower(paleodb_clean_rock_unit_basic[two_plus]);
two_plus_names <- tolower(unique(find_me));
tp <- 0;
#for (tp in 1:length(two_plus_names))	{
while (tp < length(two_plus_names))	{
	tp <- tp+1;
	xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
	yyy <- unique(xxx[,1]);
	vvv <- c();
	if (length(yyy[tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==""])==1)	{
		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])=="")];
		}	else if (length(yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])])==1)	{
		# failing that, look to see if the group name is a stand-alone formation name
		vvv  <- yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])];
		} else if (length(yyy[(tolower(wagner_rocks$rock_unit_clean_no_rock_formal[yyy])==two_plus_names[tp])])>1)	{
		vvv <- yyy[wagner_rocks$member[yyy]==""];
		}
	# get all the collections with this rock unit
	if (length(vvv)==1)	{
		tried <- (1:n_i_s)[tolower(paleodb_clean_rock_unit_no_rock_formal[need_info_still])==two_plus_names[tp]];
		external_rock_rows[tried] <- vvv[1];
		}
	}

emended_rocks <- (1:n_i_s)[external_rock_rows>0];
updated_collections <- need_info_still[emended_rocks];
if (length(updated_collections)>0)	{
	ma_lb[updated_collections] <- wagner_rocks$ma_lb[external_rock_rows[emended_rocks]];
	ma_ub[updated_collections] <- wagner_rocks$ma_ub[external_rock_rows[emended_rocks]];
	rock_no[updated_collections] <- wagner_rocks$rock_no[external_rock_rows[emended_rocks]];
	rock_no_sr[updated_collections] <- wagner_rocks$rock_no_sr[external_rock_rows[emended_rocks]];
	formation_no[updated_collections] <- wagner_rocks$formation_no[external_rock_rows[emended_rocks]];
	#rock_unit_senior[updated_collections] <- as.character(wagner_rocks$rock_unit_senior[emended_rocks]);
	rock_unit_senior[updated_collections] <- wagner_rocks$rock_unit_senior[external_rock_rows[emended_rocks]];
	interval_lb[updated_collections] <- wagner_rocks$interval_lb[external_rock_rows[emended_rocks]];
	interval_ub[updated_collections] <- wagner_rocks$interval_ub[external_rock_rows[emended_rocks]];
	clean_match[updated_collections] <- 0;
	}

# I am getting stabby, so let's cut up some units with 2+ names.....
need_info_still <- paleodb_coll_w_formation_or_member[ma_lb[paleodb_coll_w_formation_or_member]==0];
n_i_s <- length(need_info_still);
find_me <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
#funky_text <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
#find_me <- sapply(funky_text,transmogrify_diacritics);
cc <- 0;
#for (cc in 1:n_i_s)	{
while (cc < n_i_s)	{
	cc <- cc+1;
	named_rock_unit <- collections$formation[need_info_still[cc]];
	named_rock_unit <- diffindo_lumped_rock_units(named_rock_unit);
	if (length(named_rock_unit)>1)	{
		named_rock_unit <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=TRUE);
		find_me2 <- tolower(named_rock_unit);
		external_rock_rows <- sapply(find_me2,matrix_row_match,examined_matrix_p);
		if (sum(external_rock_rows==-2)>0)	{
			two_plus_names <- find_me2[external_rock_rows==-2];
			for (tp in 1:length(two_plus_names))	{
				xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
				yyy <- unique(xxx[,1]);
				zzz <- yyy[wagner_rocks$rock_no[yyy] %in% wagner_rocks$formation_no[yyy]];
				if (length(zzz)==0)
					zzz <- yyy[wagner_rocks$rock_no[yyy] %in% wagner_rocks$rock_no_sr[yyy]];
				if (length(zzz)==1)	{
					external_rock_rows[match(two_plus_names[tp],find_me2)] <- zzz[1];
					} else	{
					vvv <- yyy[wagner_rocks$formation[yyy]=="" && wagner_rocks$member[yyy]==""];
					if (length(vvv)==1)	{
						external_rock_rows[match(two_plus_names[tp],find_me2)] <- vvv[1];
						} else	{
						uuu <- yyy*(wagner_rocks$formation[yyy]==wagner_rocks$group[yyy]) * (wagner_rocks$member[yyy]=="");
						uuu <- uuu[uuu>1];
						if (length(uuu)==1)	{
							external_rock_rows[match(two_plus_names[tp],find_me2)] <- uuu[1];
							}
						}
					}
				}
			}
		if (sum(external_rock_rows>0)==1)	{
			# do a thing
			external_rock_rows <- external_rock_rows[external_rock_rows>0];
			rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[1]]
			rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[1]]
			formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows[1]]
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[1]],wagner_rocks$rock_no)];
			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[1]];
			interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows[1]];
			interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows[1]];
			ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows[1]];
			ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows[1]];
			clean_match[need_info_still[cc]] <- 1;
			} else if (sum(external_rock_rows>0)>1)	{
			# do another thing
			external_rock_rows <- external_rock_rows[external_rock_rows>0];
			if (length(external_rock_rows)>2)	{
				e_r_r <- 1:length(external_rock_rows);
				aa <- match(max(wagner_rocks$ma_lb[external_rock_rows]),wagner_rocks$ma_lb[external_rock_rows]);
				r_e_r_r <- e_r_r[!e_r_r %in% aa];
				zz <- r_e_r_r[match(min(wagner_rocks$ma_ub[external_rock_rows[r_e_r_r]]),wagner_rocks$ma_ub[external_rock_rows[r_e_r_r]])];
				mm <- r_e_r_r[!r_e_r_r %in% zz];
				external_rock_rows <- external_rock_rows[c(aa,zz,mm)];
				}
			ma_lb[need_info_still[cc]] <- max(wagner_rocks$ma_lb[external_rock_rows[1:2]]);
			ma_ub[need_info_still[cc]] <- min(wagner_rocks$ma_ub[external_rock_rows[1:2]]);
			rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[1]];
			rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[1]];
			formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows[1]];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[1]];
			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[1]];
			rock2_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[2]];
			rock2_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[2]];
			formation2_no[need_info_still[cc]] <-  wagner_rocks$formation_no[external_rock_rows[2]];
#			rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[2]];
#			rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[2]],wagner_rocks$rock_no)];
			rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[2]];
			interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows[match(max(wagner_rocks$ma_lb[external_rock_rows[1:2]]),wagner_rocks$ma_lb[external_rock_rows[1:2]])]];
			interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows[match(min(wagner_rocks$ma_ub[external_rock_rows[1:2]]),wagner_rocks$ma_ub[external_rock_rows[1:2]])]];
			clean_match[need_info_still[cc]] <- 1;
			}
		} else	{
		named_rock_unit <- collections$member[need_info_still[cc]];
		named_rock_unit <- diffindo_lumped_rock_units(named_rock_unit);
		if (length(named_rock_unit)>1)	{
			named_rock_unit <- sapply(named_rock_unit,scourgify_rock_unit_names,dehyphenate=TRUE,delete_rock_type=TRUE,delete_informal=TRUE);
			find_me2 <- tolower(named_rock_unit);
			external_rock_rows <- sapply(find_me2,matrix_row_match,examined_matrix_p);
			if (sum(external_rock_rows==-2)>0)	{
				two_plus_names <- find_me2[external_rock_rows==-2];
				for (tp in 1:length(two_plus_names))	{
					xxx <- which(examined_matrix_p==two_plus_names[tp],arr.ind=TRUE);
					yyy <- unique(xxx[,1]);
					zzz <- yyy[wagner_rocks$rock_no[yyy] %in% wagner_rocks$formation_no[yyy]];
					if (length(zzz)==0)
						zzz <- yyy[wagner_rocks$rock_no[yyy] %in% wagner_rocks$rock_no_sr[yyy]];
					if (length(zzz)==1)	{
						external_rock_rows[match(two_plus_names[tp],find_me2)] <- zzz[1];
						} else	{
						vvv <- yyy[wagner_rocks$formation[yyy]=="" && wagner_rocks$member[yyy]==""];
						if (length(vvv)==1)	{
							external_rock_rows[match(two_plus_names[tp],find_me2)] <- vvv[1];
							} else	{
							uuu <- yyy*(wagner_rocks$formation[yyy]==wagner_rocks$group[yyy]) * (wagner_rocks$member[yyy]=="");
							uuu <- uuu[uuu>1];
							if (length(uuu)==1)	{
								external_rock_rows[match(two_plus_names[tp],find_me2)] <- uuu[1];
								}
							}
						}
					}
				}
			if (sum(external_rock_rows>0)==1)	{
				# do a thing
				external_rock_rows <- external_rock_rows[external_rock_rows>0];
				ma_lb[need_info_still[cc]] <- max(wagner_rocks$ma_lb[external_rock_rows[1]]);
				ma_ub[need_info_still[cc]] <- min(wagner_rocks$ma_ub[external_rock_rows[1]]);
				rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[1]]
				rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[1]]
				formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows[1]]
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[1]];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[1]],wagner_rocks$rock_no)];
				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[1]];
				interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows[1]];
				interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows[1]];
				clean_match[need_info_still[cc]] <- 1;
				} else if (sum(external_rock_rows>0)>1)	{
				# do another thing
				external_rock_rows <- external_rock_rows[external_rock_rows>0];
				if (length(external_rock_rows)>2)	{
					e_r_r <- 1:length(external_rock_rows);
					aa <- match(max(wagner_rocks$ma_lb[external_rock_rows]),wagner_rocks$ma_lb[external_rock_rows]);
					r_e_r_r <- e_r_r[!e_r_r %in% aa];
					zz <- r_e_r_r[match(min(wagner_rocks$ma_ub[external_rock_rows[r_e_r_r]]),wagner_rocks$ma_ub[external_rock_rows[r_e_r_r]])];
					mm <- r_e_r_r[!r_e_r_r %in% zz];
					external_rock_rows <- external_rock_rows[c(aa,zz,mm)];
					}
				ma_lb[need_info_still[cc]] <- max(wagner_rocks$ma_lb[external_rock_rows[1:2]]);
				ma_ub[need_info_still[cc]] <- min(wagner_rocks$ma_ub[external_rock_rows[1:2]]);
				rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[1]];
				rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[1]];
				formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows[1]];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[1]];
				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[1]],wagner_rocks$rock_no)];
				rock2_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows[2]];
				rock2_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows[2]];
				formation2_no[need_info_still[cc]] <-  wagner_rocks$formation_no[external_rock_rows[2]];
#				rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows[2]];
#				rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows[2]],wagner_rocks$rock_no)];
				rock_unit_senior2[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows[2]];
				interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows[match(max(wagner_rocks$ma_lb[external_rock_rows[1:2]]),wagner_rocks$ma_lb[external_rock_rows[1:2]])]];
				interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows[match(min(wagner_rocks$ma_ub[external_rock_rows[1:2]]),wagner_rocks$ma_ub[external_rock_rows[1:2]])]];
				clean_match[need_info_still[cc]] <- 1;
				}
			}
		}
	}
#cbind(need_info_still,collections$collection_no[need_info_still],paleodb_clean_rock_unit_basic[need_info_still],ma_lb[need_info_still])

#### last of it ???
pcwfom <- length(paleodb_coll_w_formation_or_member);
#(1:pcwfom)[is.na(ma_lb[paleodb_coll_w_formation_or_member])]
need_info_still <- paleodb_coll_w_formation_or_member[ma_lb[paleodb_coll_w_formation_or_member]==0];
n_i_s <- length(need_info_still);
find_me <- tolower(paleodb_clean_rock_unit_basic[need_info_still]);
#for (cc in 1:n_i_s)	{
cc <- 0;
while (cc<n_i_s)	{
	cc <- cc+1;
	if (paleodb_clean_formation_no_rock_formal[need_info_still[cc]]!="" && paleodb_clean_member_no_rock_formal[need_info_still[cc]]!="")	{
		find_me2 <- paleodb_clean_member_no_rock_formal[need_info_still[cc]];
		external_rock_rows_m <- sapply(tolower(find_me2),matrix_row_match,examined_matrix_p);
		find_me3 <- paleodb_clean_formation_no_rock_formal[need_info_still[cc]];
		external_rock_rows_f <- sapply(tolower(find_me3),matrix_row_match,examined_matrix_p);
		if (external_rock_rows_f>0 && external_rock_rows_m==-1)	{
#			print("One possible formation");
			external_rock_rows <- external_rock_rows_f;
			rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows];
			rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows];
			formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows],wagner_rocks$rock_no)];
			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows];
			interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows];
			interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows];
			ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows];
			ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows];
			} else if (external_rock_rows_f==-2 && external_rock_rows_m==-1)	{
#			print("Multiple possible formations");
			xxx <- which(examined_matrix_p==tolower(find_me3),arr.ind=TRUE);
			yyy <- unique(xxx[,1]);
			yyy <- yyy[wagner_rocks$rock_no[yyy]==wagner_rocks$formation_no[yyy]];
			if (length(yyy)==0)
				yyy <- yyy[wagner_rocks$rock_no_sr[yyy]==wagner_rocks$formation_no[yyy]];
			if (length(yyy)>1)
				yyy <- yyy[tolower(wagner_rocks$formation[yyy])==tolower(find_me3)]
			if (length(yyy)==1)	{
				external_rock_rows <- yyy[1];
				rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows];
				rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows];
				formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows],wagner_rocks$rock_no)];
				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows];
				interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows];
				interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows];
				ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows];
				ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows];
				}
			} else if (external_rock_rows_m>0)	{
#			print("One possible member");
			external_rock_rows <- external_rock_rows_m;
			rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows];
			rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows];
			formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows];
#			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows],wagner_rocks$rock_no)];
			rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows];
			interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows];
			interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows];
			ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows];
			ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows];
			} else if (external_rock_rows_m==-2)	{
#			print("Multiple possible members");
			xxx <- which(examined_matrix_p==tolower(find_me2),arr.ind=TRUE);
			yyy <- unique(xxx[,1]);
			yyy <- yyy[wagner_rocks$rock_no[yyy]==wagner_rocks$wagner_rocks$rock_no_sr[yyy]];
			if (length(yyy)==0)	{
				yyy <- unique(xxx[,1]);
				yyy <- yyy[wagner_rocks$rock_unit_clean_no_rock_formal[yyy]==find_me2];
				}
			if (length(yyy)>1)
				yyy <- yyy[tolower(wagner_rocks$member[yyy])==tolower(find_me2)]
			if (length(yyy)==1)	{
				external_rock_rows <- yyy[1];
				rock_no[need_info_still[cc]] <- wagner_rocks$rock_no[external_rock_rows];
				rock_no_sr[need_info_still[cc]] <- wagner_rocks$rock_no_sr[external_rock_rows];
				formation_no[need_info_still[cc]] <- wagner_rocks$formation_no[external_rock_rows];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[external_rock_rows];
#				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$full_name[match(wagner_rocks$rock_no_sr[external_rock_rows],wagner_rocks$rock_no)];
				rock_unit_senior[need_info_still[cc]] <- wagner_rocks$rock_unit_senior[external_rock_rows];
				interval_lb[need_info_still[cc]] <- wagner_rocks$interval_lb[external_rock_rows];
				interval_ub[need_info_still[cc]] <- wagner_rocks$interval_ub[external_rock_rows];
				ma_lb[need_info_still[cc]] <- wagner_rocks$ma_lb[external_rock_rows];
				ma_ub[need_info_still[cc]] <- wagner_rocks$ma_ub[external_rock_rows];
				}
			}
#		print(cbind(find_me[1:cc],ma_lb[need_info_still[1:cc]]));
		}
	}

# I've had it: I'm going Arya on any last undated collections!!!!
ma_lb[ma_lb==0] <- as.numeric(collections$max_ma[ma_lb==0]);
ma_ub[ma_ub==0] <- as.numeric(collections$min_ma[ma_ub==0]);
interval_lb[interval_lb=="?"] <- as.character(collections$early_interval[interval_lb=="?"]);
interval_ub[interval_ub=="?"] <- as.character(collections$late_interval[interval_ub=="?"]);
redone_collections <- cbind(collections,rock_no_sr,rock_no,formation_no,rock_unit_senior,rock2_no_sr,rock2_no,formation2_no,rock_unit_senior2,ma_lb,ma_ub,interval_lb,interval_ub,clean_match);
return(redone_collections);
}

# refine PaleoDB dates given finer time scale information than the PaleoDB uses.
redate_paleodb_collections_with_time_scale <- function(paleodb_collections,time_scale,zone_database)	{
# paleodb_collections: dataframe where:
#	paleodb_collections$early_interval gives oldest possible interval
#	paleodb_collections$late_interval gives youngest possible interval (this is filled if blank)
#	paleodb_collections$max_ma gives early_interal onset (initially from PaleoDB, then from Gradstein)
#	paleodb_collections$min_ma gives late_interal end (initially from PaleoDB, then from Gradstein)
# time_scale: dataframe where:
#	time_scale$interval gives interval name
#	time_scale$ma_lb gives interval onset (given Gradstein et al. 2012)
#	time_scale$ma_ub gives interval end (given Gradstein et al. 2012)
# zone_database: dataframe where:
#	we search for the entered interval names in one of several versions of the zone name
#	zone_database$interval_lb gives onset interval
#	zone_database$interval_ub gives end interval
#	zone_database$ma_lb gives interval onset (given Gradstein et al. 2012)
#	zone_database$ma_ub gives interval end (given Gradstein et al. 2012)

ncolls <- nrow(paleodb_collections);
no_late <- (1:ncolls)[paleodb_collections$late_interval==""];
paleodb_collections$late_interval[no_late] <- paleodb_collections$early_interval[no_late];
early_intervals <- match(paleodb_collections$early_interval,time_scale$interval);
late_intervals <- match(paleodb_collections$late_interval,time_scale$interval);
#paste(paleodb_collections$collection_no[paleodb_collections$early_interval %in% "Freboldi"],collapse=",");
if (sum(is.na(early_intervals))>0 || sum(is.na(late_intervals))>0)	{
	poss_zones <- sort(unique(c(paleodb_collections$early_interval[is.na(early_intervals)],paleodb_collections$late_intervals[is.na(late_intervals)])));
	trouble <- trouble_pz <- c();
	for (pz in 1:length(poss_zones))	{
#		pz <- pz+1;
		xxx <- unique(c(which(zone_database==poss_zones[pz],arr.ind = T)[,1],which(zone_database==tolower(poss_zones[pz]),arr.ind = T)[,1]));
		this_trouble <- c();
		if (length(xxx)==1)	{
			paleodb_collections$early_interval[paleodb_collections$early_interval %in% poss_zones[pz]] <- zone_database$interval_lb[xxx];
			paleodb_collections$late_interval[paleodb_collections$late_interval %in% poss_zones[pz]] <- zone_database$interval_ub[xxx];
			} else if (length(xxx)>1)	{
			if(length(unique(zone_database$interval_lb[xxx]))==1)	{
				temp_zones <- paleodb_collections$zone[paleodb_collections$early_interval %in% poss_zones[pz]]
				temp_zones[temp_zones==""] <- zone_database$zone_sr[xxx[1]];
				paleodb_collections$zone[paleodb_collections$early_interval %in% poss_zones[pz]] <- temp_zones;
				paleodb_collections$early_interval[paleodb_collections$early_interval %in% poss_zones[pz]] <- zone_database$interval_lb[xxx[1]];
				} else	{
#				print(paste("eff me down",pz));
				#zone_database$zone[xxx];
				this_trouble <- c(this_trouble,xxx);
				}
			if(length(unique(zone_database$interval_ub[xxx]))==1)	{
				paleodb_collections$late_interval[paleodb_collections$late_interval %in% poss_zones[pz]] <- zone_database$interval_ub[xxx[1]];
#				paleodb_collections$min_ma[paleodb_collections$late_interval %in% poss_zones[pz]];
				} else	{
#				print(paste("eff me up",pz));
				#print(zone_database[xxx,]);
				this_trouble <- c(this_trouble,xxx);
				}
			}
		trouble<- c(trouble,(unique(this_trouble)));
		trouble_pz <- c(trouble_pz,rep(pz,length(unique(this_trouble))));
		}
#	zone_database[trouble,];
	}
early_intervals <- match(paleodb_collections$early_interval,time_scale$interval);
late_intervals <- match(paleodb_collections$late_interval,time_scale$interval);
#old_maxes <- collections$max_ma
#old_mins <- collections$min_ma
# sum(is.na(late_intervals))
#(1:nrow(paleodb_collections))[is.na(late_intervals)]
#paleodb_collections$early_interval[(1:nrow(paleodb_collections))[is.na(late_intervals)]]
#paleodb_collections$late_interval[(1:nrow(paleodb_collections))[is.na(late_intervals)]]
paleodb_collections$max_ma[!is.na(early_intervals)] <- time_scale$ma_lb[early_intervals[!is.na(early_intervals)]];
paleodb_collections$min_ma[!is.na(late_intervals)] <- time_scale$ma_ub[late_intervals[!is.na(late_intervals)]];
misentered <- (1:ncolls)[paleodb_collections$max_ma==paleodb_collections$min_ma];
if (length(misentered)>0)	{
	dummy <- paleodb_collections$early_interval[misentered];
	paleodb_collections$early_interval[misentered] <- paleodb_collections$late_interval[misentered];
	paleodb_collections$late_interval[misentered] <- dummy;
	paleodb_collections$max_ma[misentered] <- time_scale$ma_lb[match(paleodb_collections$early_interval[misentered],time_scale$interval)];
	paleodb_collections$min_ma[misentered] <- time_scale$ma_ub[match(paleodb_collections$late_interval[misentered],time_scale$interval)];
	}
return(paleodb_collections);
}

# match zones from numerous collections to an external database
match_collections_zones_to_zone_database <- function(collections,zone_database)	{
# reduce zone database to feasible zones
#zone_database <- subset(zone_database,as.numeric(zone_database$ma_lb)<(max(collections$max_ma)+25));
#zone_database <- subset(zone_database,as.numeric(zone_database$ma_ub)>(min(collections$min_ma)-25));

if (is.na(match("zone_epithet",colnames(zone_database))))	{
	taxon_name <- zone <- as.character(zone_database$zone);
	zone_epithet <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
	zone <- as.character(zone_database$zone_sr);
	zone_epithet_sr <- sapply(zone,transmogrify_full_zone_names_to_species_names_only);
	zone_database <- cbind(zone_database,zone_epithet,zone_epithet_sr);
	}

ncolls <- nrow(collections);
colls_w_zones <- (1:ncolls)[collections$zone!=""];
zone <- collections$zone[colls_w_zones];
zone <- sapply(zone,turgio_zone);
zone_matches <- rep("",ncolls);
zone_matches[colls_w_zones] <- sapply(zone,match_one_collection_zone_to_zone_database,zone_database);
unmatched_zones <- sort(unique(zone[zone_matches[colls_w_zones]==""]));

output <- list(zone_matches,unmatched_zones);
names(output) <- c("zone_matches","unmatched_zones");
return(output);
}

# compare zones from PaleoDB collections to a zone database
match_paleodb_collections_to_possible_zones <- function(paleodb_collections,zone_database,paleodb_finds)	{
# zone_database: database of biozones with
#	zone_database$zone: eponymous taxon
#	zone_database$zone_sr: senior synonym of zone taxon (or whole zone)
#	zone_database$genus_species_combo: name omitting subgenus (=Zone if no subgenus)
#	zone_database$subgenus_species_combo: name treating subgenus as genus (=Zone if no subgenus)
#	zone_database$ma_lb: onset of zone (with 485 meaning 485 million years ago)
#	zone_database$ma_ub: end of zone (with 443 meaning 443 million years ago)
ncolls <- nrow(paleodb_collections);
paleodb_zone_info <- match_collections_zones_to_zone_database(collections=paleodb_collections,zone_database);
colls_w_zones <- (1:ncolls)[paleodb_collections$zone!=""];
zone_matches <- paleodb_zone_info$zone_matches;
#unmatched_colls <- paleodb_collections$collection_no[colls_w_zones[zone_matches[colls_w_zones]==""]];
unmatched_zone_info <- data.frame(coll_no = as.numeric(paleodb_collections$collection_no[colls_w_zones[zone_matches[colls_w_zones]==""]]),zone=as.character(paleodb_collections$zone[colls_w_zones[zone_matches[colls_w_zones]==""]]),stringsAsFactors=hell_no);

# look for zone taxa in assemblages
colls_wo_zones <- (1:ncolls)[zone_matches==""];
cwoz <- length(colls_wo_zones);
coll_id <- paleodb_collections$collection_no[colls_wo_zones];
if (is.null(paleodb_finds) || paleodb_finds=="")	{
	zone_taxa_present <- sapply(coll_id,revelio_zone_taxa_in_one_paleodb_collection,zone_database);
	zone_matches[colls_wo_zones] <- zone_taxa_present;
	}	else	{
	coll_id <- coll_id[coll_id %in% unique(paleodb_finds$collection_no)];
	zone_taxa_present <- sapply(coll_id,revelio_zone_taxa_in_paleodb_finds,paleodb_finds,zone_database);
	colls_wo_zones <- colls_wo_zones[paleodb_collections$collection_no[colls_wo_zones] %in% coll_id];
	zone_matches[colls_wo_zones] <- zone_taxa_present;
	}

output <- list(zone_matches,unmatched_zone_info);
names(output) <- c("zone_matches","unmatched_zone_info");
return(output);
}

# match a particular zone to an external zone database
match_one_collection_zone_to_zone_database <- function(zone,zone_database)	{
#print(zone);
test <- which(zone_database==zone,arr.ind=T);

# if nothing, the check to see if it is 2+ zones
if (nrow(test)==0)	{
	zones <- diffindo_zone(zone);
	zones <- zones[zones!=""];
	if (length(zones)>1)	{
		test <- matrix("",nrow=0,ncol=2);
		for (zz in 1:length(zones))
			test <- rbind(test,which(zone_database==zones[zz],arr.ind=T));
		}
	}

# if nothing, then check to see if it is a mismatch because of subgenus assignment
if (nrow(test)==0)	{
	genus <- diffindo_genus_names_from_species_names(zone);
	species <- diffindo_species_epithets(taxon_name=zone);
	if(length(strsplit(genus,split=" ")[[1]])==2 && species!="")	{
		gs <- strsplit(genus,split=" ")[[1]];
		jj <- strsplit((gs)[2],split="")[[1]];
		if (jj[1]=="(")	{
			gs[2] <- paste(jj[2:(length(jj)-1)],collapse="");
			zone_a <- paste(gs[1],species);
			zone_b <- paste(gs[2],species);
			test <- matrix("",nrow=0,ncol=2);
			test <- rbind(test,which(zone_database==zone_a,arr.ind=T));
			test <- rbind(test,which(zone_database==zone_b,arr.ind=T));
			}
		}
	}

zone_rows <- unique(test[,1]);
zone_data_reduced <- zone_database[zone_rows,];
zone_data_reduced[,(1:ncol(zone_database))[!colnames(zone_database) %in% c("zone_sr","Regional_Scale","ma_lb","ma_ub")]] <- NULL;
rownames(zone_data_reduced) <- zone_rows;
zone_data_reduced <- unique(zone_data_reduced);
zdb <- paste(unique(rownames(zone_data_reduced)),collapse=";");
return(zdb);
}

# download occurrences from a collection & look for zone taxa
revelio_zone_taxa_in_one_paleodb_collection <- function(coll_id,zone_database)	{
# zone_database: database of biozones with
#	zone_database$zone: eponymous taxon
#	zone_database$zone_sr: senior synonym of zone taxon (or whole zone)
#	zone_database$genus_species_combo: name omitting subgenus (=Zone if no subgenus)
#	zone_database$subgenus_species_combo: name treating subgenus as genus (=Zone if no subgenus)
#	zone_database$ma_lb: onset of zone (with 485 meaning 485 million years ago)
#	zone_database$ma_ub: end of zone (with 443 meaning 443 million years ago)
this_coll_finds <- accio_occurrences_from_one_paleodb_collection(coll_id);
# prepare to separate out just species occurrences
this_coll_finds$identified_rank[this_coll_finds$identified_rank=="subspecies"] <- "species"
this_coll_finds$accepted_rank[this_coll_finds$accepted_rank=="subspecies"] <- "species"
this_coll_finds <- subset(this_coll_finds,this_coll_finds$identified_rank=="species");
if (nrow(this_coll_finds) > 0)	{
	taxon_name <- this_coll_finds$identified_name;
	# replace identified name with species name if we have taxonomic data
	taxon_name[this_coll_finds$accepted_rank=="species"] <- this_coll_finds$accepted_name[this_coll_finds$accepted_rank=="species"];
	taxon_name <- sapply(taxon_name,scourgify_taxon_names);
	poss_zones <- c();
	for (tn in 1:length(taxon_name))
		poss_zones <- c(poss_zones,match_one_collection_zone_to_zone_database(zone=taxon_name[tn],zone_database));
	poss_zones <- unique(poss_zones)
	return(paste(poss_zones[poss_zones!=""],collapse=";"));
	} else	{
	return("");
	}
}

# look for zone taxa in already downloaded PaleoDB collections
revelio_zone_taxa_in_paleodb_finds <- function(coll_id,paleodb_finds,zone_database)	{
# zone_database: database of biozones with
#	zone_database$zone: eponymous taxon
#	zone_database$zone_sr: senior synonym of zone taxon (or whole zone)
#	zone_database$genus_species_combo: name omitting subgenus (=Zone if no subgenus)
#	zone_database$subgenus_species_combo: name treating subgenus as genus (=Zone if no subgenus)
#	zone_database$ma_lb: onset of zone (with 485 meaning 485 million years ago)
#	zone_database$ma_ub: end of zone (with 443 meaning 443 million years ago)
this_coll_finds <- subset(paleodb_finds,paleodb_finds$collection_no==coll_id);
# prepare to separate out just species occurrences
this_coll_finds$identified_rank[this_coll_finds$identified_rank=="subspecies"] <- "species"
this_coll_finds$accepted_rank[this_coll_finds$accepted_rank=="subspecies"] <- "species"
this_coll_finds <- subset(this_coll_finds,this_coll_finds$identified_rank=="species");
this_coll_finds <- this_coll_finds[!this_coll_finds$flags %in% c("uncertain species","uncertain genus, uncertain species"),];
if (nrow(this_coll_finds) > 0)	{
	taxon_name <- this_coll_finds$identified_name;
	# replace identified name with species name if we have taxonomic data
	taxon_name[this_coll_finds$accepted_rank=="species"] <- this_coll_finds$accepted_name[this_coll_finds$accepted_rank=="species"];
	taxon_name <- sapply(taxon_name,scourgify_taxon_names);
	poss_zones <- c();
	for (tn in 1:length(taxon_name))
		poss_zones <- c(poss_zones,match_one_collection_zone_to_zone_database(zone=taxon_name[tn],zone_database));
	poss_zones <- unique(poss_zones);
	return(paste(poss_zones[poss_zones!=""],collapse=";"));
	} else	{
	return("");
	}
}

##### BIOSTRATIGRAPHY 101 #####
# get temporal overlap between two spans (ranges, intervals, etc.)
accio_temporal_overlap <- function(lb1,ub1,lb2,ub2)	{
lb1 <- round(lb1,3);
lb2 <- round(lb2,3);
ub1 <- round(ub1,3);
ub2 <- round(ub2,3);
overlap <- data.frame(ma_lb=as.numeric(0),ma_ub=as.numeric(0));
if ((lb1 <= lb2 && lb1 > ub2) || (lb2 <= lb1 && lb2 > ub1)) {
#if (fa1 > la2 || fa2 > la1)	{
	overlap$ma_lb <- min(lb1,lb2);
	overlap$ma_ub <- max(ub1,ub2);
	}
return (overlap);
}

# find age assignments that minimize gaps and range extensions among species
optimo_paleodb_collection_and_occurrence_stratigraphy <- function(paleodb_finds,paleodb_collections,hierarchical_chronostrat,zone_database,update_search=T)	{
# rescore collections if there is any lumping of reported stages into useful stages
ncolls <- nrow(paleodb_collections);
nstages <- max(hierarchical_chronostrat$bin_last);

# delete redundant occurrences of species in localities.  (This usually reflects two co-occuring
# 	species being synonymized).
paleodb_finds <- evanesco_duplicate_occurrences_paleodb(occurrences=paleodb_finds);
	
# make sure that key collections fields are characters, not factors
paleodb_collections$early_interval <- as.character(paleodb_collections$early_interval);
paleodb_collections$late_interval <- as.character(paleodb_collections$late_interval);
paleodb_collections$interval_lb <- as.character(paleodb_collections$interval_lb);
paleodb_collections$interval_ub <- as.character(paleodb_collections$interval_ub);

# update names of un-entered species
noccr <- nrow(paleodb_finds);
to_fix <- (1:noccr)[paleodb_finds$accepted_rank %in% c("genus","subgenus")];
accepted_genus <- paleodb_finds$genus[to_fix];
identified_name <- paleodb_finds$identified_name[to_fix];
fixed_names <- c();
for (fn in 1:length(to_fix))
	fixed_names <- c(fixed_names,transmogrify_accepted_species_name(identified_name[fn],accepted_genus[fn]))
paleodb_finds$accepted_name[to_fix] <- fixed_names;

# cull out "sp.
paleodb_finds <- evanesco_indeterminate_species(paleodb_finds);
taxon_names <- sort(unique(paleodb_finds$accepted_name));

# update overall info
noccr <- nrow(paleodb_finds);
ntaxa <- length(taxon_names);
taxa <- paleodb_finds$accepted_name;

# number taxa
taxon_no <- match(taxa,taxon_names);
paleodb_finds <- cbind(paleodb_finds,taxon_no);

# add collection data to occurrence data
coll_to_find_key <- match(paleodb_finds$collection_no,paleodb_collections$collection_no);
paleodb_finds$ma_lb <- paleodb_collections$ma_lb[coll_to_find_key];
paleodb_finds$ma_ub <- paleodb_collections$ma_ub[coll_to_find_key];
paleodb_finds$interval_lb <- as.character(paleodb_collections$interval_lb[coll_to_find_key]);
paleodb_finds$interval_ub <- as.character(paleodb_collections$interval_ub[coll_to_find_key]);
#paleodb_finds <- cbind(paleodb_finds,ma_lb,ma_ub,interval_lb,interval_ub);
#paleodb_finds[105,]
#paleodb_finds[109,]
finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
finest_chronostrat <- finest_chronostrat[match(finest_chronostrat$bin_first,finest_chronostrat$bin_first),];

for (tx in 1:ntaxa)	{
	taxon_finds <- subset(paleodb_finds,paleodb_finds$accepted_name==taxon_names[tx]);
	if (sum(taxon_finds$flags %in% c("uncertain species","uncertain genus, uncertain species")) < nrow(taxon_finds))
		taxon_finds <- taxon_finds[!taxon_finds$flags %in% c("uncertain species","uncertain genus, uncertain species"),];
	if (tx==1)	{
		minimum_range_data <- tally_fuzzy_stratigraphic_ranges(ma_lb=as.numeric(taxon_finds$ma_lb),ma_ub=as.numeric(taxon_finds$ma_ub),interval_lb=as.character(taxon_finds$interval_lb),interval_ub=as.character(taxon_finds$interval_ub),hierarchical_chronostrat);
		} else {
		minimum_range_data <- rbind(minimum_range_data,tally_fuzzy_stratigraphic_ranges(ma_lb=as.numeric(taxon_finds$ma_lb),ma_ub=as.numeric(taxon_finds$ma_ub),interval_lb=as.character(taxon_finds$interval_lb),interval_ub=as.character(taxon_finds$interval_ub),hierarchical_chronostrat));
		}
	}
rownames(minimum_range_data) <- taxon_names;
#which(is.na(minimum_range_data),arr.ind = T)
# redo things for zone taxa: their latest FA cannot be after the zone starts and earliest LA cannot be before zone ends.
zone_taxa <- taxon_names[taxon_names %in% zone_database$zone];
for (zt in 1:length(zone_taxa))	{
	zn <- match(zone_taxa[zt],taxon_names);
#	if (minimum_range_data$ma_fa_mn[zn]==0)
#		print(zt)
	this_zone <- subset(zone_database,zone_database$zone==zone_taxa[zt]);
	if (minimum_range_data$ma_fa_mn[zn]==0)	{
		minimum_range_data$ma_fa_mx[zn] <- minimum_range_data$ma_fa_mn[zn] <- max(this_zone$ma_lb);
		minimum_range_data$ma_la_mx[zn] <- minimum_range_data$ma_la_mn[zn] <- min(this_zone$ma_ub);
		minimum_range_data$interval_lb[zn] <- this_zone$interval_lb[match(max(this_zone$ma_lb),this_zone$ma_lb)];
		minimum_range_data$interval_ub[zn] <- this_zone$interval_ub[match(min(this_zone$ma_ub),this_zone$ma_ub)];
		} else	{
		if (minimum_range_data$ma_fa_mx[zn]<max(this_zone$ma_lb))	{
			minimum_range_data$ma_fa_mx[zn] <- max(this_zone$ma_lb);
			minimum_range_data$interval_lb[zn] <- finest_chronostrat$interval[max(1,sum(max(this_zone$ma_lb)<=finest_chronostrat$ma_lb))];
			}
		if (minimum_range_data$ma_fa_mn[zn]<max(this_zone$ma_lb)){
			minimum_range_data$ma_fa_mn[zn] <- max(this_zone$ma_lb);
			}
		if (minimum_range_data$ma_la_mx[zn]>min(this_zone$ma_ub)){
			minimum_range_data$ma_la_mx[zn] <- min(this_zone$ma_ub);
			minimum_range_data$interval_ub[zn] <- finest_chronostrat$interval[max(1,sum(max(this_zone$ma_ub)<=finest_chronostrat$ma_ub))];
			}
		if (minimum_range_data$ma_la_mn[zn]>min(this_zone$ma_ub))	{
			minimum_range_data$ma_la_mn[zn] <- min(this_zone$ma_ub);
			}
		}
	}
#ttt <- (1:ntaxa)[minimum_range_data$ma_fa_mx<minimum_range_data$ma_la_mx];
#which(is.na(minimum_range_data),arr.ind = T)
bin_lb <- hierarchical_chronostrat$bin_first[match(paleodb_collections$interval_lb,hierarchical_chronostrat$interval)];
bin_ub <- hierarchical_chronostrat$bin_last[match(paleodb_collections$interval_ub,hierarchical_chronostrat$interval)];
# kluge: eliminate the need for this;
xx <- (1:ncolls)[bin_lb>bin_ub];
if (length(xx)>0)	{
	dummy <- bin_lb[xx];
	bin_lb[xx] <- bin_ub[xx];
	bin_ub[xx] <- dummy;
	paleodb_collections$interval_lb[xx] <- finest_chronostrat$interval[bin_lb[xx]];
	paleodb_collections$interval_ub[xx] <- finest_chronostrat$interval[bin_ub[xx]];
	}
#paleodb_collections$early_interval[xx]
#paleodb_collections$late_interval[xx]
#paleodb_collections$interval_lb <- finest_chronostrat$interval[bin_lb];
#paleodb_collections$interval_ub <- finest_chronostrat$interval[bin_ub];

#unfixed_collections <- paleodb_collections$collection_no[as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub)];
bin_fuzz <- bin_ub-bin_lb;
unfixed_collections <- paleodb_collections$collection_no[bin_fuzz>0];

#paste(unfixed_collections,collapse = ",")
orig_unfixed <- unfixed <- length(unfixed_collections);
improved <- fixed <- ncolls - unfixed;
min_gap_to_set <- 0;
min_taxa_to_set <- 3;
attempt <- 1;
reboots <- 0;
# paste(unfixed_collections,collapse=",");
fuzz <- 0:max(bin_ub-bin_lb);
#hist(bin_ub-bin_lb)
fcolls <- hist(bin_fuzz,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
progress <- cbind(fuzz,fcolls);
if (update_search)
	print(progress);
while (improved > 0)	{
	newly_fixed <- new_and_improved <- 0;
	uc <- 0;	# uc <- match(problems[12],unfixed_collections)
 	while (uc < unfixed)	{
		uc <- uc + 1;
		coll_no <- match(unfixed_collections[uc],paleodb_collections$collection_no);
		index_species <- paleodb_finds$accepted_name[paleodb_finds$collection_no==unfixed_collections[uc]];
		index_ranges_all <- minimum_range_data[match(index_species,rownames(minimum_range_data)),];
		index_ranges <- subset(index_ranges_all,index_ranges_all$ma_fa_mn!=0);
		if (nrow(index_ranges)>0)	{
			ranges <- cbind(hierarchical_chronostrat$bin_first[match(index_ranges$interval_lb,hierarchical_chronostrat$interval)],hierarchical_chronostrat$bin_last[match(index_ranges$interval_ub,hierarchical_chronostrat$interval)]);
			# something is resetting bin_lb & bin_ub to NA; might be ma_lb &/or ma_ub doing it
			poss_range <- bin_lb[coll_no]:bin_ub[coll_no];
			
			bin_gaps <- array(0,dim=c(length(poss_range)));
			for (rr in 1:nrow(ranges))
				bin_gaps <- bin_gaps+procrustean_binning_one_taxon(range=ranges[rr,],poss_range=c(min(poss_range),max(poss_range)),debug=F);

			# case where some of the previously assigned range induces no gaps in known species
			if ((min(bin_gaps)<=min_gap_to_set && length(poss_range)>sum(bin_gaps<=min_gap_to_set)) && nrow(index_ranges)>=min_taxa_to_set)	{
				bin_lb[coll_no] <- min(poss_range[bin_gaps==min(bin_gaps)]);
				bin_ub[coll_no] <- max(poss_range[bin_gaps==min(bin_gaps)]);
				new_and_improved <- new_and_improved+1;
				paleodb_collections$interval_lb[coll_no] <- finest_chronostrat$interval[match(bin_lb[coll_no],finest_chronostrat$bin_first)];
				paleodb_collections$interval_ub[coll_no] <- finest_chronostrat$interval[match(bin_ub[coll_no],finest_chronostrat$bin_last)];
				paleodb_collections$ma_lb[coll_no] <- min(paleodb_collections$ma_lb[coll_no],finest_chronostrat$ma_lb[match(bin_lb[coll_no],finest_chronostrat$bin_first)]);
				paleodb_collections$ma_ub[coll_no] <- max(paleodb_collections$ma_ub[coll_no],finest_chronostrat$ma_ub[match(bin_ub[coll_no],finest_chronostrat$bin_last)]);
				# case where we've narrowed it down to one bin
				if (sum(bin_gaps==min_gap_to_set)==1)	{
					# the earliest possible first occurrences must be at least as old as the oldest possible first occurrences
					#	for those taxa first known from this interval
					sub_index_ranges <- subset(index_ranges,index_ranges$interval_lb==paleodb_collections$interval_lb[coll_no]);
					if (nrow(sub_index_ranges)>0)	{
						sub_index_ranges$ma_fa_mx[sub_index_ranges$ma_fa_mx<max(sub_index_ranges$ma_fa_mn)] <- max(sub_index_ranges$ma_fa_mn);
						index_ranges[match(rownames(sub_index_ranges),rownames(index_ranges)),] <- sub_index_ranges
						}
					
					# the latest possible last occurrences must be at least as young as the youngest possible last occurrences
					#	for those taxa last known from this interval
					sub_index_ranges <- subset(index_ranges,index_ranges$interval_ub==paleodb_collections$interval_ub[coll_no]);
					if (nrow(sub_index_ranges)>0)	{
						sub_index_ranges$ma_la_mn[sub_index_ranges$ma_la_mn>min(sub_index_ranges$ma_la_mx)] <- min(sub_index_ranges$ma_la_mx);
						index_ranges[match(rownames(sub_index_ranges),rownames(index_ranges)),] <- sub_index_ranges
						}
					indexed_species <- rownames(index_ranges);
					i_s <- match(indexed_species,rownames(minimum_range_data));
					minimum_range_data[i_s,] <- index_ranges;
					
					if (nrow(index_ranges_all) > nrow(index_ranges))	{
					# update species that do not have constrained finds yet
						unindexed_species <- rownames(index_ranges_all)[index_ranges_all$ma_fa_mn==0];
						u_s <- match(unindexed_species,rownames(minimum_range_data));
						sub_index_ranges <- subset(index_ranges,index_ranges$interval_lb==paleodb_collections$interval_lb[coll_no]);
						if (nrow(sub_index_ranges)>0)	{
							minimum_range_data$ma_fa_mx[u_s] <- min(sub_index_ranges$ma_fa_mx);
							minimum_range_data$ma_fa_mn[u_s] <- max(sub_index_ranges$ma_fa_mn);
							} else	{
							minimum_range_data$ma_fa_mx[u_s] <- paleodb_collections$ma_lb[coll_no];
							minimum_range_data$ma_fa_mn[u_s] <- paleodb_collections$ma_ub[coll_no];
							}
						sub_index_ranges <- subset(index_ranges,index_ranges$interval_ub==paleodb_collections$interval_ub[coll_no]);
						if (nrow(sub_index_ranges)>0)	{
							minimum_range_data$ma_la_mx[u_s] <- min(sub_index_ranges$ma_la_mx);
							minimum_range_data$ma_la_mn[u_s] <- max(sub_index_ranges$ma_la_mn);
							} else	{
							minimum_range_data$ma_la_mx[u_s] <- paleodb_collections$ma_lb[coll_no];
							minimum_range_data$ma_la_mn[u_s] <- paleodb_collections$ma_ub[coll_no];
							}
						minimum_range_data$interval_lb[u_s] <- paleodb_collections$interval_lb[coll_no];
						minimum_range_data$interval_ub[u_s] <- paleodb_collections$interval_ub[coll_no];
						}
					newly_fixed <- newly_fixed+1;
					}
				}
 			}
 		if (sum(bin_lb==Inf)>0 || sum(is.na(bin_lb))>0)
 			print(uc);
#		print(cbind(uc,minimum_range_data[2309,]));
#		print(c(paleodb_collections$ma_lb[195],paleodb_collections$ma_ub[195]));
		}
#	(1:ncolls)[is.na(as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub))];
	unfixed_collections <- paleodb_collections$collection_no[as.character(paleodb_collections$interval_lb)!=as.character(paleodb_collections$interval_ub)];
	unfixed <- length(unfixed_collections);
	improved <- new_and_improved;
	fixed <- newly_fixed;
	unfixed <- length(unfixed_collections);
	if (improved==0 && min_taxa_to_set>0)	{
		reboots <- reboots+1;
		if ((reboots %% 2)==0)	{
			min_taxa_to_set <- min_taxa_to_set-1;
			} else	{
			min_gap_to_set <- min_gap_to_set+1;
			}
		if (min_taxa_to_set>0)
			improved <-max(1,improved);
		}
	attempt <- attempt+1;
	fcolls <- hist(bin_ub-bin_lb,breaks=c(min(fuzz)-1,fuzz),plot=F)$counts;
#	fcolls <- hist(bin_ub-bin_lb,breaks=c(fuzz,max(fuzz)+1),plot=F)$counts;
	progress <- cbind(progress,fcolls);
	if (update_search)
		print(progress);
	}
return(paleodb_collections);
}

# find collections that might be in 2+ chronostratigraphic bins
tally_fuzzy_stratigraphic_ranges <- function(ma_lb,ma_ub,interval_lb,interval_ub,hierarchical_chronostrat)	{
fa_latest <- min(unique(hierarchical_chronostrat$bin_last[match(interval_ub,hierarchical_chronostrat$interval)]));
la_earliest <- max(unique(hierarchical_chronostrat$bin_first[match(interval_lb,hierarchical_chronostrat$interval)]));

#ttl_finds <- length(coll_no);
bins_early <- hierarchical_chronostrat$bin_first[match(interval_lb,hierarchical_chronostrat$interval)];
bins_late <- hierarchical_chronostrat$bin_last[match(interval_ub,hierarchical_chronostrat$interval)];
bins_early[bins_early<fa_latest] <- fa_latest;
bins_late[bins_late>la_earliest] <- la_earliest;

definite_bins <- unique(bins_early[bins_early==bins_late]);
output <- data.frame(ma_fa_mx=0,ma_fa_mn=0,ma_la_mx=0,ma_la_mn=0,interval_lb=as.character(""),interval_ub=as.character(""),stringsAsFactors=hell_no);
if (length(definite_bins)>0)	{
	finest_chronostrat <- hierarchical_chronostrat[hierarchical_chronostrat$bin_first==hierarchical_chronostrat$bin_last,];
	# extremes are easy
	output$ma_fa_mx <- as.numeric(max(ma_lb));						# earliest possible FA
	output$ma_la_mn <- as.numeric(min(ma_ub));						# latest possible LA
	# latest FA must precede or coincide with latest possible LA & vice versa
	output$ma_fa_mn <- as.numeric(max(ma_ub));						# latest possible FA
	output$ma_la_mx <- as.numeric(min(ma_lb));						# earliest possible LA
	output$interval_lb <- as.character(finest_chronostrat$interval[match(min(definite_bins),finest_chronostrat$bin_first)]);
	output$interval_ub <- as.character(finest_chronostrat$interval[match(max(definite_bins),finest_chronostrat$bin_first)]);
	}
return(output);
}

# score how many sum of range extensions required to make
procrustean_binning_one_taxon <- function(range,poss_range,debug=FALSE)	{
# range: taxon range in bin numbers (with 1 older than 2)
# poss_range: lower and upper bounds on collection age in bin numbers
if (debug)	print(range);
obs_bins <- range[1]:range[2];
poss_bins <- poss_range[1]:poss_range[2];
#ok_bins <- poss_bins[poss_bins %in% obs_bins];
#if (length(ok_bins)==0)
#	ok_bins <- poss_bins;
bin_score <- c();
for (pb in 1:length(poss_bins))	{
	bin_score <- c(bin_score,min(abs(poss_bins[pb]-obs_bins)));
	}
return(bin_score);
}

##### CAL3 #####
### Routines for Phylogeny Prior Probabilities
prob_sampling_clade_bapst <- function(p,q,r)	{
# modified from Bapst (2013) equation 2
# p: origination (lambda in Bapst 2013)
# q: extinction  (mu in Bapst 2013)
# r: sampling (psi in Bapst 2013)
s <- 1		# richness of clade (K in Bapst's paper)
pmiss <- mxp <- 0
num <- den <- 1
while ((mxp/100000)<(num/den))	{
	num <- (p^(s-1))*(q^s)*choose((2*s)-2,s-1)
	den <- s*((p+q+r)^((2*s)-1))
	pmiss <- pmiss+(num/den)
	if (mxp<(num/den))	mxp <- (num/den)
#	print(c(s,num,den,num/den,pmiss))
	s <- s+1
	}
return(1-pmiss)
}

##### SAMPLING DISTRIBUTION ROUTINES ####
# second-order (modified) AIC
modified_AIC <- function(lnL,k,n)	{
# L: log-likelihood; # k: parameters; # n: data points
#if (n==k)	n <- n+2
if (is.na(n / (n - k - 1)) || (n / (n - k - 1)<0) || (n / (n - k - 1))==Inf)	{
	aic_c <- AIC(lnL,k);
	}	else	aic_c <- (-2*lnL) + (2*k)*(n / (n - k - 1));
return(aic_c);
}

# routine to estimate most likely uniform rate of sampling/occupancy for hS taxa given finds & possible finds
# convert abundances to "Fisher plot" giving # taxa with 1…N Finds
fisher_plot <- function(finds)	{
unique_finds <- sort(unique(finds),decreasing=FALSE)
observed <- vector(length=max(unique_finds))
for (s in 1:length(unique_finds))	observed[unique_finds[s]] <- length(finds[finds==unique_finds[s]])
return(observed)
}

# Get Chao2 Richness estimate from vector giving # taxa with 1…N finds
Chao2_Fisher <- function(observed)	{
ntaxa <- sum(observed)
S <- round(ntaxa + ((observed[1]*observed[1])/(2*(observed[2]+1))) - ((observed[1]*observed[2])/((2*(observed[2]+1)*(observed[2]+1)))))
return(S)
}

# Get Jackknife 2 Richness estimate from vector giving # taxa with 1…N finds
jack2_Fisher <- function(observed)	{
ntaxa <- sum(observed)
ss <- 0
for (i in 1:length(observed))	ss <- ss+(i*observed[i])
S <- round(ntaxa + (observed[1]*(((2*ss)-3)/ss)) - observed[2]*(((ss-2)*(ss-2))/(ss*(ss-1))))
return(S)
}

# get expected occurrences given hypothesized occupancy distribution and numbers of collections
#	this generates an "expected" analog to fisher_plot above
expected_occurrences <- function(rocd, ncoll, S)	{
# find expectations of this gamma at this sample size
# exp will give the expected number of taxa found 1…ncoll times
b <- 1:ncoll;
expected <- dbinom(b,ncoll,rocd[1]);
for (j in 2:S)
	expected <- expected+dbinom(b,ncoll,rocd[j]);
min <- 1.0e-323;
for (i in 1:ncoll)	if (expected[i]==0)	expected[i] <- min;
return(expected);
}

## Uniform Distrributions
# get multinomial lognormal including P[observed total | expected total] as substitute for P[0]
distribution_loglikelihood_mul <- function(observed,expected,oS,hS)	{
#print(c(oS,hS))
mxfind <- length(observed)							# maximum finds observed
prop_expected <- expected[1:mxfind]/sum(expected)		# convert expected species to proportions
prop_expected[prop_expected==0] <- MINNO
lnlo <- observed*log(prop_expected)	# exact probability of observing # taxa with 1, 2, 3, etc. finds
lnlo[is.na(lnlo)] <- 0
#for (i in 1:mxfind)	if (is.na(lnlo[i]))	lnlo[i] <- 0
sobs <- sum(lnlo)
# log probability of observing X taxa from hypothesized hS taxa given expected # taxa with 1, 2, 3, etc. finds
eS <- sum(expected)									# get expected sampled species
if (eS==hS)	eS <- 0.99999999999*hS				# this is a kluge to get around rounding error of eS->hS
if (eS<hS)	{
	lnls <- lfactorial(hS)-(lfactorial(hS-oS)+lfactorial(oS))+(oS*log(eS/hS))+((hS-oS)*log((hS-eS)/hS))
	}	else if (oS<hS)	{
	lnls <- oS*log(MINNO)
	}	else	{
	lnls <- 0	# if we expect to observe all of the species and oS==hS, then P=1.0
	}
return(sobs+lnls)
}

# routine to get loglikelihood of a uniform rate of sampling among hS taxa given finds and possible finds
loglikelihood_uniform_occupancy_rates <- function(sc, hS, observed, oS, ncoll)	{
rocd <- rep(sc,hS)
#print(sc)
if (rocd[1]<=1)	{
	expected <- expected_occurrences(rocd,ncoll,hS)
		# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
#	print(c(round(sc,20),round(lnl,3)))
	if (is.na(lnl))	lnl <- oS*log(MINNO)
	}	else {
	lnl <- oS*log(MINNO)
	}
return(lnl)
}

# routine to estimate most likely uniform rate of sampling/occupancy for hS taxa given finds & possible finds
optimize_uniform_occupancy_for_hS <- function(hS,observed,oS,ncoll)	{
#print(hS)	# for debugging
cl <- list(fnscale=-1)	# no idea what this means.....
sc <- (sum((1:length(observed))*observed)/hS)/ncoll
minsc <- sc/100
maxsc <- length(observed)/ncoll
w <- optim(sc,fn=loglikelihood_uniform_occupancy_rates,method="L-BFGS-B",hS=hS,oS=oS,observed=observed,ncoll=ncoll,control=cl,lower=minsc,upper=maxsc);
bH <- c(w$par,hS,w$value)
names(bH) <- c("scale","richness","loglikelihood")
return(bH)
}

# routine to estimate most likely uniform rate of sampling/occupancy given finds & possible finds
optimize_uniform_occupancy <- function(finds,ncoll)	{
observed <- fisher_plot(finds)
oS <- length(finds[finds>0])	# observed taxa
minS <- stS <- length(finds)
pa <- 1
pz <- span <- 4		# span of richnesses to consider
enS <- stS+((span-1)*minS)	# highest richness to consider
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
last_hS <- c(0,1,2,3)
while (incr>0)	{
	hS <- round(seq(stS,enS,by=incr),0)
	if (pz>length(hS))	pz <- length(hS)
	if (sum(hS %in% last_hS) == length(hS))	break
	results <- sapply(hS[pa:pz],optimize_uniform_occupancy_for_hS,observed=observed,oS=oS,ncoll=ncoll)
#	results <- sapply(hS[pa:pz],optimize_uniform_occupancy_for_hS,observed=observed,oS=oS,ncoll=ncoll)
	if (pa==2)
		results <-cbind(old_result_1,results)
	if (pz<length(hS))
		results <-cbind(results,old_result_2)
	mlnl <- max(results[3,])
	mlSc <- match(mlnl,results[3,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- length(results[3,])
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				pa <- 2		# make note of best hypothesis
				old_result_1 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				stS <- hS[spn]
				enS <- stS+((spn-1)*minS)
				}	else	{
				pa <- 2				# make note of 2nd best hypothesis
				pz <- (spn-1)		# make note of best hypothesis
				old_result_1 <- results[,spn-1]	# set aside stats so as to avoid recalculating
				old_result_2 <- results[,spn]	# set aside stats so as to avoid recalculating
				enS <- hS[spn]
				stS <- hS[spn-1]
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- (enS-stS)
					pz <- span-1
					incr <- 1
					}
#				stS <- hS[spn]-(span*incr)
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]
			old_result_1 <- results[,1]
			pa <- 2
			old_result_2 <- results[,2]
			pz <- length(hS)-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- (enS-stS)/(span-1)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_uniform_occupancy_for_hS,observed=observed,oS=oS,ncoll=ncoll)
			if (results2[3,1]>mlnl && results2[3,1]>results2[3,2])	{
			# start just above 2nd best richness and go up to best
				stS <- hS[mlSc-1]
				enS <- hS2[1]
				old_result_1 <- results[,mlSc-1]
				old_result_2 <- results2[,1]
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else if (results2[3,2]>mlnl && results2[3,2]>results2[3,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS2[2]
				enS <- hS[mlSc+1]
				old_result_1 <- results2[,2]
				old_result_2 <- results[,mlSc+1]
				if (incr>span)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
	last_hS <- hS
	}
bH[3]<-round(bH[3],2)
uniform_AICc <- round(modified_AIC(bH[3],2,sum(finds)),2)
names(uniform_AICc) <- "AICc"
bH <- c(bH,uniform_AICc)
return(bH)
}

## Exponential Distributions
# upper and lower bounds for plausible exponentials
accio_param_bounds_exponential <- function(finds,observed,ncoll)	{
# step 1: get numbers of species with 1, 2, 3, etc. finds
prop_drops <- unique(accio_prop_drops(finds))
oS <- length(finds)
decay_min <- 1+((min(subset(prop_drops,prop_drops>1))-1)/5)
decay_max <- min(max(prop_drops),1/exp(log(MINEXPN)/oS))
# get the exponential distributions from the two extreme slopes
dmin <- exponential_distribution(decay_min)
dmax <- exponential_distribution(decay_max)
# set scaling so that the most common taxon will be within one half to twice the prop
#	of the observed
scale_min <- (finds[1]/ncoll)/dmax[1]
scale_max <- (finds[1]/ncoll)/dmin[1]
return(c(scale_min,scale_max,decay_min,decay_max))
}

accio_richness_with_prop_greater_than_x_given_exponential <- function(decay,x)	{
d <- min(decay,1/decay)
return(round(1+((log(x)-log(1-d))/log(d))))
}

# exponential distribution with decay rate 1/ev; note that this is truncated at some minimum number
# generate geometric (= exponential) distribution
exponential_distribution <- function(decay)	{
#print(decay);
d <- min(decay,1/decay)
S <- accio_richness_with_prop_greater_than_x_given_exponential(decay,MINEXPN)
ranks <- (1:S)
prop <- (1-d)*(d^(ranks-1))
return(prop)
}

# find the decays from one rank to the next
accio_prop_drops <- function(finds)	{
oS <- length(finds)
s1 <- (1:(oS-1))
s2 <- (2:oS)
return(finds[s1]/finds[s2])
}

# rough initial estimate of an exponential distribution to seed searches
accio_initial_occupancy_exponential <- function(finds,observed,ncoll)	{
p0 <- vector(length=2)	# pO[1]=rate adjuster; p0[2]=decay
names(p0) <- c("scale", "decay")
# step 1: estimate average shifts
prop_drops <- accio_prop_drops(finds[finds>0]);
p0[2] <- mean(prop_drops)
d <- exponential_distribution(p0[2])
p0[1] <- (finds[1]/ncoll)/d[1]
names(p0) <- c("scale","decay")
return(p0)
}

# get loglikelihood of exponential distribution of sampling rates given finds and collections
loglikelihood_exponential_occupancy_rates <- function(p0, observed, oS, ncoll)	{
r <- p0[1]
decay <- p0[2]
#print(p0)
# p0[1]: r	# p0[2]: ev	# p0[3]: S
rocd <- scaled_exponential_distribution(r,decay)		# basic exponential distribution
if (length(rocd)<oS)	{
	hS <- Chao2_Fisher(observed)
	dummy <- vector(length=(hS-length(rocd)))
	dummy[1] <- rocd[length(rocd)]/decay
	if (length(dummy)>1)		for (i in 2:length(dummy))	dummy[i] <- dummy[i-1]/decay
	rocd <- c(rocd,dummy)
#	rocd <- scaled_exponential_distribution_min_rich(r,decay,hS)
	}	else	{
	minp <- ((10^-7)/ncoll)
	hS <- min(Chao2_Fisher(observed),length(subset(rocd,rocd>minp)))
	}
#rocd <- rocd*max(finds)/ncoll
if (rocd[1]<=1 && hS>=oS)	{
	expected <- expected_occurrences(rocd,ncoll,hS)[1:length(observed)]
		# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
	}	else {
	lnl <- oS*log(MINNO)
	}
#print(c(r,decay,lnl))
return(round(lnl,2))
}

# rough initial estimate of an exponential distribution to seed searches
accio_initial_occupancy_exponential <- function(finds,observed,ncoll)	{
p0 <- vector(length=2)	# pO[1]=rate adjuster; p0[2]=decay
names(p0) <- c("scale", "decay")
# step 1: estimate average shifts
prop_drops <- accio_prop_drops(finds[finds>0]);
p0[2] <- mean(prop_drops)
d <- exponential_distribution(p0[2])
p0[1] <- (finds[1]/ncoll)/d[1]
names(p0) <- c("scale","decay")
return(p0)
}

# exponential distribution with decay rate that is rescaled to some "average" rate; note that this is truncated at some minimum number
scaled_exponential_distribution <- function(sc,decay)	{
return(sc*exponential_distribution(decay))
}

# optimize exponential distribution
optimize_exponential_occupancy <- function(finds,ncoll)	{
observed <- fisher_plot(finds);
oS <- length(finds);	# observed taxa
finds <- finds[finds>0];
#prop_finds <- finds/ncoll
p0 <- accio_initial_occupancy_exponential(finds,observed,ncoll);	# pO[1]=; p0[3]=richness
bnds <- accio_param_bounds_exponential(finds,observed,ncoll);
cl <- list(fnscale=-1)
w <- optim(p0,fn=loglikelihood_exponential_occupancy_rates,method="L-BFGS-B",oS=oS,observed=observed,ncoll=ncoll,lower=c(bnds[1],bnds[3]),upper=c(bnds[2],bnds[4]),control=cl)

d <- 1/w$par[2]
hS <- 1+round((log((10^-7)/(1-d)))/log(d),0)
bH <- c(w$par,hS,w$value,round(modified_AIC(w$value,2,sum(finds)),2))
names(bH) <- c("scale","decay","richness","loglikelihood","AICc")
return(bH)
}

## Beta Distribution
# get beta distribution
beta_distribution <- function(shape1, shape2, S)	{
# NOTE:I get some wonky assed results using dbeta & qbeta
ranks <- (1:S)/(S+1)
numer <- (ranks^(shape1-1))*((1-ranks)^(shape2-1))
#numer[numer==0] <- ZERO
# use log gammas, as gamma does not like numbers much past 100
Bab <- exp((lgamma(shape1)+lgamma(shape2))-lgamma(shape1+shape2))
cdf <- numer/Bab
pdf <- cdf/sum(cdf)
return(pdf)
}

# get rough estimate
accio_initial_occupancy_beta <- function(obs_mean, obs_var)	{
shape1 <- 0.1
shape2 <- (shape1/obs_mean)-shape1
combo_var <- (shape1*shape2)/(((shape1+shape2)^2)*(shape1+shape2+1))

if (combo_var>obs_var)	{
	while (combo_var>obs_var)	{
		shape1 <- 1.01*shape1
		shape2 <- (shape1/obs_mean)-shape1
		combo_var <- (shape1*shape2)/(((shape1+shape2)^2)*(shape1+shape2+1))
		}
	}	else if (combo_var<obs_var)	{
	while (combo_var<obs_var)	{
		shape1 <- shape1/1.01
		shape2 <- (shape1/obs_mean)-shape1
		combo_var <- (shape1*shape2)/(((shape1+shape2)^2)*(shape1+shape2+1))
		}
	}
return(c(shape1,shape2))
}

# calculate loglikelihood of beta for occupancy
loglikelihood_beta_occupancy_rates <- function(p0, hS, observed, oS, ncoll)	{
# p0[1]: shape1	# p0[2]: beta	
shape1 <- p0[1]
shape2 <- p0[2]
#print(c(shape1,shape2,hS,oS))
mxobs <- length(observed)
if (shape2>=shape1)	{
	rocd <- beta_distribution(shape1, shape2, hS)
	#print(c(r,mag))
	if (rocd[1]<=1 && length(rocd[rocd>0])==length(rocd))	{
	# expected number of taxa with 1, 2, 3, etc. finds
	#	expected_prop <- expected_prop(rocd,ncoll,p0[3])
		expected <- expected_occurrences(rocd,ncoll,hS)
		# log likelihood
		lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
		}	else {
		lnl <- 10000*-oS
		}
	}	else lnl <- 10000*-oS
return(round(lnl,3))
}

# find best beta occupancy distribution
optimize_beta_occupancy_given_hS <- function(hS,oS,finds,observed,ncoll)	{
#print(hS)
obs_mean <- mean(c(finds,rep(0,hS-oS)))/ncoll
obs_var <- var(c(finds/ncoll,rep(0,hS-oS)))
p0 <- accio_initial_occupancy_beta(obs_mean,obs_var)
names(p0) <- c("beta_alpha","beta_beta")
cl <- list(fnscale=-1)
w <- optim(p0,fn=loglikelihood_beta_occupancy_rates,method="L-BFGS-B",oS=oS,hS=hS,observed=observed,ncoll=ncoll,lower=p0/100, upper=100*p0,control=cl)
bH <- c(w$par,hS,w$value);
names(bH) <- c("alpha", "beta", "richness","loglikelihood");
return(bH);
}

# find best beta occupancy distribution
optimize_beta_occupancy <- function(finds,ncoll)	{
observed <- fisher_plot(finds)
oS <- length(finds[finds>0])		# observed taxa
minS <- length(finds)				# observed + gap taxa
#bS <- min(Chao2_Fisher(observed),jack2_Fisher(observed))	# lowest richness estimator we will consider
pa <- 1
pz <- span <- 4
stS <- minS
enS <- stS+((span-1)*minS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
last_hS <- c(0,1,2,3)
while (incr>0)	{
	hS <- round(seq(stS,enS,by=incr),0)
	if (pz>length(hS))	pz <- length(hS)
	if (sum(hS %in% last_hS) == length(hS))	break
	results <- sapply(hS[pa:pz],optimize_beta_occupancy_given_hS,oS=oS,finds=finds,observed=observed,ncoll=ncoll)
	if (pa==2)
		results <-cbind(old_result_1,results)
	if (pz<length(hS))
		results <-cbind(results,old_result_2)
	mlnl <- max(results[4,])
	mlSc <- match(mlnl,results[4,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- length(results[1,])
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				pa <- 2		# make note of best hypothesis
				old_result_1 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				stS <- hS[spn]
				enS <- stS+((spn-1)*minS)
				}	else	{
				stS <- hS[mlSc-1]
				enS <- hS[mlSc]
				old_result_1 <- results[,mlSc-1]	# set aside stats so as to avoid recalculating
				old_result_2 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
#				stS <- hS[spn]-(span*incr)
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]
			old_result_1 <- results[,1]
			pa <- 2
			old_result_2 <- results[,2]
			pz <- length(hS)-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- (enS-stS)/(span-1)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_beta_occupancy_given_hS,oS=oS,finds=finds,observed=observed,ncoll=ncoll)
			if (results2[4,1]>mlnl && results2[4,1]>results2[4,2])	{
			# start just above 2nd best richness and go up to best
				stS <- hS[mlSc-1]
				enS <- hS2[1]
				old_result_1 <- results[,mlSc-1]
				old_result_2 <- results2[,1]
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else if (results2[4,2]>mlnl && results2[4,2]>results2[4,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS2[2]
				enS <- hS[mlSc+1]
				old_result_1 <- results2[,2]
				old_result_2 <- results[,mlSc+1]
				if (incr>span)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			}
		# end case where we have a better richness in middle somewhere.  
		}	else	{
		incr <- 0
		}
#	print(hS)		# for debugging
#	if (incr>0)	print(round(seq(stS,enS,by=incr),0))	# for debugging
	last_hS <- hS
	}
beta_AICc <- round(modified_AIC(mlnl,3,sum(finds)),2)
bH <- c(bH[1],bH[2],bH[3],round(bH[4],2),beta_AICc)
names(bH)[5] <- "AICc"
return(bH)
}

## Lognormal distributions
# generate lognormal distribuiton for S entities with log-variance = ev
lognormal_distribution <- function(mag, S)	{
fi <- seq(1/(S+1),S/(S+1),by=1/(S+1))
prop <- mag^(qnorm(fi,0,1))/sum(mag^(qnorm(fi,0,1)))
return(prop)
}

# generate lognormal distribuiton for S entities with log-variance = mag and rescaled by r
scaled_lognormal_distribution <- function(sc, mag, S)	{
return(sort(sc*lognormal_distribution(mag, S),decreasing=TRUE))
}

# get lower and upper bounds for lognormal magnitude and rescaling
accio_min_max_lognormal_occupancy_params <- function(finds,observed,oS,ncoll)	{
mmm <- finds[1]/finds[oS]
# use a maximum richness estimate to get maximum magnitude
#	this works because the Nth taxon is a higher overall rank given high richness
#	we therefore go from obs. min -> obs. max over the smallest shift in rel. ranks
smax <- max(Chao2_Fisher(observed),jack2_Fisher(observed))
smb <- 1-(oS/(smax+1))
s1b <- 1-(1/(smax+1))
qmb <- qnorm(smb,0,1)
q1b <- qnorm(s1b,0,1)

mag_max <- 2*exp(log(mmm)/(q1b-qmb))

### use observed richness to get lowest estimate.
smin <- oS
sma <- 1-(oS/(smin+1))
s1a <- 1-(1/(smin+1))
qma <- qnorm(sma,0,1)
q1a <- qnorm(s1a,0,1)
mag_min <- exp(log(mmm)/(q1a-qma))/2

# use fewest taxa and biggest magnitude shift to get the most common possible #1
ru <- sort(lognormal_distribution(mag_max,smin),decreasing=TRUE)
# use most taxa and smallest magnitude shift to get the least common possible #1
rl <- sort(lognormal_distribution(mag_min,smax),decreasing=TRUE)

sc_raw <- finds[1]/ncoll
sc_min <- sc_raw/ru[1]
sc_max <- sc_raw/rl[1]
#ru2 <- sort(lognormal_distribution(mag_max,smax),decreasing=TRUE)
#print(c(rl1[1],rl2[1],ru1[1],ru2[1]))
boundary_params <- c(mag_min,mag_max,sc_min,sc_max)
names(boundary_params) <- c("magn_min","magn_max","scale_min","scale_max")
return(boundary_params)
}

# get rough estimate of lognormal given variance in log occupancy 
accio_initial_occupancy_lognormal <- function(finds,observed,oS,minS,ncoll)	{
# editted 2019-08-15
p0 <- vector(length=3)	# pO[1]=; p0[3]=richness
names(p0) <- c("scale", "magnitude", "richness")
	# step 1: get rough estimate of lognormal based on variance above the median
p0[3] <- min(Chao2_Fisher(observed),jack2_Fisher(observed))	# lowest richness estimator we will consider
if (p0[3]<=minS)
	p0[3] <- minS+(p0[3]-oS);

if (p0[3]<=oS)	{
	mid <- round(p0[3]/2);
	for (i in 1:min(p0[3],10))	{			# do not let number exceed richness!
		s <- p0[3]-(i-1)
		x1 <- qnorm(s/(s+1),0,1)			# get q-value associated with dominant taxon
		p0[2] <- p0[2]+exp(log(finds[i]/finds[mid])/x1)/10
		}	# get rough estimate of magnitude parameter based on Top 10 taxa
	d <- lognormal_distribution(p0[2],p0[3])
	p0[1] <- (finds[mid]/ncoll)/d[mid]	# set median rate to equal rate of extrapolated median
	}	else	{
	dc <- length(subset(observed,observed>0))
	k <- 1
	uniqfnds <- sort(unique(finds),decreasing=TRUE)
	uo <- observed[observed>0]
	uniqobs <- uo[dc:1]
	midpts <- vector(length=dc)
	midpts[1] <- (1+uniqobs[1])/2
	for (k in 2:dc)	{
		midpts[k] <- midpts[k-1]+((1+uniqobs[k])/2)
		}
	quants <- qnorm((p0[3]-(midpts-1))/(p0[3]+1))	# mag^quant[i] is how many times more common a species is than the median
	### each quant[i] gives the quantile of the median ranked taxon with uniqfnds[i] specimens
	for (s1 in 1:(dc-1))	{
		f1 <- uniqfnds[s1]
		for (s2 in (s1+1):dc)	{
			f2 <- uniqfnds[s2]
			if (s1==1 && s2==2)	{
				avem <- (f1/f2)^(1/(quants[s1]-quants[s2]))
				}	else	{
				avem <- c(avem,(f1/f2)^(1/(quants[s1]-quants[s2])))
				}
			}
		}
	# use the minimum of the median, arithmetic mean and geometric mean.
	p0[2] <- min(median(avem),mean(avem),(exp(sum(log(avem))))^(1/length(avem)))
	d <- sort(lognormal_distribution(p0[2],p0[3]),decreasing=TRUE)
	p0[1] <- median((finds[1:oS]/ncoll)/d[1:oS]);
	if ((p0[1]*d[1])>1)	p0[1] <- 0.99/d[1]
	}
return(p0)
}

# calculate loglikelihood of lognormal for occupancy
loglikelihood_lognormal_occupancy_rates <- function(p0, hS, observed, oS, ncoll)	{
# p0[1]: r	# p0[2]: ev	# p0[3]: S		print(p0)
r <- p0[1]
mag <- p0[2]
rocd <- sort(scaled_lognormal_distribution(r,mag,hS),decreasing=TRUE)		# basic lognormal distribution
#print(c(r,mag))
if (rocd[1]<=1 && length(rocd[rocd>0])==length(rocd))	{
	# expected number of taxa with 1, 2, 3, etc. finds
	#	expected_prop <- expected_prop(rocd,ncoll,p0[3])
	expected <- expected_occurrences(rocd,ncoll,hS)
	# log likelihood
	lnl <- distribution_loglikelihood_mul(observed,expected,oS,hS)
	}	else {
	lnl <- 10000*-oS
	}
return(lnl)
}

# find best lognormal occupancy distribution for particular richness
optimize_lognormal_occupancy_given_hS <- function(hS,oS,observed,ncoll,p0,sc,mm,debug=0)	{
if (debug==1)	print(hS)
#print(hS)
cl <- list(fnscale=-1)
w <- optim(c(br=p0[1],mag=p0[2]),fn=loglikelihood_lognormal_occupancy_rates,method="L-BFGS-B",oS=oS,hS=hS,observed=observed,ncoll=ncoll,lower=c(min(sc),min(mm)),upper=c(max(sc),max(mm)),control=cl)
bH <- c(w$par,hS,w$value)
names(bH) <- c("scale","mag_var","richness","loglikelihood")
return(bH)
}

# find best lognormal occupancy distribution
optimize_lognormal_occupancy <- function(finds,ncoll)	{
observed <- fisher_plot(finds)
oS <- length(finds[finds>0])				# observed taxa
minS <- length(finds)						# observed + inferred taxa
observed <- fisher_plot(finds)
p0 <- accio_initial_occupancy_lognormal(finds,observed,oS,minS,ncoll)	# pO[1]=r; p0[1]=magnitude; p0[3]=richness
b0 <- accio_min_max_lognormal_occupancy_params (finds,observed,oS,ncoll)
mm <- b0[1:2]
sc <- b0[3:4]
# make sure that lower & upper bounds contain starting values!
if (mm[1] > p0[2])	{
	mm[1] <- 1+((p0[2]-1)/2)
	}	else if (mm[2] < p0[2])	{
	mm[2] <- 2*p0[2]
	}
if (sc[1] > p0[1])	{
	sc[1] <- p0[1]/2
	}	else if (sc[2] < p0[1])	{
	sc[2] <- p0[1]
	}
pa <- 1
pz <- span <- 4		# span of richnesses to consider
stS <- minS
enS <- stS+((span-1)*minS)
incr <- floor(enS/span)
cl <- list(fnscale=-1)
peak <- 0
last_hS <- c(0,1,2,3)
while (incr>0)	{
	hS <- round(seq(stS,enS,by=incr),0)
	if (sum(hS %in% last_hS) == length(hS))	break
	if (pz>length(hS))	pz <- length(hS)
	results <- sapply(hS[pa:pz],optimize_lognormal_occupancy_given_hS,oS=oS,observed=observed,ncoll=ncoll,p0=p0,sc=sc,mm=mm,debug=0)
	if (pa==2)
		results <-cbind(old_result_1,results)
	if (pz<length(hS))
		results <-cbind(results,old_result_2)
	mlnl <- max(results[4,])
	mlSc <- match(mlnl,results[4,])
	mlS <- hS[mlSc]
	bH <- results[,mlSc]
	spn <- length(results[1,])
	if (incr>1)	{
		# if runs so far are still producing higher likelihoods at higher richnesses:
		if (mlSc==spn)	{
			if (peak==0) {
				# if we have not yet found a peak, then keep increasing richness
				pa <- 2		# make note of best hypothesis
				old_result_1 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				stS <- hS[spn]
				enS <- stS+((spn-1)*minS)
				}	else	{
				stS <- hS[mlSc-1]
				enS <- hS[mlSc]
				old_result_1 <- results[,mlSc-1]	# set aside stats so as to avoid recalculating
				old_result_2 <- results[,mlSc]	# set aside stats so as to avoid recalculating
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
#				stS <- hS[spn]-(span*incr)
				}	# end case where last number is best, but this is after finding a peak.
			} else if (mlSc==1)	{
			# if first richness is the best
			peak <- 1
			stS <- hS[mlSc]
			enS <- hS[mlSc+1]
			old_result_1 <- results[,1]
			pa <- 2
			old_result_2 <- results[,2]
			pz <- length(seq(stS,enS,by=incr))-1
			if (incr==1)	incr <- 0
			if (span<incr)	{
				incr <- (enS-stS)/(span-1)
				}	else {
				span <- enS-stS
				incr <- 1
				}
			} else {
			# if a richness in the middle is best & we still can find a better one
			peak <- 1
			hS2 <- c(mlS-1,mlS+1)
			results2 <- sapply(hS2,optimize_lognormal_occupancy_given_hS,oS=oS,observed=observed,ncoll=ncoll,p0=p0,sc=sc,mm=mm,debug=0)
			if (results2[4,1]>mlnl && results2[4,1]>results2[4,2])	{
			# start just above 2nd best richness and go up to best
				stS <- hS[mlSc-1]
				enS <- hS2[1]
				old_result_1 <- results[,mlSc-1]
				old_result_2 <- results2[,1]
				if (span<incr)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else if (results2[4,2]>mlnl && results2[4,2]>results2[4,1]) {
			# start at best richness and go just below 2nd best
				stS <- hS2[2]
				enS <- hS[mlSc+1]
				old_result_1 <- results2[,2]
				old_result_2 <- results[,mlSc+1]
				if (incr>span)	{
					incr <- (enS-stS)/(span-1)
					}	else {
					span <- enS-stS
					incr <- 1
					}
				pa <- 2
				pz <- length(seq(stS,enS,by=incr))-1
				if (pz < pa)	incr <- 1
				}	else {
				# we already had the best, so just end it
				incr <- 0
				}
			pa <- 2
			pz <- spn-1
			}
		# end case where we have a better richness in middle somewhere.  
		last_hS <- hS
		}	else	{
		incr <- 0
		}
#	print(hS)		# for debugging
#	if (incr>0)	print(round(seq(stS,enS,by=incr),0))	# for debugging
	}
lognormal_AICc <- round(modified_AIC(mlnl,3,sum(finds)),2)
bH <- c(bH[1],bH[2],bH[3],round(bH[4],2),lognormal_AICc)
names(bH)[5] <- "AICc"
return(bH)
}

##### DIVERSIFICATION #####
accio_best_diversification_given_sampling <- function(pmiss,S1,two_timer,gap_filler,continuous)	{
# rate gives the proportion NOT shared between bins
#		originations or extinctions for historical data
# pmiss: typical probability of failing to sample a taxon;
# S1: observed standing richness (no unsampled range-throughs);
# two_timer: number shared with the prior interval;
# gap_filler: number of unsampled taxa spanning a gap (synoptic - observed richness);
# continuous: if T, then continuous diversification assumed; otherwise, pulsed.
if (continuous==TRUE)	{
	if ((two_timer+gap_filler)>0)	{
		rate <- -log((two_timer+gap_filler)/S1);
		}	else {
		rate <- -log(pmiss/S1);
		}
	}	else {
	rate <- 1-(two_timer/S1);
	}
# we expect to sample (1-pmiss) of the taxa.
minrate <- 0;
maxrate <- rate;	# no point in considering a rate higher than face value

cl <- list(fnscale=-1);
if ((two_timer+gap_filler)<S1 && rate>0)	{
	accio <- optim(rate,fn=likelihood_diversification_rate_given_sampling,method="L-BFGS-B",pmiss=pmiss,S1=S1,two_timer=two_timer,gap_filler=gap_filler,continuous=continuous,lower=0,upper=maxrate,control=cl)
	best_diversification <- max(0,accio$par)
	L_best_diversification <- accio$value
	}	else {
	best_diversification <- 0.0
	L_best_diversification <- 1.0
	}
return(c(best_diversification,log(L_best_diversification)))
}

likelihood_diversification_rate_given_sampling <- function(rate,pmiss,S1,two_timer,gap_filler,continuous)	{
# rate: per-lineage rate
#	poisson if continuous==TRUE; binomial if continuous==FALSE
# pmiss: probability that an unobserved taxon is present-but-unsampled rather than non-existent
# S1: observed richness in "test" interval
# two_timer: shared richness
# gap_filler: number of taxa inferred to be present because they are there before and after
# continuous: TRUE for continuous turnover, FALSE for discrete turnover
if (continuous==TRUE)	{
	freq <- Poisson_rate_to_probability(rate)
	}	else {
	freq <- rate
	}	# get expected frequency of shared taxa
hS <- seq(two_timer+gap_filler,S1,by=1)	# range of possible shared taxa; max=S1, min=observed+directly inferred
pt <- prob_observing_n_shared_given_diversification_sampling_and_hypothesized_shared(n=two_timer,S1=S1,S2=hS,freq_turn=freq,pmiss=pmiss)

return(max(10^-320,sum(pt)))
}

# function to calculate probability of observing n of N shared taxa give sampling
#	x probability of N survivors given extinction/origination/beta & S1 original taxa
prob_observing_n_shared_given_diversification_sampling_and_hypothesized_shared <- function(n,S1,S2,freq_turn,pmiss)	{
# S1 = observed in bin
# S2 = hypothesized in other bin
# n = observed shared
# freq_turn=expected n=(S1-S2)S1
# pfind = prob finding taxon
#return(dbinom(S1-S2,S1,freq_turn)*dbinom(n,S2,pfind))
return(dbinom(S1-S2,S1,freq_turn)*pmiss^(S2-n))
}

##### ROUTINES TO WRITE CHARACTER MATRIX INFORMATION IN NEXUS FILE #######
scribio_rev_bayes_nexus_file_from_character_matrix <- function(ch_matrix,new_file_name,no_states,UNKNOWN=-11,INAP=-22)	{
# ch_matrix: character matrix being printed to file.
# new_file_name: name of nexus file, including the directory to which it should be printed.
# no_states: number of states; note that if only a portion of the taxa are used, then sometimes no otus have some states
if (!is.matrix(ch_matrix))	{
	ch_matrix <- data.frame(ch=ch_matrix);
	}
notu <- nrow(ch_matrix);
taxon_names <- rownames(ch_matrix);
nchars <- ncol(ch_matrix);

nexus_file_content <- c();
nexus_file_content <- rbind("#NEXUS","","BEGIN DATA;")
nexus_file_content <- rbind(nexus_file_content,paste("	DIMENSIONS  NTAX=",notu," NCHAR=",nchars,";",sep=""));
if (no_states<10) {
	state_symbols <- " ";
	for (st in 1:no_states)
		state_symbols <- paste(state_symbols,st-1,sep=" ");	
	} else	{
	mxl <- no_states-10;
	letter_states <- LETTERS[!LETTERS %in% c("I","O")][1:mxl]
	all_states <- c(0:9,letter_states);
	state_symbols <- paste(all_states,collapse=" ");
	}
nexus_file_content <- rbind(nexus_file_content,paste("	FORMAT DATATYPE = STANDARD RESPECTCASE GAP = - MISSING = ? SYMBOLS = \"",state_symbols,"\";"));
nexus_file_content <- rbind(nexus_file_content,"	MATRIX");

string_to_count <- taxon_names;
name_lengths <- sapply(string_to_count,count_characters_in_string);
max_name_length <- max(name_lengths);
need_quotes <- c(".","(",")","[","]");
for (nn in 1:notu)	{
	test_name <- strsplit(taxon_names[nn],split="",fixed=TRUE)[[1]]
	if (sum(test_name %in% need_quotes)==0)	{
		taxon <- gsub(" ","_",taxon_names[nn]);
		} else	{
		taxon <- paste("\"",taxon_names[nn],"\"",sep="");
		name_lengths[nn] <- name_lengths[nn]+2;
		}
	this_line <- paste("\t",taxon,paste(rep(" ",(5+(max_name_length-name_lengths[nn]))),collapse=""),sep="");
	otu_code <- c();
	for (ch in 1:nchars)	{
		if (ch_matrix[nn,ch]>=0 && ch_matrix[nn,ch]<=9)	{
			otu_code <- paste(otu_code,ch_matrix[nn,ch],sep="");
			} else if (ch_matrix[nn,ch]>9)	{
			otu_code <- paste(otu_code,all_states[1+ch_matrix[nn,ch]],sep="");	# note: we need +1 because of state 0
			} else if (ch_matrix[nn,ch]==UNKNOWN)	{
			otu_code <- paste(otu_code,"?",sep="");
			} else if (ch_matrix[nn,ch]==INAP)	{
			otu_code <- paste(otu_code,"-",sep="");
			} else if (ch_matrix[nn,ch]<0)	{
			polystates <- strsplit(as.character(ch_matrix[nn,ch]),split="",fixed=TRUE)[[1]];
			polystates <- as.numeric(polystates[polystates!="-"]);
			otu_code <- paste(otu_code,ravel_polymorph_for_file(polystates),sep="");
			}
		}
	nexus_file_content <- rbind(nexus_file_content,paste(this_line,otu_code,sep=""));
	}

nexus_file_content <- rbind(nexus_file_content,";");
nexus_file_content <- rbind(nexus_file_content,"END;");
#nexus_file_content <- rbind(nexus_file_content,"begin mrbayes;");
#nexus_file_content <- rbind(nexus_file_content,"	set autoclose=yes nowarn=yes;");
#nexus_file_content <- rbind(nexus_file_content,"	lset nst=6 rates=invgamma;");
#nexus_file_content <- rbind(nexus_file_content,"	unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all); ");
#nexus_file_content <- rbind(nexus_file_content,"	prset applyto=(all) ratepr=variable;");
#nexus_file_content <- rbind(nexus_file_content,"	mcmcp ngen= 100000000 relburnin=yes burninfrac=0.25 printfreq=10000  samplefreq=10000 nchains=4 savebrlens=yes;");
#nexus_file_content <- rbind(nexus_file_content,"	mcmc;");
#nexus_file_content <- rbind(nexus_file_content,"	sumt;");
#nexus_file_content <- rbind(nexus_file_content,"end;");
write(nexus_file_content,file=new_file_name);
}

## Deal with Polymorphs	##
switch_letter_state_to_numeric <- function(state)  {
# 2017-10-09: now will pass numeric characters through unchanged
# 2019-01-25: simplified greatly!
if (state > 9)	{
	state <- toupper(state)
	poss_letter_states <- toupper(letters[!letters %in% c("i","o")]);
	return(9+match(state,poss_letter_states));
	} else	{
	return(state);
	}
}

ravel_polymorph <- function(polystates)	{
polystates <- sort(polystates,decreasing = TRUE);
polym <- polystates[1];
for (st in 2:length(polystates))	polym <- (10*polym)+polystates[st]
return(-1*polym)
}

ravel_polymorph_for_file <- function(polystates)	{
polystates <- sort(polystates,decreasing = FALSE);
return(paste("(",paste(polystates,collapse=""),")",sep=""));
}

##### ELEMENTARY DIVERGENCE TIME ESTIMATION	#####
simple_likelihood_divergence <- function(bound_1,bound_2,psi,cutoff=exp(-4),precision=0.01)	{
# bound_1: oldest taxon;
# bound_2: 2nd oldest taxon;
# psi: probability of sampling a taxon per time unit 1;
# cutoff: area under likelihood curve at which we fix the divergence (null: 1 unit support)
# precision: temporal precision for estimating divergence, with expected finds = psi*precision
bound_1 <- abs(bound_1);
bound_2 <- abs(bound_2);
min_divergence_date <- divergence_date <- max(bound_1,bound_2);
gap_ldf <- c();				# likelihood density function for gap
pgap <- 1;
while (pgap > 10^-10)	{
	pgap <- dpois(0,psi*((divergence_date-bound_1)+(divergence_date-bound_2)));
	gap_ldf <- c(gap_ldf,pgap);
	divergence_date <- divergence_date+precision;
	}
return(min_divergence_date+(precision*(1+sum(gap_ldf>cutoff))));
#return(min_divergence_date+(precision*(1+sum(cumsum(gap_ldf/sum(gap_ldf))<0.5))));
#return(min_divergence_date+(precision*(1+sum(cumsum(gap_pdf)<0.5))));
}

simple_probability_divergence <- function(bound_1,bound_2,phi,psi,cutoff=exp(-8),precision=0.01)	{
# bound_1: oldest taxon;
# bound_2: 2nd oldest taxon;
# phi: probability of sampled sister taxon (lineage or clade) arising;
# psi: probability of sampling a taxon per time unit 1;
# cutoff: area under likelihood curve at which we fix the divergence.  (Default: logP = -1)
# precision: temporal precision for estimating divergence, with expected finds = psi*precision
bound_1 <- abs(bound_1);
bound_2 <- abs(bound_2);
min_divergence_date <- divergence_date <- max(bound_1,bound_2);
gap_pdf <- c();
pgap <- lgap <- 1;
while (exp(log(pgap)+log(lgap)) > 10^-10)	{
#	gap_pdf <- c(gap_pdf,dpois(0,phi*((divergence_date-bound_1)+(divergence_date-bound_2))));
	pgap <- (1-phi)^((divergence_date-bound_1)+(divergence_date-bound_2));
	lgap <- dpois(0,psi*((divergence_date-bound_1)+(divergence_date-bound_2)));
	gap_pdf <- c(gap_pdf,exp(log(pgap)+log(lgap)));
	divergence_date <- divergence_date+precision;
	}
return(min_divergence_date+(precision*(1+(sum(gap_pdf>cutoff)))));
}

##### ROUTINES FOR A PRIOR EVALUATION OF POSSIBLE CHARACTER MODELS #####
revelio_possible_driven_trends <- function(analysis_name,nexus_file_name,nexus_file_directory,time_scale,strat_data_file,UNKNOWN=-11,INAP=-22)	{
new_file_name <- paste(nexus_file_directory,nexus_file_name,sep="");
basic_data <- accio_data_from_nexus_file(nexus_file_name=new_file_name);
ch_matrix <- basic_data$Matrix;
outgroup <- as.numeric(basic_data$Outgroup);
otus <- basic_data$OTUs;
ch_matrix_red <- ch_matrix[!rownames(ch_matrix) %in% otus[outgroup],];
finds_per_bin_data <- read.csv(strat_data_file,header=T,stringsAsFactors=hell_no);
if (is.na(match("Species",colnames(finds_per_bin_data))))	{
	ccns <- colnames(finds_per_bin_data);
	ccns[1] <- "Species";
	colnames(finds_per_bin_data) <- ccns;
	}
rounded_finds_per_bin <- finds_per_bin_data[!finds_per_bin_data$Species %in% otus[outgroup],];
rownames(rounded_finds_per_bin) <- rounded_finds_per_bin$Species;
rounded_finds_per_bin$Species <- NULL;
notu <- nrow(rounded_finds_per_bin);
for (i in 1:notu)	{
	rounded_finds_per_bin[i,] <- round(cumsum(as.numeric(rounded_finds_per_bin[i,]))+0.001,0);
	}
bin_names <- colnames(rounded_finds_per_bin);
bin_names <- bin_names[colSums(rounded_finds_per_bin)>0]
rounded_finds_per_bin <- rounded_finds_per_bin[,colnames(rounded_finds_per_bin) %in% bin_names];
nbins <- length(bin_names);
fa_bin_no <- fa_bin <- c();
for (i in 1:notu)	{
	fa_bin_no <- c(fa_bin_no,min((1:nbins)[rounded_finds_per_bin[i,]>0]));
	}
fa_bin <- bin_names[fa_bin_no];

nchars <- ncol(ch_matrix_red);
nstates <- basic_data$States;
assessed_characters <- significant_cases <- kendall_pval <- kendall_tau <- c();
for (nch in 1:nchars)	{
	char_states <- ch_matrix_red[,nch];
	char_binning <- fa_bin_no*((char_states != UNKNOWN)*(char_states != INAP));	# first appearances of states in different taxa
	char_binning <- char_binning[char_binning!=0];								# cull out species with unknown or inapplicable
	unique_bins <- sort(unique(char_binning));
	char_states <- char_states[char_states!=UNKNOWN];							# cull out species with unknown
	char_states <- char_states[char_states!=INAP];								# cull out species with inapplicable

	if (sum(char_states<0)>0)	{
		poly_cases <- (1:length(char_states))[char_states<0];
		poly_fas <- char_binning[char_states<0];
		for (pf in 1:length(poly_fas))	{
			polystates <- unravel_polymorph(char_states[poly_cases[pf]]);
			char_states[poly_cases[pf]] <- polystates[1];
			char_states <- c(char_states,polystates[2]);
			char_binning <- c(char_binning,poly_fas[pf]);
			}
		}
	names(char_states) <- names(char_binning);
	state_nos <- sort(unique(char_states[char_states>=0]));
	this_chars_fas <- fa_bin_no[match(names(char_states),rownames(ch_matrix_red))];
	state_fas <- c();
	for (nst in 1:nstates[nch])
		state_fas <- c(state_fas,min(char_binning[char_states==state_nos[nst]]));
	state_nos <- state_nos[order(state_fas)];
	
	if (length(unique(char_states))==1)	{
		kendall_pval <- c(kendall_pval,1.0);
		kendall_tau <- c(kendall_tau,0.00);
		} else	{
		kendall_pval <- c(kendall_pval,cor.test(x=char_binning,y=char_states,method="kendall")$p.value);
		kendall_tau <- c(kendall_tau,cor.test(x=char_binning,y=char_states,method="kendall")$estimate);
		}
	
	time_state_matrix <- array(0,dim=c(nbins,nstates[nch]));
	for (nn in 1:length(char_states))	{
		if (char_states[nn]<0)	{
			polystates <- unravel_polymorph(poly=char_states[nn]);
			time_state_matrix[this_chars_fas[nn],match(polystates,state_nos)] <- time_state_matrix[this_chars_fas[nn],match(polystates,state_nos)]+(1/length(polystates));
			} else	{
			time_state_matrix[this_chars_fas[nn],match(char_states[nn],state_nos)] <- time_state_matrix[this_chars_fas[nn],match(char_states[nn],state_nos)]+1;
			}
		}
	bin_cases <- rowSums(time_state_matrix);
	mann_whitney_results <- c();
	comparisons <- 0;
	for (bn in 1:(nrow(time_state_matrix)-1))	{
		if (bin_cases[bn] > 0 && bin_cases[bn+1] > 0)	{
			data_v <- category <- c();
			for (bnb in 0:1)	{
				for (bnst in 1:nstates[nch])	{
					data_v <- c(data_v,rep(state_nos[bnst],round(time_state_matrix[bn+bnb,bnst])));
					category <- c(category,rep(bnb,round(time_state_matrix[bn+bnb,bnst])));
					}
				}
			comparisons <- comparisons + 1;
			mann_whitney_results <- rbind(mann_whitney_results,mann_whitney(data_v,category));
			}
		}
	if (comparisons>0)	{
		assessed_characters <- c(assessed_characters,nch);
		significant_cases <- c(significant_cases,sum(mann_whitney_results$pval<=0.10));
		} else	{
		significant_cases <- c(significant_cases,0);
		}
	}
kendall_summary <- data.frame(nch=as.numeric(1:nchars),tau=as.numeric(kendall_tau),pval=as.numeric(kendall_pval),stringsAsFactors=hell_no);
kendall_summary <- kendall_summary[order(kendall_summary$pval),];
trend_summary <- data.frame(nch=as.numeric(1:nchars),states=as.numeric(nstates),tau=as.numeric(kendall_tau),pval=as.numeric(kendall_pval),MW_sign=as.numeric(significant_cases),stringsAsFactors=hell_no);
trend_summary[order(trend_summary$states,trend_summary$pval),];

#hist(kendall_summary$pval,breaks=((0:100)/100))

cbind(assessed_characters[significant_cases==max(significant_cases)],significant_cases[significant_cases==max(significant_cases)])
#(1:nchars)[significant_cases==max(significant_cases)]
}

#median(c(0.136, 0.008, 0.078, 0.054, 0.156, 0.235, 0.855, 0.133, 0.065, 0.016, 0.169, 0.903, 0.013, 0.160, 0.591, 0.021, 0.286, 0.124, 0.006, 0.727, 0.142, 0.297,0.115, 0.034, 0.025, 0.074, 0.810, 0.100, 0.373, 0.019, 0.095, 0.590, 0.168, 0.137, 0.088, 0.021, 0.339, 0.350, 0.158, 0.294, 0.132, 0.163, 0.084, 0.034,0.112, 0.058, 0.731, 0.162, 0.359, 0.182, 1.093, 0.066, 0.074, 0.359, 0.488, 0.166, 0.385, 0.450, 0.504, 0.300, 0.293, 0.261, 0.654, 0.025, 0.219, 0.012,0.282, 0.099, 0.184, 0.126, 0.072, 0.082, 0.287, 0.684, 0.020, 0.398, 0.298, 0.444, 0.153, 0.126, 0.048, 0.032, 0.305, 0.015, 0.232, 0.328, 0.121, 0.031,0.748, 0.001, 0.023, 0.113, 0.312, 0.018, 0.636, 0.032, 0.260, 0.116, 0.033, 0.642))
rescore_states_by_first_appearances <- function(otu_states,otu_fas,UNKNOWN=-11,INAP=-22)	{
# remove taxa with inapplicable & unknown conditions
otus <- (1:length(otu_states))[!otu_states %in% c(UNKNOWN,INAP)];
otu_fas_relv <- otu_fas[!otu_states %in% c(UNKNOWN,INAP)];
otu_states_relv <- otu_states[!otu_states %in% c(UNKNOWN,INAP)];

# deal with polymorphics:
poly_otus <- (1:length(otu_states_relv))[otu_states_relv<0];
poly_otus_orig <- otus[otu_states_relv<0];
poly_scores <- otu_states_relv[otu_states_relv<0];
if (length(poly_otus)>0)	{
	for (po in 1:length(poly_otus))	{
		psts <- unravel_polymorph(poly_scores[po]);	# observed_polymorphic states
		otu_states_relv[poly_otus[po]] <- psts[1];
		otu_states_relv <- c(otu_states_relv,psts[2:length(psts)]);
		otu_fas_relv <- c(otu_fas_relv,rep(otu_fas_relv[poly_otus[po]],length(psts)-1));
		}
	}
relv_states <- sort(unique(otu_states_relv));
state_appearance <- c();
for (rs in 1:length(relv_states))
	state_appearance <- c(state_appearance,min(otu_fas_relv[otu_states_relv==relv_states[rs]]));

new_coding <- otu_states;
if (sum(state_appearance==state_appearance[order(state_appearance)]) < length(state_appearance))	{
	states_by_appearances <- relv_states[order(state_appearance)];
	for (rs in 1:length(relv_states))
		new_coding[otu_states==relv_states[rs]] <- states_by_appearances[rs];
	po <- 0;
	while (po < length(poly_otus))	{
		po <- po+1;
		new_coding[poly_otus_orig[po]] <- ravel_polymorph(relv_states[match(unravel_polymorph(otu_states[poly_otus_orig[po]]),states_by_appearances)])
		}
	}
return(new_coding);
}

##### PROBABILITY 101 ####
probability_to_Poisson_rate <- function(pn)	{
return(-1*log(1-pn))
}

Poisson_rate_to_probability <- function(expectation)	{
return(1-exp(-expectation))
}

##### KLUGES! #####
# clean NA from matrix
evanesco_na_from_matrix <- function(data, replacement="")  {
for (i in 1:ncol(data))	{
	if(sum(is.na(data[,i]))>0)	{
		duds <- (1:nrow(data))[is.na(data[,i])]
		data[duds,i] <- replacement
		}
	}
return(data)
}

# clean NA from vector
evanesco_na_from_vector <- function(data, replacement="")	{
if(sum(is.na(data))>0)	{
	duds <- (1:length(data))[is.na(data)]
	data[duds] <- replacement
	}
return(data)
}

# Change String To Title Case 
transmogrify_to_title_case <- function(name) {
name_part <- strsplit(name, " ")[[1]]
return (paste(toupper(substring(name_part, 1,1)), substring(name_part, 2),sep="", collapse=" "));
}

# Change ü to ue, etc.
transmogrify_diacritics <- function(funky_text)	{
j <- strsplit(as.character(funky_text),split="",fixed=TRUE)[[1]];
eek <- c("à","á","â","ã","ā","ă","ȧ","ä","ả","å","ǎ","ȁ","ȃ","ą","ạ","ḁ","ẚ","ầ","ấ","ẫ","ẩ","ằ","ắ","ẵ","ẳ","ǡ","ǟ","ǻ","ậ","ặ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "a";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "A";
eek <- c("æ","ǽ","ǣ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "ae";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Ae";
eek <- c("ć","ĉ","ċ","č","ƈ","ç","ḉ","ȼ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "c";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "C";
eek <- c("ḋ","ɗ","ḍ","ḏ","ḑ","ḓ","ď");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "d";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "D";
eek <- c("è","é","ê","ẽ","ē","ĕ","ė","ë","ẻ","ě","ȅ","ȇ","ẹ","ȩ","ę","ḙ","ḛ","ề","ế","ễ","ể","ḕ","ḗ","ệ","ḝ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "e";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "E";
eek <- c("ǵ","ĝ","ḡ","ğ","ġ","ǧ","ɠ","ģ","ǥ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "g";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "G";
eek <- c("ĥ","ḣ","ȟ","ḥ","ḩ","ḫ","ẖ","ħ","ⱨ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "h";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "H";
eek <- c("ì","í","î","ĩ","ī","ĭ","ı","ï","ỉ","ǐ","ị","į","ȉ","ȋ","ḭ","ɨ","ḯ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "i";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "I";
eek <- c("ĵ","ǰ","ȷ","ɉ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "j";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "J";
eek <- c("ḱ","ǩ","ḵ","ƙ","ḳ","ĸ","ⱪ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "k";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "K";
eek <- c("ĺ","ľ","ŀ","ł","ƚ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "l";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "L";
eek <- c("ḿ","ṁ","ṃ","ɱ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "m";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "M";
eek <- c("ǹ","ń","ñ","ṅ","ň");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "n";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "N";
eek <- c("ò","ó","ô","õ","ō","ŏ","ȯ","ö","ỏ","ő","ǒ","ȍ","ȏ","ơ","ǫ","ọ","ɵ","ø");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "o";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "O";
eek <- c("ṕ","ṗ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "p";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "P";
eek <- c("ŕ","ṙ","ř","ȑ","ȓ","ṛ","ŗ","ṟ","ṝ","ɍ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "r";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "R";
eek <- c("ś","ŝ","ṡ","š","ṣ","ș","ş","ȿ","ṥ","ṧ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "s";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "S";
eek <- c("ß");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "ss";
eek <- c("ṫ","ẗ","ť","ƫ","ṭ","ț","ţ","ṱ","ṯ","ŧ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "t";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "T";
eek <- c("ù","ú","û","ũ","ū","ŭ","ü","ủ","ű","ǔ","ȕ","ȗ","ụ","ṳ","ų","ṷ","ṵ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "u";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "U";
eek <- c("ṽ","ṿ","ⱱ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "v";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "V";
eek <- c("ẁ","ẃ","ŵ","ẇ","ẅ","ẘ","ẉ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "w";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "W";
eek <- c("ẋ","ẍ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "x";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "X";
eek <- c("ỳ","ý","ŷ","ȳ","ẏ","ÿ","ỷ","ẙ","ɏ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "y";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Y";
eek <- c("ź","ẑ","ż","ž","ȥ","ẓ","ẕ","ƶ");
eekeek <- (1:length(j))[j %in% eek];
j[eekeek] <- "z";
eekeek <- (1:length(j))[j %in% toupper(eek)];
j[eekeek] <- "Z";
j <- j[j!="̈"];
j <- j[j!="̀"];
j <- j[j!="̌"];
return(paste(j,collapse=""));
}

scourgify_web_text <- function(web_text)	{
web_text <- gsub("–","-",web_text);
web_text <- gsub("“","\"",web_text);
web_text <- gsub("”","\"",web_text);
web_text <- gsub("‘","\'",web_text);
web_text <- gsub("’","\'",web_text);
web_text <- gsub("‚Äì","–",web_text);
web_text <- gsub("‚Äî","-",web_text);
web_text <- gsub("‚Äú","“",web_text);
web_text <- gsub("‚Äù","”",web_text);
web_text <- gsub("‚Äò","‘",web_text);
web_text <- gsub("‚Äô","’",web_text);
web_text <- gsub("é","e",web_text);
web_text <- gsub("√©","e",web_text);
#web_text <- gsub("ö","&ouml;",web_text);
#web_text <- gsub("ü","&uuml;",web_text);
web_text <- gsub("ü","u",web_text);
web_text <- gsub("√º","u",web_text);
web_text <- gsub("ý","y",web_text);
web_text <- gsub("√Ω","y",web_text);
return(web_text);
}

scourgify_web_text_vanilla <- function(web_text)	{
web_text <- gsub("–","-",web_text);
web_text <- gsub("“","\"",web_text);
web_text <- gsub("”","\"",web_text);
web_text <- gsub("‘","\'",web_text);
web_text <- gsub("’","\'",web_text);
web_text <- gsub("‚Äì","–",web_text);
web_text <- gsub("‚Äî","-",web_text);
web_text <- gsub("‚Äú","“",web_text);
web_text <- gsub("‚Äù","”",web_text);
web_text <- gsub("‚Äò","‘",web_text);
web_text <- gsub("‚Äô","’",web_text);
web_text <- gsub("â€™","’",web_text);
web_text <- gsub("√°","a",web_text);
web_text <- gsub("√©","e",web_text);
web_text <- gsub("&#321;","L",web_text);
web_text <- gsub("≈ç","o",web_text);
web_text <- transmogrify_diacritics(funky_text=web_text);
web_text <- gsub("√º","u",web_text);
web_text <- gsub("√∫","u",web_text);
web_text <- gsub("√Ω","y",web_text);
web_text <- gsub("\t","     ",web_text);
return(web_text);
}

#find_me <- paleodb_clean_member_basic_no_rock[need_info_member[3]]
matrix_row_match <- function(find_me,examined_matrix)	{
xxx <- which(examined_matrix==find_me,arr.ind=TRUE);
#examined_matrix[unique(xxx[,1]),]
if (nrow(xxx)==0)	{
	return(-1);
	} else if (length(unique(xxx[,1]))==1)	{
	return(xxx[1,1])
	} else {
	vvv <- unique(xxx[,1]);
	mx <- 0;
	mj <- 0;
	for (j in 1:length(vvv))	{
		if (sum(vvv[j]==xxx[,1])>mx)	{
			mj <- j;
			mx <- sum(vvv[j]==xxx[,1])
			} else if (sum(vvv[j]==xxx[,1])==mx)	{
			mj <- c(mj,j)
			}
		}
	if (length(mj)==1)	{
		return(vvv[mj]);
		} else {
		return(-2);
		}
	}
}

#character_array <- strsplit(rock_unit_name,split = "");
simplify2vector <- function(character_array)	{
return(simplify2array(character_array)[,1])
}

specialis_revelio_file_type <- function(filename)	{
file_info <- stringr::str_split(filename,pattern="")[[1]];
f_i <- length(file_info);
f_i_s <- 1+(1:f_i)[file_info %in% "."];
file_type <- c();
for (i in f_i_s:f_i)	file_type <- paste(file_type,file_info[i],sep="");
return(file_type);
}

insert_cell_into_vector_x <- function(x,new_value,cell_no)	{
# routine to put take vector x and add a new value at some point (cell_no)
# x = c(0,1,2,3,4,5), new_value=10, cell_no=3 gives x=c(0,1,10,2,3,4,5)
if (cell_no==1)	{
	return(c(new_value,x));
	} else if (cell_no==(length(x)+1))	{
	return(c(x,new_value));
	} else	{
	return(c(x[1:(cell_no-1)],new_value,x[cell_no:length(x)]));
	}
}

insert_row_into_matrix_x <- function(x,new_row,row_no)	{
# routine to put take matrix x and add a new row at that bumps other rows down one
# x = c(0,1,2,3,4,5), new_value=10, cell_no=3 gives x=c(0,1,10,2,3,4,5)
if (row_no==1)	{
	xx <- rbind(new_row,x);
	colnames(xx) <- colnames(x);
	return(xx);
	} else if (row_no==(nrow(x)+1))	{
	return(rbind(x,new_row));
	} else	{
	return(rbind(x[1:(row_no-1),],new_row,x[(row_no:nrow(x)),]));
#	return(c(x[1:(cell_no-1)],new_value,x[cell_no:length(x)]));
	}
}
