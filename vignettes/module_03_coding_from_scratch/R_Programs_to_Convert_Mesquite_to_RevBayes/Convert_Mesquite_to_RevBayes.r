common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
source(paste(common_source_folder,"RevBayes_Setup.r",sep=""));

# This will take the Mesquite file and break it down into three matrices for 2-, 3- & 4-state characters.
# None of the notes or information about taxa, characters & states are retained; these choke RevBayes
analysis_name <- "Cinctans_PS_2019";
mesquite_to_RevBayes <- diffindo_character_matrix_by_state_numbers_and_other_partitions(analysis_name);

set_wdir <- "/Users/peterjwagner/Documents/RevBayes_Projects/";
scribio_RevBayes_scripts_from_chosen_nexus_file_and_existing_FBD_script_and_data(analysis_name=analysis_name,set_wdir=set_wdir);

end()