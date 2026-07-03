# VCF Table Viewer - CCDI & CGC Platform
# October 2025 

library(shiny) 
library(shinydashboard)
library(shinyFiles)
library(shinyalert)
library(bslib)
library(fs)

library(tidyr)
library(dplyr)
library(ggplot2)

library(DT)
library(data.table) # for melt
library(stringr)

library(vcfR)
library(igvShiny)
library(rtracklayer)
library(GenomicAlignments)

##############################################################
####-------------Customization Section--------------------####
##############################################################

# needed for reading in the legend HTML file vcf_field_descriptions.html
addResourcePath("tmpuser", getwd()) 

global <- reactiveValues(sarekDir = "./",  
                         vcfDir = "./sbgenomics/project-files/", 
                         bamDir = "/mnt/BW-Data/recalibrated/")

# Manifest file with the following headers:
#VCFFileName, caller, ParticipantID, SampleID, BAMFileName (optional), StudyID (optional)
#manifest <- read.csv("sbgenomics/project-files/CCDI_Manifest_Example.csv")
#manifest <- manifest |> arrange(StudyID, ParticipantID, SampleID)
# get caller column from file names
#manifest <- manifest |> mutate(caller = str_select(FileName))
# a6a77776-f50a-4630-bdcf-631b7e7e51d0.vardict_somatic.norm.annot.public.vcf.gz
# 65377817-5b14-4314-a87b-5eb4bae3757c.mutect2_somatic.norm.annot.public.vcf.gz
# 9bf1f6d4-29c9-4f68-88e1-4246d9ce16e0.consensus_somatic.norm.annot.public.vcf.gz
# cc060cd2-3f50-4e33-95bb-27d81619d808.lancet_somatic.norm.annot.public.vcf.gz
# 5d9a45fe-a6ed-4619-8a2b-aa69614e8e03.strelka2_somatic.norm.annot.public.vcf.gz

## Dropdown items
#study_list <- unique(manifest$StudyID)
#subject_list <- unique(manifest$ParticipantID)
#sample_list <- unique(manifest$SampleID)

# list of callers
#callers <- c("consensus", "strelka2", "mutect2", "lancet", "vardict")
#callers <- unique(manifest$Caller)

# list of filtering levels
# Left out due to size constraints:
#  - Annotation: full annotated VCF produced by sarek
#  - Region: filter VCF by GIAB mappable region
# 
# 1. Population:        filter by population allele frequency < 0.01
# 2. Mutation:          filter by significant mutations
# 3. ML Driver Genes:   filter by myeloid cancer driver genes
# 4. Genes of Interest: filter all genes of interest out of region-filtered VCF
filters <- c("Population", "Mutation", "ML Driver Genes", "Genes of Interest")
filterNames <- c("ann.rtgfilt.popfilt", "ann.rtgfilt.popfilt.sigmut", 
                 "ann.rtgfilt.popfilt.sigmut.genesmut", "ann.rtgfilt.allgenes")
names(filterNames) <- filters

## lists of important genes to highlight in the table
gene_lists <- read.csv("./sbgenomics/project-files/Gene_lists.txt", header = T, sep = "\t")

##############################################################

printf <- function(...) print(noquote(sprintf(...))) # used with igvShiny

# for highlighting mutation severity columns
color_gradient <- function(dt, column_name, gradient_colors = c("#FF6666", "#DDDDDD")) {
  col_func <- colorRampPalette(gradient_colors)
  dt |> 
    formatStyle(column_name, 
                backgroundColor = styleEqual(
                  sort(unique(dt$x$data[[column_name]]), decreasing = TRUE),
                  col_func(length(unique(dt$x$data[[column_name]])))
                )
    ) 
}

# user interface 
ui <- dashboardPage(
      
      skin = "blue",
      dashboardHeader(title = "VCF Table Viewer - CCDI", titleWidth = 300),

      dashboardSidebar(
            collapsed = TRUE,
            width = 300,
            tags$head(tags$style(HTML('
                   /* Sidebar overall padding */
                   .main-sidebar, .left-side {
                   padding-top: 20px;
                   padding: 20px; /* Adjust as needed for padding around content */
                   }'
           ))),
           p("text "), # just used for spacer, otherwise "Input Files" overlaps with header for some reason
           h3("Input files"),
           shinyDirButton("file_dir", "File directory", "Select a folder", style="width:200px"),
           verbatimTextOutput("fileDir", placeholder = TRUE),
           hr(style = "border-top: 1px solid #ccc; margin: 10px 0;"), # Add a styled horizontal rule
           
           p("Load the manifest file.
             (Default: VCF_Table_Viewer_CCDI_manifest.csv in the data directory. 
              Required columns: VCFFileName, StudyID, ParticipantID, SampleID, caller)"),
           fileInput("manifest", label = "Select a manifest file (.csv format):", 
                     accept = ".csv"),
           #verbatimTextOutput("stylesheet_file", placeholder = TRUE),
           hr(style = "border-top: 1px solid #ccc; margin: 10px 0;"), # Add a styled horizontal rule
           
           #p("Specify the directory with the post-pipeline filtered VCF files."),
           #shinyDirButton("vcf_dir", "Select the VCF file directory", "Select a folder", style="width:200px"), # 
           #verbatimTextOutput("vcfDir", placeholder = TRUE),
           #hr(style = "border-top: 1px solid #ccc; margin: 10px 0;"), # Add a styled horizontal rule
           
           #p("Specify the directory with the BAM files. Default: preprocessing/recalibrated in the sarek output directory."),
           #shinyDirButton("bam_dir", "Select the BAM file directory", "Select a folder", style="width:200px"),
           #verbatimTextOutput("bamDir", placeholder = TRUE),
           #hr(style = "border-top: 1px solid #ccc; margin: 10px 0;"), # Add a styled horizontal rule
           
           p("Load the lists of genes of interest to highlight."),
           fileInput("gene_list", label = "Select the gene lists (tab-delimited, with headers):")
           #         accept = c(".txt", ".tsv")#,
           

      ), # dashboardSidebar
                
      dashboardBody(
                # change menus to studyID, subjectID, sampleID, caller
                # filters are Population:
                #             Field: (default MAX_AF), value (default 0.01)
                # and Mutation:
                #             Field: (default significant), value (list)            
                fluidRow(column(3, selectInput("studyID", "Study",
                                               choices = NULL)),
                         column(3, selectInput("participantID", "Participant", 
                                               choices = NULL),
                                               selected = "PT_00G007DM"),
                         column(3, selectInput("sampleID", "Sample",
                                               choices = NULL,
                                               multiple = TRUE)),
                         column(3, selectInput("caller", "Caller",
                                               choices = NULL,
                                               multiple = TRUE))
                         ),
                fluidRow(column(12, selectInput("vcfFile", "VCF File",
                                                choices = NULL,
                                                width = "600px"))),
                
                fluidRow(column(12, align= "right", htmlOutput("numberOfVariants"))),

                fluidRow(column(12, 
                                
                        tabBox( title = "",
                                id = "main",
                                width = 12,
                                    tabPanel("README",
                                         htmlOutput("readme")
                                    ),
                                    tabPanel("Table",
                                        div(dataTableOutput("dataTable"))
                                    ),
                                    tabPanel("Legend",
                                        htmlOutput("legend")
                                    ),
                                    tabPanel("BAM Viewer",
                                        fluidRow(
                                                      column(3, selectInput("variantList", "Selected Variants",
                                                                            choices = c(""))),
                                                      column(2, actionButton("removeTracks", "Remove Tracks"),
                                                             align = "left",
                                                             style = "margin-top: 25px;"),
                                                      column(7, htmlOutput("bamfile"))
                                             ),
                                             igvShinyOutput('igvShiny_0')
                                    ),
                                    tabPanel("Plots",
                                             h5("Allele Frequency from Mutect2 VCFs"),
                                             fluidRow(
                                                      column(4, actionButton("createPlot", "Generate plot for selected rows")),
                                                      column(4, downloadButton("downloadPlot", "Download Plot"))),
                                             plotOutput("plot")
                                    )
                                             #htmlOutput("viewer"))
                        ) # tabBox
                )) # fluidRow(column(12,
      ) # dashboard body
) # ui

# server is where all calculations are done, tables are pre-rendered
server <- function(input, output, session) {
  
  # Get the VCF and BAM file directories
  roots=c(wd='.', vol='/Volumes', mnt='/mnt')
  
  #get_vcf_dir <- reactive({
  #  shinyDirChoose(input, 'vcf_dir', roots=roots)  
  #  output$vcfDir <- renderText(as.character(parseDirPath(roots=roots, input$vcf_dir)))
  #    
  #  vcfpath <- as.character(parseDirPath(roots=roots, input$vcf_dir))
  #  print(vcfpath)
  #  
  #  return(vcfpath)
  #})
  
  #get_bam_dir <- reactive({
  #  shinyDirChoose(input, 'bam_dir', roots=roots) 
  #  output$bamDir <- renderText(as.character(parseDirPath(roots=roots, input$bam_dir)))
  #    
  #  bampath <- as.character(parseDirPath(roots=roots, input$bam_dir))
  #  print(bampath)
  #    
  #  return(bampath)
  #})

  
  # sampleSheet <- reactive({
  #     req(global$sampleSheet_path)
  #     read.csv(global$sampleSheet_path)
  #     #sampleSheet <-ss |> arrange(patient)
  #     #return(sampleSheet)
  # })
  
  #observeEvent(input$samplesheet, {
  #  sampleSheet <- read.csv(input$samplesheet$datapath)
  #  sampleSheet <- sampleSheet |> arrange(patient)
  #  updateSelectInput(session, "subjectID", choices = unique(sampleSheet$patient))
  #})
  
  #observeEvent(input$vcf_dir, {
  #  global$vcfDir <- input$vcf_dir
  #})
  
  # First we read in the manifest.  Then we populate the studyID dropdown with the available 
  # studyIDs.  The user selects a participantID, then a studyID, then a caller, and can 
  # optionally select FileAccess (Open or Controlled - both will be shown if not selected)
  
  #EDITING this to change the default to be 
  df_manifest <- reactive({   
    if (is.null(input$manifest)) {
      manifest_file <- "sbgenomics/project-files/VCF_Table_Viewer_CCDI_manifest.csv"
    } else { 
      manifest_file <- input$manifest$datapath
    }
    read.csv(manifest_file, header=TRUE, sep=",")
  })
  
  ccdi_study <- reactive({
    req(df_manifest())
    unique(df_manifest() %>% pull(StudyID)) 
  })
  
  ccdi_participant <- reactive({
    req(df_manifest())
    req(input$studyID)
    unique(df_manifest() %>% filter(StudyID==input$studyID) %>% 
             pull(ParticipantID))
  })
  
  ccdi_sample <- reactive({
    req(df_manifest())
    req(input$studyID)
    req(input$participantID)

    unique(df_manifest() %>% filter(StudyID==input$studyID, ParticipantID==input$participantID) %>%
                             pull(SampleID))
  })
  
  ccdi_caller <- reactive({
    req(df_manifest())
    req(input$studyID)
    req(input$participantID)
    req(input$sampleID)
    
    unique(df_manifest() %>% filter(StudyID==input$studyID, 
                                    ParticipantID==input$participantID,
                                    SampleID==input$sampleID) %>%
                             pull(caller))
  })
  
  ccdi_vcf_options <- reactive({
    req(df_manifest())
    req(input$studyID)
    req(input$participantID)
    req(input$sampleID)
    req(input$caller)
    # caller
    # FileAccess
    
    unique(df_manifest() %>% filter(StudyID==input$studyID) %>%
                             filter(ParticipantID==input$participantID) %>%
                             filter(SampleID==input$sampleID) %>%
                             #filter(FileAccess==input$fileAccess) %>%
                             filter(caller==input$caller) %>%
                             pull(VCFFileName))
  })
  #Warning: There was 1 warning in `filter()`.
  # In argument: `caller == input$caller`.
  #Caused by warning in `caller == input$caller`:
  #  ! longer object length is not a multiple of shorter object length
  
  observe({
    updateSelectInput(session, "studyID", choices=ccdi_study())
  })
  
  observe({
    updateSelectInput(session, "participantID", choices=ccdi_participant())
  })
  
  observe({
    updateSelectInput(session, "sampleID", choices=ccdi_sample())
  })
  
  observe({
    updateSelectInput(session, "caller", choices=ccdi_caller(), selected = 'consensus')
  })
  
  observe({
    updateSelectInput(session, "vcfFile", choices=ccdi_vcf_options())
  })
  
  observeEvent(input$vcfFile, {
    req(input$vcfFile)
    
    
  })
  
  #observeEvent(input$bam_dir, {
  #})

  #-----------------------------------------------------------------------------
  #  generate variant dataframe
  #-----------------------------------------------------------------------------
  inputTable <- reactive({
    
    print("Getting input data frame")

    #vcfpath <- get_vcf_dir()
    #bampath <- get_bam_dir()
    #req(input$samplesheet)
    
    #sampleSheet <- get_samplesheet()
    #vcfpath <- get_vcf_dir()
    
    #print(paste("samplesheet:", global$sampleSheet))
    #print(paste("vcfpath:", global$vcfDir))
    
    # if (is.integer(input$vcf_dir)) { # no selection made
    #   vcfpath <- vcfDir_default
    #   output$vcfDir <- renderText(vcfDir_default)
    # } else {
    #   vcfpath <- get_vcf_dir()
    # }
    
    studyID <- input$studyID
    subjID <- input$subjectID
    caller <- input$caller
    sampleID <- input$sampleID
    
    #filterLevel <- input$filterLevel # region, population, mutation, driver, genes of interest
    
    #inputDir <- paste0(global$vcfDir, caller) # sep = "/"
    inputDir <- "sbgenomics/project-files/"
    print(inputDir)

    inFile <- "NULL"
    
    #filterName <- filterNames[filterLevel]
    infile_df <- df_manifest() |> filter(StudyID==input$studyID) |>
                                  filter(ParticipantID==input$participantID) |>
                                  filter(SampleID==input$sampleID) |>
                                  filter(caller==input$caller)
                                  #filter(FileAccess==input$fileAccess) |>
    
    # select infile from manifest df
    #infile_df <- manifest |> filter(StudyID == studyID, ParticipantID == subjID, 
    #                             SampleID == sampleID, Caller == caller) |> select(VCFFileName)
    inFile <- paste0(inputDir, input$vcfFile)
    
    # TODO: check VCF file size - if GB, print warning that only 100k(?) lines will be read in, 
    # suggest filtering VCF first
    get_vcf_file_size <- function(size) {
      x <- unlist(strsplit(size, " "))
      if (x[2] == "GB" | (x[2] == "MB" & as.numeric(x[1]) > 5 )) {
        return(1)
      } else {
        return(0)
      }
    }
    VCFFileSize <- infile_df$VCFFileSize[match(input$vcfFile, infile_df$VCFFileName)]
    print(VCFFileSize)
    
    if (get_vcf_file_size(VCFFileSize)) {
      shinyalert("Large VCF File", "This VCF file is large - only the first 100,000 variants will be loaded.", 
                 type = "warning")
    }
    
    print(paste("Reading vcf file: ", inFile))
    req(inFile)
    
    vcf <- read.vcfR(inFile, checkFile = TRUE)
    # Note: using nrows = 100000L creates a dataframe with 100,000 rows even if
    # there are fewer variants than that in the file.
    print(paste("vcf dimensions: ", dim(vcf)))
    
    # check that we have variants
    #  - getFIX returns a character vector if there is 1 variant, dataframe otherwise
    numVariants <- 0
    if (!is.vector(getFIX(vcf))) {
      if (dim(getFIX(vcf))[1] == 0) {
        output$numberOfVariants <- renderText({ paste("Number of variants in VCF: ", as.character(numVariants)) })
        return(NULL) # return null table instead of error message
      } else { numVariants <- dim(getFIX(vcf))[1] }
    } else {
      numVariants <- 1
    }
    # prints number of variants for given selections above the table
    output$numberOfVariants <- renderText({ paste("Number of variants in VCF: ", as.character(numVariants)) })

    
    # mutect2 FORMAT:
    # GT:AD:AF:DP:F1R2:F2R1:FAD:SB    0/1:124,17:0.065:141:47,0:27,2:100,12:51,73,0,17
    # 
    # haplotypecaller FORMAT:
    # GT:AD:DP:GQ:PL  0/1:53,41:94:99:1246,0,1560
    #
    # extract.gt(vcf, "DP", as.numeric=T) gets just the DP field
    # -> extract = F collects the rest of the genotype info except for the named field (e.g. "DP")
    # vcf@gt gets the whole gt field 
    
    print("getting fixed and info...")
    # create data frame from VCF (fixed + genotype + INFO)
    if (!is.vector(getFIX(vcf))) {
      fixed <- getFIX(vcf)
    } else {
      fixed <- t(data.frame(getFIX(vcf)))
    }
    info <- INFO2df(vcf)
    print(colnames(info))
    my.vcf.df <- cbind(as.data.frame(fixed), vcf@gt, info)
    print("...done")
    
    # split the INFO annotations into columns
    ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. 
    ##Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|
    ##HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|
    ##ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|
    ##CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|RefSeq|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|
    ##GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|HGVSg|AFR_AF|AMR_AF|
    ##EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|
    ##gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|
    ##PUBMED|CHECK_REF|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">    
    
    #cols <- gsub(pattern='\\n',replacement="",x=cols)
    #cols <- gsub(pattern='\\s',replacement="",x=cols)
    #newcols <- strsplit(cols, "\\|")
    #newcols <- unlist(newcols)
    
    # get columns from vcf file directly
    # ##INFO=<ID=CSQ
    print("Getting annotation columns...")
    ann <- grepl("ID=CSQ", vcf@meta)
    txt <- vcf@meta[ann]
    print(txt)
    
    cols <- unlist(strsplit(vcf@meta[ann], "Format: ")) # get rid of leading text
    cols2 <- gsub("\\\">","" , cols[2]) # get rid of trailing text
    newcols <- unlist(strsplit(cols2, "\\|"))
    # print(newcols)

    # print(colnames(my.vcf.df))
    # if multiple annotations, only take the 1st one (TODO: deal with multiple alleles)
    my.vcf.ANN.df <- my.vcf.df |> separate_wider_delim(CSQ, delim=",", names = c("CSQ"), too_many="drop")
    my.vcf.ANN.df <- my.vcf.ANN.df |> separate_wider_delim(CSQ, delim = "|", names = newcols,  
                                                           too_many = "debug", too_few = "debug", 
                                                           names_repair = "universal") # 
    # rename AF columns
    if (caller == "haplotypecaller") {
      my.vcf.ANN.df <- my.vcf.ANN.df |> rename("AF...11" = "AF", "AF...70" = "AF_TG") # requires AF to be in columns 11 and 69
      #write.table(my.vcf.ANN.df, file="my.vcf.ANN.df.txt", quote=F)
    }
    #print(dim(my.vcf.ANN.df))
    #Warning: Expected 93 pieces. Additional pieces discarded in 7 rows [1, 2, 3, 4, 5, 6, 7].

    #### filter based on population frequency
    # Default: MAX_AF
    # TODO: allow user to select field instead of MAX_AF
    # pop_freq_cutoff (default 0.01)
    # 
    pop_freq_cutoff <- 1
    my.vcf.ANN.df <- my.vcf.ANN.df |> mutate(across(matches(c("gnomAD", "_AF")), \(x) as.numeric(x) )) 
                                   # |> filter(gnomAD_AF <= pop_freq_cutoff)
    
    #### filter based on mutation severity
    significant_mutation_list = c("start_lost", "stop_lost", "stop_gained", "missense_variant",
                                  "frameshift_variant", "disruptive_inframe_insertion", "disruptive_inframe_deletion",
                                  "conservative_inframe_deletion", "conservative_inframe_insertion", 
                                  "splice_donor_variant", "splice_acceptor_variant", "splice_region_variant")
                                  
    less_significant_mutation_list = c("stop_retained_variant",
                                       "5_prime_UTR_premature_start_codon_gain_variant",
                                       "5_prime_UTR_variant",
                                       "3_prime_UTR_variant",
                                       "synonymous_variant")
    # user selects mutation_filter (default "significant")
    mutation_filter = "significant"
    if (mutation_filter == "less_significant") {
      mutation_list = c(significant_mutation_list, less_significant_mutation_list)
    } else {
      mutation_list = significant_mutation_list
    }
    my.vcf.ANN.df <- my.vcf.ANN.df |> filter(Consequence %in% mutation_list)
    
    #### create columns flagging genes of interest
    for (i in colnames(gene_lists)) {
      my.vcf.ANN.df[[i]] <- ifelse(my.vcf.ANN.df$SYMBOL %in% gene_lists[[i]] & 
                                   my.vcf.ANN.df$SYMBOL != '', 'Y', 'N')
    }
    
    #### create an index from chr,pos,ref,alt
    my.vcf.ANN.df$index <- paste(my.vcf.ANN.df$CHROM, my.vcf.ANN.df$POS, my.vcf.ANN.df$REF, my.vcf.ANN.df$ALT, sep='.')
    
    #### create aa mutation e.g. I255T
    #  Protein_position 70/393
    #  Amino_acids I/T
    
    my.vcf.ANN.df <- my.vcf.ANN.df |>
      mutate(prot_pos = str_extract(Protein_position, "(^\\d+)/", group=1)) |>
      mutate(AA1 = str_extract(Amino_acids, "(^\\w+)/", group=1)) |> 
      mutate(AA2 = str_extract(Amino_acids, "/(\\w+$)", group=1)) |>
      mutate(AA_mut = case_when(Consequence == "missense_variant" | Consequence == "frameshift_variant" 
                                ~ paste0(AA1, prot_pos, AA2))) |>
      select(!prot_pos, !AA1, !AA2)
      
    #### Extract 1 value from predictors that give values for each transcript:
    # PROVEAN_pred, PolyPhen2_HDIV_pred, PolyPhen2_HVAR_pred, REVEL_score
    # my.vcf.ANN.df <- my.vcf.ANN.df |> 
    #                  mutate(MT_pred = str_extract(MutationTaster_pred, "\\w"), .keep="unused", .after="MetaSVM_pred") |>
    #                  mutate(FATHMM_pred = str_extract(FATHMM_pred, "\\w"), .keep="unused", .after="DANN_score") |>
    #                  mutate(PROVEAN_pred = str_extract(PROVEAN_pred, "\\w"), .keep="unused", .after="MT_pred") |>
    #                  mutate(PP2_HDIV_pred = str_extract(Polyphen2_HDIV_pred, "\\w"), .keep="unused", .after="PROVEAN_pred") |>
    #                  mutate(PP2_HVAR_pred = str_extract(Polyphen2_HVAR_pred, "\\w"), .keep="unused", .after="PP2_HDIV_pred") |>
    #                  mutate(REVEL_score = str_extract(REVEL_score, "\\d*\\.?\\d+"), .keep="unused")
    # 
    #### order the columns logically  "Leudrive", "ACMG", 
    my.vcf.ANN.df <- my.vcf.ANN.df |> 
      relocate(c("Allele", "FILTER", "SYMBOL", "AA_mut", "IMPACT", "Consequence",  "Existing_variation"), .after=QUAL) |>
      relocate(c("CLIN_SIG", "SIFT", "PolyPhen")) #|> # .after=SOMATIC
      #relocate(c("REVEL_score"), .after=DANN_score)
    my.vcf.ANN.df <- my.vcf.ANN.df |> relocate(index)
    
    #### make columns numeric  
    # my.vcf.ANN.df <- my.vcf.ANN.df |> mutate(across(c('DANN_score', 'CADD_phred', 'REVEL_score'), \(x) as.numeric(x))) |>  
    #                                    mutate(across(c('DANN_score', 'CADD_phred', 'REVEL_score'), \(x) round(x, 3)))   |>
    #                                    mutate_if(is.numeric, ~replace(., is.na(.), 0))
    #print(head(my.vcf.ANN.df))
    
    # Deterrmine which columns to hide upon initial table display
    hide_columns <- function(df, vcf) {
      # get column indices for columns to begin the display hidden
      # TODO: create groupings of columns to hide/show with a click
      # TODO: make this easier to modify for default
      
      # select columns to show:
      show_cols_text <- "index,CHROM,POS,REF,ALT,FILTER,DP,GERMQ,
                         POPAF,Allele,Consequence,IMPACT,SYMBOL,AA_mut,Leudrive,ACMG,Gene,
                         Feature_type,Feature,BIOTYPE,EXON,INTRON,cDNA_position,CDS_position,
                         SWISSPROT,SIFT,PolyPhen,AF,gnomADe_AF,MAX_AF,FREQS,
                         CLIN_SIG,SOMATIC,CADD_phred,DANN_score,FATHMM_pred,LRT_pred,MetaSVM_pred,
                         MT_pred,PROVEAN_pred,PP2_HDIV_pred,PP2_HVAR_pred,PrimateAI_pred,
                         REVEL_score,gnomAD_AF"
      
      # # group columns (not used):
      # show_cols_fixed <- c("index", "CHROM", "POS", "REF", "ALT")
      # show_cols_gt <- c()
      # 
      # # need a named list with vectors for each caller
      # show_cols_info <- list(haplotypecaller = list(),
      #                        mutect2 = list(),
      #                        lancet = list(),
      #                        manta = list(),
      #                        svaba = list(),
      #                        strelka2 = list(),
      #                        deepvariant = list())
      # # population group
      # show_cols_pop <- c("AF", "POPAF", "gnomADe_AF", "gnomAD_AF", "MAX_AF")
      # # gene group
      # show_cols_gene <- c("Allele,Consequence,IMPACT,SYMBOL,AA_mut,Leudrive,ACMG,Gene,Feature_type,
      #                      Feature,BIOTYPE,EXON,INTRON,cDNA_position,CDS_position,Protein_position,
      #                      Amino_acids,Existing_variation,VARIANT_CLASS,APPRIS,CCDS,ENSP,SWISSPROT, 
      #                      miRNA")
      # # impact group
      # show_cols_impact <- c("SIFT,PolyPhen,CLIN_SIG,SOMATIC,CADD_phred,DANN_score,FATHMM_pred,
      #                       LRT_pred,MetaSVM_pred,MT_pred,PROVEAN_pred,PP2_HDIV_pred,PP2_HVAR_pred,
      #                       PrimateAI_pred,REVEL_score")
      
      # replace line returns and spaces
      show_cols_text <- gsub(pattern='\\n',replacement="",x=show_cols_text)
      show_cols_text <- gsub(pattern='\\s',replacement="",x=show_cols_text)
      TEMP <- scan(text=show_cols_text, what="", sep=",", quiet=TRUE) # split into separate items
      show_cols = dput(TEMP, file="tmp.txt") # add quotes
      
      # get the sampleIDs 
      gt <- vcf@gt
      #print(colnames(gt))
      #[1] "FORMAT"      "BS_A1DV9T7G" "BS_GJWAV3E5"
      show_cols <- append(show_cols, colnames(vcf@gt))
      
      # get columns with genes of interest
      show_cols <- append(show_cols, colnames(gene_lists))
      
      hide_cols <- which(!(colnames(df) %in% show_cols)) - 1
      return(hide_cols)
    }
    
    hide_cols <- hide_columns(my.vcf.ANN.df, vcf) # list of indices to hide initially
    
    return(list(df = my.vcf.ANN.df, hide_cols = hide_cols))
  })
  
  
  #-----------------------------------------------------------------------------
  #  render data table
  #-----------------------------------------------------------------------------
  
  output$dataTable <- renderDT({
    
    inputTable <- inputTable()
    df <- inputTable$df
    hide_cols <- inputTable$hide_cols

    if (is.null(df)) {
        validate("No variants found, select another filter")
        #output$numberOfVariants <- renderText({ "Number of variants: 0" })
    } else {
      # Prints the number of variants above the table
      #output$numberOfVariants <- renderText({ paste("Number of variants: ", dim(df)[1]) })
      print("df has variants")
    }
    
    #hide_cols <- hide_columns(df, vcf) # list of indices to hide initially

    deleterious_columns <- c('MT_pred', 'PROVEAN_pred', 'PrimateAI_pred','FATHMM_pred',
                             'LRT_pred','MetaSVM_pred', 'PP2_HDIV_pred', 'PP2_HVAR_pred')
    color_gradient_columns <- c("CADD_phred", "DANN_score")
    
    # get index of FILTER column
    FILTER_index <- grep("FILTER", names(df))
    filter_list <- rep(list(NULL),FILTER_index)
    filter_list[[FILTER_index]] <- list(search="PASS") # start with only variants that PASS  
    
    dt <- DT::datatable(
      df, 
      rownames=FALSE, 
      extensions = c('FixedColumns','FixedHeader','Buttons'), 
      options = list(
        dom = 'Blfrtip',
        fixedHeader=TRUE,
        fixedColumns=TRUE,
        pageLength = 100,
        lengthMenu = list(c(50, 100, -1), c('50','100','All')),
        autoWidth = FALSE,
        scrollX = TRUE,
        columnDefs = list(list(visible = FALSE, targets = hide_cols)
                          #list(targets = 0, width = '50px') # c(0,11) doesn't work for some reason
        ),
        buttons = list('colvis',
                    list(
                      extend = "excel",
                      text = 'Excel',
                      exportOptions = list(rows = '.selected') # only export selected rows
                    )),
        searchCols = list(filter_list)  # initialize filters on each column
      ), # options
      class = "display nowrap compact", # style
      filter = "top" # location of column filters
    # format columns
    ) |>
    formatStyle('IMPACT',
                backgroundColor = styleEqual(c("HIGH", "MODERATE"), c('#FA5F55', 'yellow'))
    ) |> 
    formatStyle(colnames(gene_lists),
                backgroundColor = styleEqual(c("Y"), c('lightgreen'))
    ) |>
    formatStyle('SYMBOL', colnames(gene_lists), fontWeight = styleEqual("Y", "bold"))  
    
    for (s in color_gradient_columns) {
      if (s %in% colnames(df)) { dt <- dt |> color_gradient(s) } # "gnomAD_AF"
    }
    
    for (p in deleterious_columns) {
      if (p %in% colnames(df)) {
        dt <- dt |> formatStyle(p, backgroundColor = styleEqual(c("D"), c('#FA5F55')))  
      }
    }
    
    if ("REVEL_score" %in% colnames(df)) {
      dt <- dt |> formatStyle('REVEL_score',
                              backgroundColor  = styleInterval(c(0.5), c('white','#FA5F55')))
    }
    
    dt
  }) # renderDT
  
  # explanation of vcf labels in table
  output$legend <- renderUI({
    tags$iframe(
      seamless="seamless",
      src="tmpuser/vcf_field_descriptions.html",
      width=800, 
      height=800)
  })
  
  # README file
  output$readme <- renderUI({
    tags$iframe(
      seamless="seamless",
      src="tmpuser/Instructions.html",
      width=800,
      height=800)
  })
  
  
  #-----------------------------------------------------------------------------#
  # genome viewer
  #-----------------------------------------------------------------------------#

  observeEvent(input$dataTable_rows_selected, {
    x <- inputTable()[input$dataTable_rows_selected, ]
    
    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)
    
    # Update the selectInput menu
    updateSelectInput(session, "variantList", "Selected Variants", choices = x["index"], selected = c(""))
    
  })
  
  observeEvent(input$variantList, {
    
    req(input$variantList)
    
    if (Sys.info()['nodename'] == "NCI-02295810-ML") {
      # Posit::Connect server
      # -> The Biowulf file system is mounted at /mnt/BW-Data/ on appshare-dev
      bampath <- "/mnt/BW-Data/recalibrated/"
    }
    
    if (Sys.info()['nodename'] == "NCI-02295810-ML") {
      global$bamDir <- "/Volumes/sierk/runx/nf/sarek/sarek_april2025/preprocessing/recalibrated/"
    } 
    output$bamDir <- renderText(global$bamDir)
    
    x <- inputTable()
    variant <- x[x$index == input$variantList, ] # input$variantList
    chrom_pos <- paste0(variant$CHROM, ":", variant$POS)
    showGenomicRegion(session, id="igvShiny_0", chrom_pos) # chr21:10,397,614-10,423,341
    
    samples <- variant |> select(starts_with("FPD_")) # FPD_0028_FPD_0028_SK211E
    #print(samples)
    sample_names <- colnames(samples)
    
    for (s in sample_names) {
      print(s)
      sampleID <- substr(s,10,24)
      print(paste0("sampleID: ", sampleID))
      
      subjID <- input$subjectID
      bamFile <- paste0(global$bamDir, sampleID, "/", sampleID, ".recal.bam")
      #bamFile <- "/Volumes/sierk/runx/nf/sarek_august2024/preprocessing/recalibrated/FPD_0271_BM221E/FPD_0271_BM221E.recal.bam"
      if (file.exists(bamFile)) {
        output$bamfile <- renderText({paste("loading bam file:", bamFile)})
        x <- readGAlignments(bamFile, param = Rsamtools::ScanBamParam(what="seq", which=GRanges(chrom_pos)))
        loadBamTrackFromLocalData(session, id="igvShiny_0", trackName=sampleID, data=x)
      } else {
        output$bamfile <- renderText({ paste("<font color=\"#FF0000\"><b>", "bam file missing: ", "</b></font>", bamFile) })
      }
    }

    #bamFileDisplay <- paste0(bamDir, "\n", sampleID, "/", sampleID, ".recal.bam")
    #output$bamfile <- renderText(bamFileDisplay)

  }, ignoreInit = TRUE) # don't display until user clicks on dropdown menu
  
  # from igvShinyDemo, potentially useful:
  observeEvent(input$removeTracks, {
    printf("---- removeUserTracks")
    removeUserAddedTracks(session, id="igvShiny_0")
  })
  
  observeEvent(input$igvReady, {
    printf("--- igvReady")
    containerID <- input$igvReady
    printf("igv ready, %s", containerID)
    # loadBedTrack(session, id=containerID, trackName="bed5 loaded on ready", tbl=tbl.bed5, color="red");
  })
  
  observeEvent(input$trackClick, {
    printf("--- trackclick event")
    x <- input$trackClick
    print(x)
  })
  
  observeEvent(input[["igv-trackClick"]], {
    printf("--- igv-trackClick event")
    x <- input[["igv-trackClick"]]
    print(x)
    attribute.name.positions <- grep("name", names(x))
    attribute.value.positions <- grep("value", names(x))
    attribute.names <- as.character(x)[attribute.name.positions]
    attribute.values <- as.character(x)[attribute.value.positions]
    tbl <- data.frame(name=attribute.names,
                      value=attribute.values,
                      stringsAsFactors=FALSE)
    dialogContent <- renderTable(tbl)
    html <- HTML(dialogContent())
    showModal(modalDialog(html, easyClose=TRUE))
  })
  
  observeEvent(input$getChromLocButton, {
    # printf("--- getChromLoc event")
    # sends message to igv.js in browser; currentGenomicRegion.<id> event sent back
    # see below for how that can be captured and displayed
    getGenomicRegion(session, id="igvShiny_0")
  })
  
  observeEvent(input$clearChromLocButton, {
    output$chromLocDisplay <- renderText({" "})
  })
  
  observeEvent(input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]], {
    newLoc <- input[[sprintf("currentGenomicRegion.%s", "igvShiny_0")]]
    #printf("new chromLocString: %s", newLoc)
    output$chromLocDisplay <- renderText({newLoc})
  })
  
  output$igvShiny_0 <- renderIgvShiny({
    cat("--- starting renderIgvShiny\n");
    genomeOptions <- parseAndValidateGenomeSpec(genomeName="hg38")
    x <- igvShiny(genomeOptions,
                  displayMode="SQUISHED",
                  tracks=list()
    )
    cat("--- ending renderIgvShiny\n");
    return(x)
  })
  
  #-----------------------------------------------------------------------------#
  
  #-----------------------------------------------------------------------------#
  # plot of allele frequencies
  #-----------------------------------------------------------------------------#
  observeEvent(input$createPlot, {

    subjID <- input$participantID #changing subjectID to participantID
    print(paste("Generating plot for", subjID, "..."))

    x <- inputTable()[input$dataTable_rows_selected, ]
    
    # Genotype formatting
    # haplotypecaller: 0/1:53,41:94:99:1178,0,1559
    #                  GT:AD:DP:GQ:PL
    # mutect2: GT:AD:AF:DP:F1R2:F2R1:FAD:SB
    #          0/1:240,8:0.005226:248:58,0:100,0:208,4:135,105,8,0
    #          SB “Per-sample component statistics which comprise the Fisher’s Exact Test to detect strand bias.”
    
    # need gene symbol + aa mutation if available, select out GT fields
    samples <- x |> mutate(index = case_when(!is.na(AA_mut) ~ paste0(SYMBOL, "(", AA_mut, ")"), 
                                              .default = index)) |> 
                     select(index, starts_with("FPD_")) 
                     
    print(samples)

    getAF <- ~as.numeric(unlist(strsplit(.x, ":"))[3])
    samples <- samples |> rowwise() |> mutate(across(starts_with("FPD_"), getAF)) |>
                       rename_with(~str_remove(., "^FPD_[\\d]{4}_")) |>
                       select(contains("SK"), sort(colnames(.)))
    #print(samples)
      
    factor(substring(x, 1, 2)) # orders the sample names
    
    samples_melt <- melt(as.data.table(samples), id = "index")
    #print(samples_melt)
    
    generatePlot <- function() {
      ggplot(data = samples_melt, aes(x = variable, y = value, color = index, group = index)) +
        geom_point(size = 2) + 
        geom_line() + 
        labs(title = paste(subjID, "Mutect2"), x = "Samples", y = "Allele Frequency", color = "Variants") +
        theme(text = element_text(size = 16),
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.text = element_text(size=10))
    }
      
    if (input$caller == "mutect2") {
        
        output$plot <- renderPlot({
          generatePlot()  
        }) # renderPlot
        
        #observeEvent(input$downloadPlot, {
        output$downloadPlot <- downloadHandler(
          filename = "allele_freqs.pdf",
          content = function(file) {
            #pdf(file)
            #generatePlot()
            #dev.off()
            ggsave(file, plot=generatePlot(), dpi = 300, 
                   width = 8.5, height = 5.5, units = "in", device="pdf")
            #filename = function(){paste("input$plot3",'.png',sep='')},
            #content = function(file){
            #  ggsave(file,plot=data$plot)
          }
        )
    } # if input$caller == mutect2
    
  }) # observeEvent
  

      
} # server


# run the app
shinyApp(ui, server) # launch.browser = TRUE, options = list(width = 1600)

# Workaround for reading in cram files:
# (from https://github.com/gladkia/igvShiny/issues/102)
# 
# Install the servr package if you don't have it already:
# 
# Start a local web server: In your R console, run the following command, replacing "path/to/your/files" with the actual path to the directory containing your CRAM and CRAI files:
#   
#   servr::httd(dir = "path/to/your/files")
# This will start a local web server, usually at http://127.0.0.1:4321.
# 
# Note: inside Docker, do this (https://github.com/yihui/servr/issues/49):
# 
#   deamon_id <- servr::httd(port = 8001, daemon = TRUE, host = '0.0.0.0')
#   cramURL <- "http://0.0.0.0:8001/your_file.cram"
#
# Load the CRAM track in igvShiny: In your Shiny app, you can now use loadCramTrackFromURL, 
# providing the local URLs for your CRAM and CRAI files.
# 
# # In your Shiny server function
# ...
# session <- shiny::getDefaultReactiveDomain()
# 
# trackName <- "Local CRAM"
# cramURL <- "http://127.0.0.1:4321/your_file.cram"
# indexURL <- "http://127.0.0.1:4321/your_file.cram.crai"
# 
# loadCramTrackFromURL(session,
#                      trackName=trackName,
#                      trackColor="blue",
#                      dataURL=cramURL,
#                      indexURL=indexURL)
# ...
# Remember to replace "your_file.cram" with the name of your CRAM file.

# all_cols <- c(CHROM,POS,ID,REF,ALT,QUAL,FILTER,AS_FilterStatus,AS_SB_TABLE,AS_UNIQ_ALT_READ_COUNT,CONTQ,DP,ECNT,GERMQ,
#                MBQ,MFRL,MMQ,MPOS,NALOD,NCount,NLOD,OCM,PON,POPAF,ROQ,RPA,RU,SEQQ,STR,STRANDQ,STRQ,TLOD,LOF,NMD,
#                Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,
#                Codons,Existing_variation,DISTANCE,STRAND,FLAGS,VARIANT_CLASS,SYMBOL_SOURCE,HGNC_ID,CANONICAL,MANE_SELECT,MANE_PLUS_CLINICAL,TSL,APPRIS,CCDS,ENSP,
#                SWISSPROT,TREMBL,UNIPARC,UNIPROT_ISOFORM,GENE_PHENO,SIFT,PolyPhen,DOMAINS,miRNA,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF,gnomADe_AF,gnomADe_AFR_AF,
#                gnomADe_AMR_AF,gnomADe_ASJ_AF,gnomADe_EAS_AF,gnomADe_FIN_AF,gnomADe_NFE_AF,gnomADe_OTH_AF,gnomADe_SAS_AF,gnomADg_AF,gnomADg_AFR_AF,gnomADg_AMI_AF,
#                gnomADg_AMR_AF,gnomADg_ASJ_AF,gnomADg_EAS_AF,gnomADg_FIN_AF,gnomADg_MID_AF,gnomADg_NFE_AF,gnomADg_OTH_AF,gnomADg_SAS_AF,MAX_AF,MAX_AF_POPS,FREQS,
#                CLIN_SIG,SOMATIC,PHENO,PUBMED,MOTIF_NAME,MOTIF_POS,HIGH_INF_POS,MOTIF_SCORE_CHANGE,TRANSCRIPTION_FACTORS,CADD_phred,DANN_score,ExAC,FATHMM_pred,
#                Interpro_domain,LRT_pred,MetaSVM_pred,MutationTaster_pred,PROVEAN_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,PrimateAI_pred,REVEL_score,SIFT_pred,dbSNP)

# You need to add "l" (small letter "L") to dom, that makes Blfrtip:
# B - Buttons
# l - Length changing input control
# f - Filtering input
# r - pRocessing display element
# t - Table
# i - Table information summary
# p - Pagination control
#

