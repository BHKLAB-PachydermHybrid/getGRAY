###############################################
#########EXTRACTION OF SELECTION DATA #########
###############################################

file<- read.csv("/Users/anthonymammoliti/Desktop/selection.txt", stringsAsFactors = FALSE, header = FALSE, sep="\t")

Salmon2 <- grep("Salmon: 11.2", file$V1)
Salmon1 <- grep("Salmon: 11.3", file$V1)
Kallisto2 <- grep("Kallisto: 0.43.1", file$V1)
Kallisto1 <- grep("Kallisto: 0.44.0", file$V1)
hg19 <- grep("GRCh37", file$V1)
hg38 <- grep("GRCh38", file$V1)
transcript1_19 <- grep("Gencode v23lift37 Transcriptome", file$V1)
transcript2_19 <- grep("Ensembl GRCh37 v67 Transcriptome", file$V1)

transcript1_38 <- grep("Gencode v23 Transcriptome", file$V1)
transcript2_38 <- grep("Ensembl GRCh38 v89 Transcriptome", file$V1)
sens_2013 <- grep("DEC-2013", file$V1)
sens_2017 <- grep("NOV-2017", file$V1)


####### sensitivity #######

if (length(sens_2013 > 0)){
  
  param_sens <- "sensitivity='2013'"
  
} else if (length(sens_2017 > 0)) {
  
  param_sens <- "sensitivity='2017'"
  
} else {
  
  param_sens <- 0
}



########## hg19 ###########

if (length(hg19 > 0)){
  
  if (length(transcript1_19 > 0)) {
    
    param_transcript <- "/Users/anthonymammoliti/Desktop/Gencode.v23.annotation.RData"
    
    if (length(Salmon2 > 0)) {
      
      Salmon2 <- "GRAY_hg19_Gencode_Salmon-0-11-2Pipeline"
      Salmon2Pipeline <- paste("/pfs/",Salmon2, sep="")
      param_s2 <- "rnaseq_tool='salmon2'"
      
    } else {
      
      param_s2 <- 0
      Salmon2Pipeline <- 0
      
    }
    
    
    if (length(Salmon1 > 0)) {
      
      Salmon1 <- "GRAY_hg19_Gencode_Salmon-0-11-3Pipeline"
      Salmon1Pipeline <- paste("/pfs/",Salmon1, sep="")
      param_s1 <- "rnaseq_tool='salmon1'"
      
    }else {
      
      param_s1 <- 0
      Salmon1Pipeline <- 0
      
    }
    
    
    if (length(Kallisto1 > 0)) { 
      
      Kallisto1 <- "kallisto"
      Kallisto1Pipeline <- paste("/Users/anthonymammoliti/Desktop/",Kallisto1, sep="")
      param_k1 <- "rnaseq_tool='kallisto1'"
      
    }else {
      
      param_k1 <- 0
      Kallisto1Pipeline <- 0
    }
    
    if (length(Kallisto2 > 0)) { 
      
      Kallisto2 <- "GRAY_hg19_Gencode_Kallisto-0-43-1Pipeline"
      Kallisto2Pipeline <- paste("/pfs/",Kallisto2, sep="")
      param_k2 <- "rnaseq_tool='kallisto2'"
      
    }else {
      
      param_k2 <- 0
      Kallisto2Pipeline <- 0
    }
    
    if(length(Kallisto1>0) && length(Kallisto2>0)){
      
      param_k1 <- "rnaseq_tool='kallistocombined'"
      param_k2 <- 0
      
    }
    
    if(length(Salmon1>0) && length(Salmon2>0)){
      
      param_s1 <- "rnaseq_tool='salmoncombined'"
      param_s2 <- 0
      
    }
    
    if ((length(Kallisto1 > 0) | length(Kallisto2 > 0)) && (length(Salmon1 > 0) | length(Salmon2 > 0))) {
      
      param_k1 <- "rnaseq_tool='kallistosalmon'"
      param_k2 <- 0
      param_s1 <- 0
      param_s2 <- 0
      
    }
    
    
  }
  
  if (length(transcript2_19 > 0)) {
    
    param_transcript <- "/pfs/getGRAY/annotation/Ensembl_hg19_annotation.RData"
    
    if (length(Salmon2 > 0)) {
      
      Salmon2 <- "GRAY_hg19_Ensembl_Salmon-0-11-2Pipeline"
      Salmon2Pipeline <- paste("/pfs/",Salmon2, sep="")
      param_s2 <- "rnaseq_tool='salmon2'"
      
    } else {
      
      param_s2 <- 0
      Salmon2Pipeline <- 0
      
    }
    
    
    if (length(Salmon1 > 0)) {
      
      Salmon1 <- "GRAY_hg19_Ensembl_Salmon-0-11-3Pipeline"
      Salmon1Pipeline <- paste("/pfs/",Salmon1, sep="")
      param_s1 <- "rnaseq_tool='salmon1'"
      
    }else {
      
      param_s1 <- 0
      Salmon1Pipeline <- 0
      
    }
    
    
    if (length(Kallisto1 > 0)) { 
      Kallisto1 <- "GRAY_hg19_Ensembl_Kallisto-0-44-0Pipeline"
      Kallisto1Pipeline <- paste("/pfs/",Kallisto1, sep="")
      param_k1 <- "rnaseq_tool='kallisto1'"
      
    }else {
      
      param_k1 <- 0
      Kallisto1Pipeline <- 0
    }
    
    if (length(Kallisto2 > 0)) { 
      Kallisto2 <- "GRAY_hg19_Ensembl_Kallisto-0-43-1Pipeline"
      Kallisto2Pipeline <- paste("/pfs/",Kallisto2, sep="")
      param_k2 <- "rnaseq_tool='kallisto2'"
      
    }else {
      
      param_k2 <- 0
      Kallisto2Pipeline <- 0
    }
    
    if(length(Kallisto1>0) && length(Kallisto2>0)){
      
      param_k1 <- "rnaseq_tool='kallistocombined'"
      param_k2 <- 0
      
    }
    
    if(length(Salmon1>0) && length(Salmon2>0)){
      
      param_s1 <- "rnaseq_tool='salmoncombined'"
      param_s2 <- 0
      
    }
    
    if ((length(Kallisto1 > 0) | length(Kallisto2 > 0)) && (length(Salmon1 > 0) | length(Salmon2 > 0))) {
      
      param_k1 <- "rnaseq_tool='kallistosalmon'"
      param_k2 <- 0
      param_s1 <- 0
      param_s2 <- 0
      
    }
    
  }
  
  
}


####### hg38 ##############

if (length(hg38 > 0)){
  
  if (length(transcript1_38 > 0)) {
    
    param_transcript <- "/pfs/getGRAY/annotation/Gencode_hg38_annotation.RData"
    
    if (length(Salmon2 > 0)) {
      
      Salmon2 <- "GRAY_hg38_Gencode_Salmon-0-11-2Pipeline"
      Salmon2Pipeline <- paste("/pfs/",Salmon2, sep="")
      param_s2 <- "rnaseq_tool='salmon2'"
      
    } else {
      
      param_s2 <- 0
      Salmon2Pipeline <- 0
      
    }
    
    
    if (length(Salmon1 > 0)) {
      
      Salmon1 <- "GRAY_hg38_Gencode_Salmon-0-11-3Pipeline"
      Salmon1Pipeline <- paste("/pfs/",Salmon1, sep="")
      param_s1 <- "rnaseq_tool='salmon1'"
      
    }else {
      
      param_s1 <- 0
      Salmon1Pipeline <- 0
      
    }
    
    
    if (length(Kallisto1 > 0)) { 
      
      Kallisto1 <- "GRAY_hg38_Gencode_Kallisto-0-44-0Pipeline"
      Kallisto1Pipeline <- paste("/pfs/",Kallisto1, sep="")
      param_k1 <- "rnaseq_tool='kallisto1'"
      
    }else {
      
      param_k1 <- 0
      Kallisto1Pipeline <- 0
    }
    
    if (length(Kallisto2 > 0)) { 
      
      Kallisto2 <- "GRAY_hg38_Gencode_Kallisto-0-43-1Pipeline"
      Kallisto2Pipeline <- paste("/pfs/",Kallisto2, sep="")
      param_k2 <- "rnaseq_tool='kallisto2'"
      
    }else {
      
      param_k2 <- 0
      Kallisto2Pipeline <- 0
    }
    
    if(length(Kallisto1>0) && length(Kallisto2>0)){
      
      param_k1 <- "rnaseq_tool='kallistocombined'"
      param_k2 <- 0
      
    }
    
    if(length(Salmon1>0) && length(Salmon2>0)){
      
      param_s1 <- "rnaseq_tool='salmoncombined'"
      param_s2 <- 0
      
    }
    
    if ((length(Kallisto1 > 0) | length(Kallisto2 > 0)) && (length(Salmon1 > 0) | length(Salmon2 > 0))) {
      
      param_k1 <- "rnaseq_tool='kallistosalmon'"
      param_k2 <- 0
      param_s1 <- 0
      param_s2 <- 0
      
    }
    
    
  }
  
  if (length(transcript2_38 > 0)) {
    
    param_transcript <- "/pfs/getGRAY/annotation/Ensembl_hg38_annotation.RData"
    
    if (length(Salmon2 > 0)) {
      
      Salmon2 <- "GRAY_hg38_Ensembl_Salmon-0-11-2Pipeline"
      Salmon2Pipeline <- paste("/pfs/",Salmon2, sep="")
      param_s2 <- "rnaseq_tool='salmon2'"
      
    } else {
      
      param_s2 <- 0
      Salmon2Pipeline <- 0
      
    }
    
    
    if (length(Salmon1 > 0)) {
      
      Salmon1 <- "GRAY_hg38_Ensembl_Salmon-0-11-3Pipeline"
      Salmon1Pipeline <- paste("/pfs/",Salmon1, sep="")
      param_s1 <- "rnaseq_tool='salmon1'"
      
    }else {
      
      param_s1 <- 0
      Salmon1Pipeline <- 0
      
    }
    
    
    if (length(Kallisto1 > 0)) { 
      Kallisto1 <- "GRAY_hg38_Ensembl_Kallisto-0-44-0Pipeline"
      Kallisto1Pipeline <- paste("/pfs/",Kallisto1, sep="")
      param_k1 <- "rnaseq_tool='kallisto1'"
      
    }else {
      
      param_k1 <- 0
      Kallisto1Pipeline <- 0
    }
    
    if (length(Kallisto2 > 0)) { 
      Kallisto2 <- "GRAY_hg38_Ensembl_Kallisto-0-43-1Pipeline"
      Kallisto2Pipeline <- paste("/pfs/",Kallisto2, sep="")
      param_k2 <- "rnaseq_tool='kallisto2'"
      
    }else {
      
      param_k2 <- 0
      Kallisto2Pipeline <- 0
    }
    
    if(length(Kallisto1>0) && length(Kallisto2>0)){
      
      param_k1 <- "rnaseq_tool='kallistocombined'"
      param_k2 <- 0
      
    }
    
    if(length(Salmon1>0) && length(Salmon2>0)){
      
      param_s1 <- "rnaseq_tool='salmoncombined'"
      param_s2 <- 0
      
    }
    
    if ((length(Kallisto1 > 0) | length(Kallisto2 > 0)) && (length(Salmon1 > 0) | length(Salmon2 > 0))) {
      
      param_k1 <- "rnaseq_tool='kallistosalmon'"
      param_k2 <- 0
      param_s1 <- 0
      param_s2 <- 0
      
    }
    
  }
  
  
}



initialoption <- "verbose=FALSE, nthread=1,isoforms=FALSE"
finalparam <- paste(initialoption,param_s2,param_s1,param_k2,param_k1, param_sens, sep=",")
finalparam <- gsub(",0", "", finalparam)


getGRAYP <-
  function (
    verbose=FALSE,
    nthread=1,
    isoforms=FALSE,
    rnaseq_tool=c("kallisto1","kallisto2","salmon1","salmon2", "kallistocombined", "salmoncombined", "kallistosalmon"),sensitivity=c("2017","2013")){
    z <- list()
    options(stringsAsFactors=FALSE)
    badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    
    matchToIDTableCELL <- function(ids,tbl, column) {
      sapply(ids, function(x) {
        myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
        if(length(myx) > 1){
          stop("Something went wrong in curating cell ids")
        }
        return(tbl[myx, "unique.cellid"])
      })
    }
    
    cell_all <- read.csv(file = "/Users/anthonymammoliti/Desktop/GRAYLOCAL/annotation/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
    print(cell_all)
    curationCell <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
    curationTissue <- cell_all[which(!is.na(cell_all[ , "GRAY.cellid"])),]
    curationCell <- curationCell[ , c("unique.cellid", "GRAY.cellid")]
    curationTissue <- curationTissue[ , c("unique.tissueid", "GRAY.tissueid")]
    
    rownames(curationTissue) <- curationCell[ , "unique.cellid"]
    rownames(curationCell) <- curationCell[ , "unique.cellid"]
    
    drug_all <- read.csv(file="/Users/anthonymammoliti/Desktop/GRAYLOCAL/annotation/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
    curationDrug <- drug_all[which(!is.na(drug_all[ , "GRAY.drugid"])),]
    curationDrug <- curationDrug[ , c("unique.drugid", "GRAY.drugid")]
    rownames(curationDrug) <- curationDrug[ , "unique.drugid"]
    
    rnaseq.sampleinfo <- read.csv(file="/Users/anthonymammoliti/Desktop/GRAYLOCAL/exp/JRGraySRRMapping.csv", stringsAsFactors=FALSE, row.names=1)
    rnaseq.sampleinfo[ , "cellid"] <- rownames(curationCell)[match(rnaseq.sampleinfo[ , "cellid"], curationCell[ , "GRAY.cellid"])]
    
    
    
    #KALLISTO + SALMON
    
    if(rnaseq_tool=="kallistosalmon"){
      
      if (Kallisto1Pipeline != 0) {
        KallistoPipeline <- Kallisto1Pipeline
        tag1 <- "kallisto0.44.0"
        
      } 
      
      if (Kallisto2Pipeline != 0){
        
        KallistoPipeline <- Kallisto2Pipeline
        tag1 <- "kallisto0.43.1"
      }
      
      if (Salmon1Pipeline != 0) {
        SalmonPipeline <- Salmon1Pipeline
        tag2 <- "salmon0.11.3"
        
      } 
      
      if (Salmon2Pipeline != 0){
        
        SalmonPipeline <- Salmon2Pipeline
        tag2 <- "salmon0.11.2"
      }
      
      kallisto <- summarizeRnaSeqKallisto(dir="/Users/anthonymammoliti/Desktop/kallisto", tool="kallisto", features_annotation="/Users/anthonymammoliti/Desktop/Gencode.v23.annotation.RData", samples_annotation=rnaseq.sampleinfo)
      salmon <- summarizeRnaSeqSalmon(dir="/Users/anthonymammoliti/Desktop/salmon", tool="salmon", features_annotation="/Users/anthonymammoliti/Desktop/Gencode.v23.annotation.RData", samples_annotation=rnaseq.sampleinfo)
      
      z <- c(z,c(
        "kallisto_rnaseq"=kallisto[["gene_exp"]],
        "kallisto_rnaseq.counts"= kallisto[["gene_count"]],
        "kallisto_isoforms"=kallisto[["transcript_exp"]],
        "kallisto_isoforms.counts"=kallisto[["transcript_count"]],
        "salmon_rnaseq"=salmon[["gene_exp"]],
        "salmon_rnaseq.counts"= salmon[["gene_count"]],
        "salmon_isoforms"=salmon[["transcript_exp"]],
        "salmon_isoforms.counts"=salmon[["transcript_count"]])
      )
    } 
    
    
    
    
    #KALLISTO 1 (Version 0.44.0) & #KALLISTO 2 (Version 0.43.1)
    
    if(rnaseq_tool=="kallistocombined"){
      rnaseq1 <- summarizeRnaSeqKallisto(dir=Kallisto1Pipeline, tool="kallisto", features_annotation=param_transcript, samples_annotation=rnaseq.sampleinfo)
      rnaseq2 <- summarizeRnaSeqKallisto(dir=Kallisto2Pipeline, tool="kallisto", features_annotation=param_transcript, samples_annotation=rnaseq.sampleinfo)
      
      z <- c(z,c(
        "kallisto0.44.0_rnaseq"=rnaseq1[["gene_exp"]],
        "kallisto0.44.0_rnaseq.counts"= rnaseq1[["gene_count"]],
        "kallisto0.44.0_isoforms"=rnaseq1[["transcript_exp"]],
        "kallisto0.44.0_isoforms.counts"=rnaseq1[["transcript_count"]],
        "kallisto0.43.1_rnaseq"=rnaseq2[["gene_exp"]],
        "kallisto0.43.1_rnaseq.counts"= rnaseq2[["gene_count"]],
        "kallisto0.43.1_isoforms"=rnaseq2[["transcript_exp"]],
        "kallisto0.43.1_isoforms.counts"=rnaseq2[["transcript_count"]])
      )
      
      
    } 
    
    #KALLISTO 1 (Version 0.44.0)
    
    
    if(rnaseq_tool=="kallisto1"){
      rnaseq1 <- summarizeRnaSeqKallisto(dir="/Users/anthonymammoliti/Desktop/kallisto", tool="kallisto", features_annotation="/Users/anthonymammoliti/Desktop/Gencode.v23.annotation.RData", samples_annotation=rnaseq.sampleinfo)
      
      z <- c(z,c(
        "kallisto0.44.0_rnaseq"=rnaseq1[["gene_exp"]],
        "kallisto0.44.0_rnaseq.counts"= rnaseq1[["gene_count"]],
        "kallisto0.44.0_isoforms"=rnaseq1[["transcript_exp"]],
        "kallisto0.44.0_isoforms.counts"=rnaseq1[["transcript_count"]])
      )
      
      
    } 
    
    
    #KALLISTO 2 (Version 0.43.1)
    
    
    if(rnaseq_tool=="kallisto2"){
      rnaseq2 <- summarizeRnaSeqKallisto(dir="/Users/anthonymammoliti/Desktop/43", tool="kallisto", features_annotation="/Users/anthonymammoliti/Desktop/Gencode.v23.annotation.RData", samples_annotation=rnaseq.sampleinfo)
      
      z <- c(z,c(
        "kallisto0.43.1_rnaseq"=rnaseq2[["gene_exp"]],
        "kallisto0.43.1_rnaseq.counts"= rnaseq2[["gene_count"]],
        "kallisto0.43.1_isoforms"=rnaseq2[["transcript_exp"]],
        "kallisto0.43.1_isoforms.counts"=rnaseq2[["transcript_count"]])
      )
    } 
    
    
    #SALMON 1 (Version 0.11.3)
    
    if(rnaseq_tool=="salmon1"){
      salmon1 <- summarizeRnaSeqSalmon(dir="/Users/anthonymammoliti/Desktop/salmon", tool="salmon", features_annotation="/Users/anthonymammoliti/Desktop/Gencode.v23.annotation.RData", samples_annotation=rnaseq.sampleinfo)
      
      z <- c(z,c(
        "salmon0.11.3_rnaseq"=salmon1[["gene_exp"]],
        "salmon0.11.3_rnaseq.counts"= salmon1[["gene_count"]],
        "salmon0.11.3_isoforms"=salmon1[["transcript_exp"]],
        "salmon0.11.3_isoforms.counts"=salmon1[["transcript_count"]])
      )
    } 
    
    
    #SALMON 2 (Version 0.11.2)
    
    if(rnaseq_tool=="salmon2"){
      salmon2 <- summarizeRnaSeqSalmon(dir=Salmon2Pipeline, tool="salmon", features_annotation=param_transcript, samples_annotation=rnaseq.sampleinfo)
      
      z <- c(z,c(
        "salmon0.11.2_rnaseq"=salmon2[["gene_exp"]],
        "salmon0.11.2_rnaseq.counts"= salmon2[["gene_count"]],
        "salmon0.11.2_isoforms"=salmon2[["transcript_exp"]],
        "salmon0.11.2_isoforms.counts"=salmon2[["transcript_count"]])
      )
    } 
    
    
    
    #SALMON 1 (Version 0.11.3) & SALMON 2 (Version 0.11.2)
    
    if(rnaseq_tool=="salmoncombined"){
      salmon1 <- summarizeRnaSeqSalmon(dir=Salmon1Pipeline, tool="salmon", features_annotation=param_transcript, samples_annotation=rnaseq.sampleinfo)
      salmon2 <- summarizeRnaSeqSalmon(dir=Salmon2Pipeline, tool="salmon", features_annotation=param_transcript, samples_annotation=rnaseq.sampleinfo)
      
      z <- c(z,c(
        "salmon0.11.3_rnaseq"=salmon1[["gene_exp"]],
        "salmon0.11.3_rnaseq.counts"= salmon1[["gene_count"]],
        "salmon0.11.3_isoforms"=salmon1[["transcript_exp"]],
        "salmon0.11.3_isoforms.counts"=salmon1[["transcript_count"]],
        "salmon0.11.2_rnaseq"=salmon2[["gene_exp"]],
        "salmon0.11.2_rnaseq.counts"= salmon2[["gene_count"]],
        "salmon0.11.2_isoforms"=salmon2[["transcript_exp"]],
        "salmon0.11.2_isoforms.counts"=salmon2[["transcript_count"]])
      )
    } 
    
    
    
    
    
    
    ## cell information
    
    cellineinfo <- read.xlsx("/Users/anthonymammoliti/Desktop/GRAYLOCAL/cellline/gb-2013-14-10-r110-s1.xlsx", sheet = 1)
    cellineinfo[!is.na(cellineinfo) & cellineinfo == ""] <- NA
    rn <- cellineinfo[-1, 1]
    cn <- t(cellineinfo[1, -1])
    cn <- gsub(badchars, ".", cn)
    cellineinfo <- cellineinfo[-1, -1]
    dimnames(cellineinfo) <- list(rn, cn)
    cellineinfo <- data.frame("cellid"=rn, "tissueid"="breast", cellineinfo[,1:10])
    cellineinfo <- cellineinfo[which(!is.na(cellineinfo$Transcriptional.subtype)), ]
    #cellineinfo$cellid <- rownames(curationCell)[match(cellineinfo$cellid, curationCell[ , "GRAY.cellid"])]
    c1 <- matchToIDTableCELL(cellineinfo$cellid, curationCell, "GRAY.cellid")
    c1 <- as.character(c1)
    cellineinfo$cellid <- c1
    cellineinfo$cellid[is.na(cellineinfo$cellid)]<-"NA"
    rownames(cellineinfo) <-  cellineinfo$cellid
    cellineinfo <- cellineinfo[rownames(curationCell), ]
    head(cellineinfo)
    
    ##drug information
    druginfo <- data.frame("drugid"=curationDrug$unique.drugid)
    rownames(druginfo) <- druginfo$drugid
    
    
    if (sensitivity == 2013) {
      
      ##sensitivity
      
      load("/Users/anthonymammoliti/Desktop/GRAYLOCAL/Output/drug_norm_post.RData")
      
      
      sensitivity.info <- raw.sensitivity[ , c(1,2, 3, grep(tt, colnames(raw.sensitivity)))]
      colnames(sensitivity.info) <- c("cellid", "drugid", "min.Dose.uM", "max.Dose.uM")
      sensitivity.info <- cbind(sensitivity.info, "nbr.conc.tested"=con_tested)
      
      matchToIDTableCELL <- function(ids,tbl, column) {
        sapply(ids, function(x) {
          myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
          if(length(myx) > 1){
            stop("Something went wrong in curating cell ids")
          }
          return(tbl[myx, "unique.cellid"])
        })
      }
      
      x <- matchToIDTableCELL(sensitivity.info[, "cellid"], curationCell, "GRAY.cellid")
      #sensitivity.info[, "cellid"] <- rownames(curationCell)[match( sensitivity.info[, "cellid"], curationCell[ , "GRAY.cellid"])]
      x <- as.character(x)
      sensitivity.info[, "cellid"] <- x
      
      
      matchToIDTableDRUG <- function(ids,tbl, column) {
        sapply(ids, function(x) {
          myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
          if(length(myx) > 1){
            stop("Something went wrong in curating drug ids")
          }
          return(tbl[myx, "unique.drugid"])
        })
      }
      
      sensitivity.info[,"drugid"] <- gsub("\\s*\\([^\\)]+\\)","",sensitivity.info[,"drugid"])
      x2 <- matchToIDTableDRUG(sensitivity.info[, "drugid"], curationDrug, "GRAY.drugid")
      x2 <- as.character(x2)
      sensitivity.info[, "drugid"] <- x2
      #sensitivity.info[, "drugid"] <- rownames(curationDrug)[match( sensitivity.info[, "drugid"], curationDrug[ , "GRAY.drugid"])]
      
      raw.sensitivity <- raw.sensitivity[ ,-c(1,2)]
      raw.sensitivity <- array(c(as.matrix(raw.sensitivity[ ,1:con_tested]), as.matrix(raw.sensitivity[ ,(con_tested+1):(2*con_tested)])), c(nrow(raw.sensitivity), con_tested, 2),
                               dimnames=list(rownames(raw.sensitivity), colnames(raw.sensitivity[ ,1:con_tested]), c("Dose", "Viability")))
      
      load("/Users/anthonymammoliti/Desktop/GRAYLOCAL/Output/GRAY_sens_recomputed.RData")
      profiles <- read.xlsx("/Users/anthonymammoliti/Desktop/GRAYLOCAL/sens_published/gb-2013-14-10-r110-s1.xlsx", sheet = 1)
      profiles[!is.na(profiles) & profiles == ""] <- NA
      rn <- profiles[-1, 1]
      cn <- t(profiles[1, -1])
      profiles <- profiles[-1, -1]
      dimnames(profiles) <- list(rn, cn)
      profiles <- profiles[which(!is.na(profiles[, "Transcriptional subtype"])), ]
      colnames(profiles)[which(colnames(profiles) == "L-779450")] <-  "L-779405"
      indices <- 11:ncol(profiles)
      GI50 <- as.numeric(array(apply(profiles, 1, function(x)(x[indices]))))
      drugs <- rownames(curationDrug)[match(colnames(profiles)[indices], curationDrug[ , "GRAY.drugid"])]
      celllines <- rownames(curationCell)[match(rownames(profiles), curationCell[ , "GRAY.cellid"])]
      x <- expand.grid(drugs,celllines)
      names(GI50) <- paste("drugid", paste(x[,1],x[,2],sep = "_"), sep="_")
      GI50 <- GI50[which(!is.na(GI50))]
      sensitivity.profiles <- matrix(NA, dimnames = list(rownames(sensitivity.info), "GI50_published"), nrow=nrow(sensitivity.info))
      for(nn in names(GI50)) {
        sensitivity.profiles[grep(nn, rownames(sensitivity.profiles)), "GI50_published"] <- GI50[nn]
      }
      
      
      
      #sensitivity.profiles <- cbind(sensitivity.profiles, "auc_recomputed"=recomputed$AUC, "ic50_recomputed"=recomputed$IC50)
      
      slope <- NULL
      for(exp in rownames(sensitivity.info)){
        slope <- c(slope, computeSlope(raw.sensitivity[exp, , "Dose"], raw.sensitivity[exp, , "Viability"])) #computeSlope (returns normalized slope of drug response curve)
      }
      
      names(slope) <- rownames(sensitivity.info)
      sensitivity.profiles <- cbind(sensitivity.profiles, "slope_recomputed"=slope)
      head(sensitivity.profiles)
      
    }
    
    
    if (sensitivity == 2017) {
      
      gray.GRvalues <- read.csv("/Users/anthonymammoliti/Desktop/DrugSensNew/DS2_datafile.tsv", header=T, stringsAsFactors=FALSE, sep="\t")
      
      dd <- which(duplicated(gray.GRvalues))
      gray.GRvalues <- gray.GRvalues[-dd,]
      gray.GRvalues <- gray.GRvalues[-which(gray.GRvalues[, "Small.Molecule.HMS.LINCS.ID.1"]==""),]
      
      load("/Users/anthonymammoliti/Desktop/drug_norm_post2017.RData")
      load("/Users/anthonymammoliti/Desktop/GRAY_sens_recomputed2017.RData")
      
      sensitivity.profiles <- cbind("auc_recomputed"=recomputed$AUC/100, "ic50_recomputed"=recomputed$IC50)
      rownames(sensitivity.profiles) <- rownames(raw.sensitivity)
      
      dim(sensitivity.info)
      
      ###First need to cross reference cells/drugs(PERTURBAGENS) to be able to do matching
      cells.cross.reference <- read.csv("/Users/anthonymammoliti/Desktop/DrugSensNew/DS0_crossreferencingCELLS.txt", header=T, stringsAsFactors=FALSE, sep="\t")
      drugs.cross.reference <- read.csv("/Users/anthonymammoliti/Desktop/DrugSensNew/DS0_crossreferencingPERTURBAGENS.txt", header=T, stringsAsFactors=FALSE, sep="\t")
      
      remove.items <- which(sensitivity.info[,"cellid"] %in% cells.cross.reference[which(cells.cross.reference$COMMENT == "REMOVE"), "HEISER.NAME"])
      sensitivity.info <- sensitivity.info[-remove.items,]
      sensitivity.profiles <- sensitivity.profiles[-remove.items,]
      raw.sensitivity <- raw.sensitivity[-remove.items,,]
      sensitivity.info[, "cellid"] <- cells.cross.reference$LINCS.NAME[match( sensitivity.info[, "cellid"], cells.cross.reference[ , "HEISER.NAME"])]
      
      ##check
      dim(sensitivity.info)
      length(intersect(unique(gray.GRvalues[, "Cell.Name"]), unique(sensitivity.info[, "cellid"])))
      length(setdiff(unique(gray.GRvalues[, "Cell.Name"]), unique(sensitivity.info[, "cellid"])))
      length(setdiff(unique(sensitivity.info[, "cellid"]), unique(gray.GRvalues[, "Cell.Name"])))
      
      remove.items.dd <- which(sensitivity.info[,"drugid"] %in% drugs.cross.reference[which(drugs.cross.reference$COMMENT == "REMOVE"), "HEISER.NAME"])
      sensitivity.info <- sensitivity.info[-remove.items.dd,]
      sensitivity.profiles <- sensitivity.profiles[-remove.items.dd,]
      raw.sensitivity <- raw.sensitivity[-remove.items.dd,,]
      #sensitivity.info[, "drugid"] <- drugs.cross.reference$HMSLID..SALT.ID..BATCH.ID[match( sensitivity.info[, "drugid"], drugs.cross.reference[ , "HEISER.NAME"])]
      
      ##check
      dim(sensitivity.info)
      length(intersect(unique(gray.GRvalues[, "Small.Molecule.HMS.LINCS.ID"]), unique(sensitivity.info[, "drugid"])))
      length(setdiff(unique(gray.GRvalues[, "Small.Molecule.HMS.LINCS.ID"]), unique(sensitivity.info[, "drugid"])))
      ###two compounds in gray.GRvalues have not been in cross reference file "10110-101-1", "10252-101-1"
      length(setdiff(unique(sensitivity.info[, "drugid"]), unique(gray.GRvalues[, "Small.Molecule.HMS.LINCS.ID"])))
      
      #Add New SUM190PT Cell line & Tissue
      
      curationCell <- rbind(curationCell, c("SUM190PT","SUM190PT")) 
      rownames(curationCell)[85] <- "SUM190PT"
      
      curationTissue <- rbind(curationTissue, c("breast","breast")) 
      rownames(curationTissue)[85] <- "SUM190PT"
      
      
      cellineinfo <- rbind(cellineinfo, c("SUM190PT", "breast", "NA", "NA", "NA","NA","NA","NA","NA","NA","NA" ))
      rownames(cellineinfo)[85] <- "SUM190PT"
      
      sensitivity.info[,"cellid"] <- gsub("-", "", sensitivity.info[,"cellid"])
      sensitivity.info[,"cellid"] <- gsub(" ", "", sensitivity.info[,"cellid"])
      sensitivity.info[,"cellid"] <- toupper(sensitivity.info[,"cellid"])
      x2 <- matchToIDTableCELL(sensitivity.info[,"cellid"], curationCell, "GRAY.cellid")
      x2 <- as.character(x2)
      sensitivity.info[, "cellid"] <- x2
      
      sensitivity.info[,"drugid"] <- gsub("\\s*\\([^\\)]+\\)","",sensitivity.info[,"drugid"])
      
      matchToIDTableDRUG <- function(ids,tbl, column) {
        sapply(ids,function(x) {
          myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
          if(length(myx) > 1){
            stop("Something went wrong in curating drug ids")
          }
          return(tbl[myx, "unique.drugid"])
        })
      }
      
      x <- matchToIDTableDRUG(sensitivity.info[, "drugid"], curationDrug, "GRAY.drugid")
      
      #i <- sapply(x, is.factor)
      #x[i] <- lapply(x[i], as.character)
      #sensitivity.info[, "drugid"] <- rownames(curationDrug)[match(sensitivity.info[, "drugid"], curationDrug[ , "GRAY.drugid"])]
      
      y <- as.character(x)
      sensitivity.info[, "drugid"] <- y
      
      replicates <- unique(paste(gray.GRvalues[, "Cell.Name"], gray.GRvalues[, "Small.Molecule.HMS.LINCS.ID"], sep="_"))
      replicates.ll <-NULL
      i=1
      for(replicate in replicates){
        i=i+1
        cell <- strsplit(replicate, "_")[[1]][1]
        drug <- strsplit(replicate, "_")[[1]][2]
        xx <- which(gray.GRvalues[, "Cell.Name"] == cell& gray.GRvalues[, "Small.Molecule.HMS.LINCS.ID"]==drug)
        yy <- which(sensitivity.info[,"cellid"] == cell & sensitivity.info[,"drugid"]==drug)
        if(length(xx) > 5){
          #browser()
          #print(length(xx))
          if(length(xx) == length(yy)){
            replicates.ll <- c(replicates.ll, length(xx))
          }
        }
      }
      
      
      
      
      
      
      
      
      
    }
    
    
    
    GRAY <- PharmacoSet(molecularProfiles=z,
                        
                        
                        name="GRAY", 
                        cell=cellineinfo, 
                        drug=druginfo, 
                        sensitivityInfo=sensitivity.info, 
                        sensitivityRaw=raw.sensitivity, 
                        sensitivityProfiles=sensitivity.profiles, 
                        sensitivityN=NULL,
                        curationCell=curationCell, 
                        curationDrug=curationDrug, 
                        curationTissue=curationTissue, 
                        datasetType="sensitivity")
    save(GRAY,file="/Users/anthonymammoliti/Desktop/GRAYLOCAL/GRAY.RData")
    
    
    return (GRAY)
    
  }


summarizeRnaSeqKallisto <- function (dir, 
                                     tool="kallisto", 
                                     features_annotation,
                                     samples_annotation) {
  
  load(features_annotation)
  tx2gene <- as.data.frame(cbind("transcript"=toil.transcripts$transcript_id, "gene"=toil.transcripts$gene_id))
  #dir = "/Users/anthony/Desktop/test
  #dir = "/Users/anthony/Desktop/CCLE/ge/kallisto"
  files <- list.files(dir, recursive = TRUE, full.names = T)
  resFiles <- grep("abundance.h5", files)
  resFiles <- files[resFiles]
  length(resFiles)
  #txi.kallisto <- tximport(resFiles, type = "kallisto", txOut = TRUE) abundance files
  #txi.kallisto <- tximport(resFiles, type = "kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE) tsv files
  
  txi <- tximport(resFiles, type="kallisto", tx2gene=tx2gene, ignoreAfterBar = TRUE)
  head(txi$counts[,1:5])
  dim(txi$counts)
  
  xx <- txi$abundance
  colnames(xx) <- samples_annotation$cellid
  gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
  fData(gene.exp) <- toil.genes[featureNames(gene.exp),]
  pData(gene.exp) <- samples_annotation
  annotation(gene.exp) <- "rnaseq"
  
  xx <- txi$counts
  colnames(xx) <- samples_annotation$cellid
  gene.count <- Biobase::ExpressionSet(log2(xx + 1))
  fData(gene.count) <- toil.genes[featureNames(gene.count),]
  pData(gene.count) <- samples_annotation
  annotation(gene.count) <- "rnaseq"
  
  txii <- tximport(resFiles, type="kallisto", txOut=T)
  
  xx <- txii$abundance
  colnames(xx) <- samples_annotation$cellid
  transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
  fData(transcript.exp) <- toil.transcripts[featureNames(transcript.exp),]
  pData(transcript.exp) <- samples_annotation
  annotation(transcript.exp) <- "isoforms"
  
  xx <- txii$counts
  colnames(xx) <- samples_annotation$cellid
  transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
  fData(transcript.count) <- toil.transcripts[featureNames(transcript.count),]
  pData(transcript.count) <- samples_annotation
  annotation(transcript.count) <- "isoforms"
  
  return(list("gene_exp"=gene.exp, 
              "gene_count"=gene.count, 
              "transcript_exp"=transcript.exp, 
              "transcript_count"=transcript.count))
}



summarizeRnaSeqSalmon <- function (dir, 
                                   tool="salmon", 
                                   features_annotation,
                                   samples_annotation) {
  
  load(features_annotation)
  
  tx2gene <- as.data.frame(cbind("transcript"=toil.transcripts$transcript_id, "gene"=toil.transcripts$gene_id))
  
  files <- list.files(dir, recursive = TRUE, full.names = T)
  resFiles <- grep("quant.sf", files)
  resFiles <- files[resFiles]
  length(resFiles)
  
  txi <- tximport(resFiles, type="salmon", tx2gene=tx2gene, ignoreAfterBar = TRUE)
  dim(txi$counts)
  
  xx <- txi$abundance
  
  gene.exp <- Biobase::ExpressionSet(log2(xx + 0.001))
  fData(gene.exp) <- toil.genes[featureNames(gene.exp),]
  pData(gene.exp) <- samples_annotation
  annotation(gene.exp) <- "rnaseq"
  
  xx <- txi$counts
  
  gene.count <- Biobase::ExpressionSet(log2(xx + 1))
  fData(gene.count) <- toil.genes[featureNames(gene.count),]
  pData(gene.count) <- samples_annotation
  annotation(gene.count) <- "rnaseq"
  
  txii <- tximport(resFiles, type="salmon", txOut=T, ignoreAfterBar = TRUE)
  
  
  rownames(txii$abundance) <- gsub("\\|.*","",rownames(txii$abundance)) 
  
  xx <- txii$abundance
  
  transcript.exp <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 0.001))
  fData(transcript.exp) <- toil.transcripts[featureNames(transcript.exp),] #ISSUE
  pData(transcript.exp) <- samples_annotation
  annotation(transcript.exp) <- "isoforms"
  
  xx <- txii$counts
  
  transcript.count <- Biobase::ExpressionSet(log2(xx[,1:length(resFiles)] + 1))
  
  pData(transcript.count) <- samples_annotation
  annotation(transcript.count) <- "isoforms"
  
  return(list("gene_exp"=gene.exp, 
              "gene_count"=gene.count, 
              "transcript_exp"=transcript.exp, 
              "transcript_count"=transcript.count))
}




library(PharmacoGxPrivate)
library(PharmacoGx)
library(readr)
library(tximport)
library(rhdf5)
library(gdata)
library(readxl)
library(openxlsx)
library(CoreGx)






eval(parse(text = paste("getGRAYP(", finalparam, ")")))




