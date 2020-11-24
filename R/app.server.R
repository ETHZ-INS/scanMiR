#' scanMiR.server
#'
#' @param modlists A named list of `CompressedKdModelList`
#' @param targetlists An optional list of aggregated targets, named as `modlists`
#' @param ensdbs A named list of `ensembldb` objects (with names corresponding to those 
#' of `modlists`)
#' @param genomes A named list of `BSgenome` objects (with names corresponding to those 
#' of `modlists`)
#' @param maxCacheSize Maximum cache size in bytes.
#'
#' @return A shiny server function
#' @importFrom digest digest
#' @export
scanMiR.server <- function(modlists, targetlists=list(), ensdbs=list(), genomes=list(),
                           scans=list(), maxCacheSize=100*10^6, 
                           BP=BiocParallel::SerialParam() ){
  library(DT)
  library(stringr)
  library(Biostrings)
  library(ggplot2)
  
  dtwrapper <- function(d, pageLength=25, ...){
    datatable( d, filter="top", class="compact", extensions=c("Buttons","ColReorder"),
               options=list(pageLength=pageLength, dom = "fltBip", rownames=FALSE,
                            colReorder=TRUE, 
                            buttons=c('copy', 'csv', 'excel', 'csvHtml5') ), ... )
  }
  
  function(input, output, session){
    
    ##############################
    ## initialize inputs
    
    updateSelectizeInput(session, "mirlist", choices=names(modlists))
    updateSelectizeInput(session, "annotation", choices=names(ensdbs))
    
    observe({
      # when the choice of collection changes, update the annotation to
      # use the same genome
      if(!is.null(input$mirlist) && input$mirlist!="")
        updateSelectizeInput(session, "annotation", selected=input$mirlist)
    })
    
    
    ##############################
    ## select collection
    
    allmods <- reactive({ # all models from collection
      if(is.null(input$mirlist)) return(NULL)
      modlists[[input$mirlist]]
    })
    
    # prints a summary of the model collection
    output$collection_summary <- renderPrint({
      if(is.null(allmods())) return(NULL)
      summary(allmods())
    })
    
    observe({ # when the selected collection changes, update the miRNA selection inputs
      updateSelectizeInput(session, "mirnas", choices=names(allmods()), server=TRUE)
      updateSelectizeInput(session, "mirna", choices=names(allmods()), server=TRUE)
    })
    
    ##############################    
    ## scan specific sequence
    
    ## transcript selection
    
    sel_ensdb <- reactive({ # the ensembldb for the selected genome
      if(is.null(input$annotation) || input$annotation=="" || 
         !(input$annotation %in% names(ensdbs))) return(NULL)
      ensdbs[[input$annotation]]
    })
    
    allgenes <- reactive({ # all genes in the selected genome
      if(is.null(sel_ensdb())) return(NULL)
      g <- genes( sel_ensdb(), columns="gene_name",
                  return.type="data.frame")
      gs <- g[,2]
      names(gs) <- paste(g[,1], g[,2])
      gs
    })
    
    selgene <- reactive({ # selected gene id
      if(is.null(input$gene) || input$gene=="") return(NULL)
      changeFlag()
      input$gene
    })
    
    output$gene_link <- renderUI({
      if(is.null(selgene())) return(NULL)
      tags$a(href=paste0("https://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=", selgene()),
             "view on ensembl", target="_blank")
    })
    
    alltxs <- reactive({ # all tx from selected gene
      if(is.null(selgene())) return(NULL)
      tx <- transcripts(sel_ensdb(), columns=c("tx_id","tx_biotype"),
                        filter=~gene_id==selgene(), return.type="data.frame")
      txs <- tx$tx_id
      names(txs) <- paste0(tx$tx_id, " (", tx$tx_biotype,")")
      changeFlag()
      txs
    })
    
    seltx <- reactive({ # the selected transcript
      if(is.null(input$transcript) || input$transcript=="" || is.na(input$transcript))
        return(NULL)
      changeFlag()
      input$transcript
    })
    
    # when the ensembldb is updated, update the gene input
    observe(updateSelectizeInput(session, "gene", choices=allgenes(), server=TRUE))
    # when the gene selection is updated, update the transcript input
    observe(updateSelectizeInput(session, "transcript", choices=alltxs()))
    
    # takes a genome package name as input, and returns the genome
    getGenome <- function(x){
      if(is.character(x)){
        library(x, character.only=TRUE)
        x <- get(x)
      }
      seqlevels(x) <- gsub("^chr","",seqlevels(x))
      x
    }
    
    seqs <- reactive({ # returns the selected sequence(s)
      if(is.null(selgene()) && is.null(seltx())) return(NULL)
      if(is.null(seltx())){
        gid <- selgene()
        filt <- ~gene_id==gid
      }else{
        txid <- seltx()
        filt <- ~tx_id==txid
      }
      if(input$utr_only){
        gr <- suppressWarnings(threeUTRsByTranscript(sel_ensdb(), filter=filt))
      }else{
        gr <- exonsBy(sel_ensdb(), by="tx", filter=filt)
      }
      get_seq(gr)
    })
    
    get_seq <- function(gr){
      if(length(gr)==0) return(NULL)
      if(is.null(genomes[[input$annotation]])) return(NULL)
      seqs <- extractTranscriptSeqs(getGenome(genomes[[input$annotation]]), gr)
      seqs <- seqs[lengths(seqs)>6]
      if(length(seqs)==0) return(NULL)
      seqs
    }
    
    output$tx_overview <- renderTable({ # overview of the selected transcript
      if(is.null(seqs()) || length(seqs())==0) return(NULL)
      w <- width(seqs())
      ss <- as.character(subseq(seqs(), 1, end=ifelse(w<40,w,40)))
      ss[w>40] <- paste0(ss[w>40],"...")
      data.frame(transcript=names(seqs()), length=w, sequence=ss)
    })
    
    ## end transcript selection

    ## custom sequence
    
    customTarget <- reactive({
      if(is.null(input$customseq) || input$customseq=="") return(NULL)
      seqtype <- suppressWarnings(scanMiR:::.guessSeqType(input$customseq))
      seq <- gsub("[^ACGTUN]","", toupper(input$customseq))
      if(input$circular) seq <- paste0(seq,substr(seq,1,min(nchar(seq),11)))
      if(seqtype=="DNA") return(DNAString(seq))
      return(RNAString(seq))
    })
    
    output$custom_info <- renderPrint({ # overview of the custom sequence
      if(is.null(input$customseq)) return("")
      out <- capture.output(customTarget())
      if(input$circular) out <- c("Circularized sequence: the first 11nt are pasted to ",
                                  "the end of the sequence.\n", out)
      cat(out)
    })

    target <- reactive({ # target subject sequence
      if(input$subjet_type=="custom"){
        return(as.character(DNAStringSet(customTarget())))
      }   
      if(is.null(seqs()) || length(seqs())>1) return(NULL)
      changeFlag()
      as.character(seqs())
    })
    
    observeEvent(input$rndseq, { # generate random sequence
      updateTextAreaInput(session, "customseq",
                          value=paste(sample(c("A","C","G","T"), size = 3000, replace=TRUE), collapse=""))
    })
    
    ## Select miRNAs for scanning
    
    observeEvent(input$mirnas_confident, {
      if(is.null(allmods())) return(NULL)
      cons <- conservation(allmods())
      if(all(is.na(cons))) return(NULL)
      updateSelectizeInput(session, "mirnas", selected=names(cons)[as.numeric(cons)>1])
    })
    observeEvent(input$mirnas_mammals, {
      if(is.null(allmods())) return(NULL)
      cons <- conservation(allmods())
      if(all(is.na(cons))) return(NULL)
      updateSelectizeInput(session, "mirnas", selected=names(cons)[as.numeric(cons)>2])
    })
    observeEvent(input$mirnas_vert, {
      if(is.null(allmods())) return(NULL)
      cons <- conservation(allmods())
      if(all(is.na(cons))) return(NULL)
      updateSelectizeInput(session, "mirnas", selected=names(cons)[as.numeric(cons)>3])
    })
    observeEvent(input$mirnas_clear, {
      updateSelectizeInput(session, "mirnas", selected="")
    })
    
    selmods <- reactive({ # models selected for scanning
      if(is.null(allmods())) return(NULL)
      if(input$mirnas_all) return(allmods())
      if(is.null(input$mirnas)) return(NULL)
      allmods()[input$mirnas]
    })
    
    ## Begin scan and results caching
    
    cached.hits <- reactiveValues() # actual and past scanning results are stored in this object
    
    changeFlag <- reactiveVal(0)
    
    cached.checksums <- reactive({
      ch <- reactiveValuesToList(cached.hits)
      ch <- ch[!sapply(ch, is.null)]
      ch <- lapply(ch, FUN=function(x){
        x[c("target","size","time","last","nsel","sel")]
      })
    })
    current.cs <- reactiveVal()
    
    hits <- reactive({ # the results currently loaded are stored in this object
      if(is.null(current.cs())) return(NULL)
      if(current.cs() %in% names(cached.checksums())) return(cached.hits[[current.cs()]])
      NULL
    })
    
    checksum <- reactive({ # generate a unique hash for the given input
      changeFlag()
      paste( digest::digest(selmods()),
             digest::digest(list(target=target(), shadow=input$shadow,
                                 keepMatchSeq=input$keepMatchSeq, 
                                 minDist=input$minDist, maxLogKd=input$maxLogKd,
                                 scanNonCanonical=input$scanNonCanonical))
      )
    })
    
    cache.size <- reactive({
      ch <- cached.checksums()
      if(is.null(ch) || length(ch)==0) return(0)
      sum(sapply(ch, FUN=function(x) as.numeric(x$size)))
    })
    
    cleanCache <- function(){ # remove last-used results when over the cache size limit
      cs <- isolate(cached.checksums())
      if(length(cs)<3 || as.numeric(cache.size())<maxCacheSize) return(NULL)
      cs <- cs[order(sapply(cs, FUN=function(x) x$last), decreasing=TRUE)]
      sizes <- sapply(cs, FUN=function(x) as.numeric(x$size))
      while(length(cs)>2 & sum(sizes)>maxCacheSize){
        cached.hits[[rev(names(cs))[1]]] <- NULL
        cs <- cs[-length(cs)]
        sizes <- sizes[-length(sizes)]
      }
    }
    
    observeEvent(input$scan, { # actual scanning
      if(is.null(selmods()) || is.null(target()) || nchar(target())==0) 
        return(NULL)
      cs <- checksum()
      cached.hits[[cs]] <- do.scan()
      current.cs(cs)
    })
    
    do.scan <- reactive({
      if(is.null(selmods()) || is.null(target()) || nchar(target())==0) return(NULL)
      tmp <- changeFlag()
      cs <- checksum()
      if(cs %in% names(cached.checksums())) return(cached.hits[[cs]])
      res <- list()
      msg <- paste0("Scanning sequence for ",length(selmods())," miRNAS")
      message(msg)
      detail <- NULL
      if(length(selmods())>4) detail <- "This might take a while..."
      if(input$circular) detail <- "'Ribosomal Shadow' is ignored when scanning circRNAs"
      withProgress(message=msg, detail=detail, value=1, max=3, {
        res$hits = findSeedMatches(target(), selmods(),
                                   keepMatchSeq=input$keepmatchseq,
                                   minDist=input$minDist, maxLogKd=input$maxLogKd,
                                   shadow=ifelse(input$circular,0,input$shadow),
                                   onlyCanonical=!input$scanNonCanonical,
                                   BP=BP )
      })
      if(length(res$hits)>0) res$hits$log_kd <- (res$hits$log_kd/1000)
      res$cs <- cs
      res$last <- res$time <- Sys.time()
      res$size <- object.size(res$hits)
      res$nsel <- nm <- length(selmods())
      res$sel <- ifelse(nm>1,paste(nm,"models"),input$mirnas)
      res$target_length <- nchar(target())
      if(input$subjet_type=="custom"){
        res$target <- "custom sequence"
      }else{
        res$target <- paste(input$gene,"-", seltx(), ifelse(input$utr_only, "(UTR)",""))
      }
      return(res)
    })
    
    output$scan_target <- renderText({
      if(is.null(current.cs()) || is.null(cached.hits[[current.cs()]])) return(NULL)
      paste("Scan results in: ", cached.hits[[current.cs()]]$target)
    })
    
    output$cache.info <- renderText({
      if(cache.size()==0) return("Cache empty.")
      paste0(length(cached.checksums()), " results cached (",
             round(cache.size()/1024^2,3)," Mb)")
    })
    
    output$cached.results <- renderUI({
      ch <- cached.checksums()
      ch2 <- names(ch)
      names(ch2) <- sapply(ch, FUN=function(x){
        paste0(x$time, ": ", x$sel, " on ", x$target, " (", format(x$size,units="Kb"),")")
      })
      radioButtons("selected.cache", "Cached results", choices=ch2)
    })
    
    observeEvent(input$loadCache, {
      if(is.null(input$selected.cache)) return(NULL)
      current.cs(input$selected.cache)
    })
    
    observeEvent(input$deleteCache, {
      if(is.null(input$selected.cache)) return(NULL)
      cached.hits[[input$selected.cache]] <- NULL
      if(current.cs()==input$selected.cache) current.cs(NULL)
    })
    
    output$hits_table <- renderDT({ # prints the current hits
      if(is.null(hits()$hits)) return(NULL)
      h <- as.data.frame(hits()$hits)
      h <- h[,setdiff(colnames(h), c("seqnames","width","strand") )]
      if(hits()$nsel == 1) h$miRNA <- NULL
      dtwrapper(h)
    })
    
    output$dl_hits <- downloadHandler(
      filename = function() {
        if(is.null(hits()$hits)) return(NULL)
        fn <- paste0("hits-",gsub("\\.[09]+","",cached.hits[[current.cs()]]$target))
        if(hits()$nsel == 1){
          fn <- paste0(fn,"-",cached.hits[[cs]]$sel,".csv")
        }else{
          fn <- paste0(fn,"-",Sys.Date(),".csv")
        }
        fn
      },
      content = function(con) {
        if(is.null(hits()$hits)) return(NULL)
        h <- as.data.frame(hits()$hits)
        h <- h[,setdiff(colnames(h), c("seqnames","width","strand") )]
        write.csv(h, con, col.names=TRUE)
      }
    )
    
    ## end scan hits and cache 
    
    output$manhattan <- renderPlotly({
      if(is.null(hits()$hits)) return(NULL)
      h <- as.data.frame(sort(hits()$hits))
      if(!is.null(h$miRNA) && length(unique(h$miRNA))>input$manhattan_n){
        tt <- sort(table(h$miRNA), decreasing=TRUE)
        mirs <- names(tt)[seq_len(min(input$manhattan_n, length(tt)))]
        h <- h[h$miRNA %in% mirs,,drop=FALSE]
      }
      if(!is.null(input$manhattan_ordinal) && input$manhattan_ordinal){
        h$position <- seq_len(nrow(h))
        xlab <- "Position (ordinal)"
      }else{
        h$position <- round(rowMeans(h[,2:3]))
        xlab <- "Position (nt) in sequence"
      }
      xlim <- c(0,nchar(hits()$target_length))
      ael <- list(x="position", y="-log_kd", type="type")
      if("sequence" %in% colnames(h)) ael$seq="sequence"
      if("miRNA" %in% colnames(h)) ael$colour="miRNA"
      p <- ggplot(h, do.call(aes_string, ael)) + 
        geom_hline(yintercept=1.5, linetype="dashed", color = "red", size=1) + 
        geom_point(size=2) + xlab(xlab) + expand_limits(x=xlim, y=0)
      ggplotly(p)
    })
    
    ##############################
    ## miRNA-centric tab
    
    mod <- reactive({ # the currently-selected KdModel
      if(is.null(allmods()) || is.null(input$mirna)) return(NULL)
      allmods()[[input$mirna]]
    })
    
    output$modconservation <- renderText({
      if(is.null(mod())) return(NULL)
      as.character(conservation(mod()))
    })
    
    output$mirbase_link <- renderUI({
      tags$a(href=paste0("http://www.mirbase.org/textsearch.shtml?q=", input$mirna), 
             "miRBase", target="_blank")
    })
    
    output$modplot <- renderPlot({ # affinity plot
      if(is.null(mod())) return(NULL)
      plotKdModel(mod())
    }, height=reactive(input$modplot_height))
    
    txs <- reactive({ # the tx to gene symbol table for the current annotation
      if(is.null(input$mirlist) || input$mirlist=="" || 
         !(input$mirlist %in% names(ensdbs))) return(NULL)
      db <- ensdbs[[input$mirlist]]
      if(is.null(db)) return(NULL)
      tx <- mcols(transcripts(db, c("tx_id","gene_id","tx_biotype")))
      tx <- merge(tx,mcols(genes(db, c("gene_id","symbol"))), by="gene_id")
      as.data.frame(tx[,c("symbol","tx_id","tx_biotype")])
    })
    
    mirtargets_prepared <- reactive({
      tl <- paste0(input$mirlist, ifelse(input$targetlist_utronly, ".utrs", ".full"))
      if(is.null(targetlists[[tl]]) || is.null(mod())) return(NULL)
      d <- targetlists[[tl]][[input$mirna]]
      d$repression <- d$repression/1000
      if(!is.null(txs())){
        d <- merge(txs(), d, by.x="tx_id", by.y="transcript")
        if(input$targetlist_gene){
          d <- d[order(d$repression),]
          d <- d[!duplicated(d$symbol),]
          ff <- c("symbol","miRNA","repression")
        }
      }
      d[order(d$repression),]
    })
    
    output$mirna_targets <- renderDT({
      d <- mirtargets_prepared()
      if(is.null(d)) return(NULL)
      colnames(d) <- gsub("^n\\.","",colnames(d))
      dtwrapper(d, selection="single", callback=JS('
      table.on("dblclick.dt","tr", function() {
        Shiny.onInputChange("dblClickSubject", table.row(this).data()[1])
      })
    '))
    })
    
    ## double-click on a transcript in miRNA targets:
    observeEvent(input$dblClickSubject, {
      sub <- input$dblClickSubject
      if(input$targetlist_gene) return(NULL)
      tx <- transcripts(sel_ensdb(), columns=c("tx_id", "gene_id"),
                          filter=~tx_id==sub, return.type="data.frame")
      updateSelectizeInput(session, "gene", selected=as.character(tx$gene_id[1]), choices=allgenes(), server=TRUE)
      updateCheckboxInput(session, "utr_only", value=input$targetlist_utronly)
      updateSelectizeInput(session, "transcript", selected=sub, choices=alltxs())
      updateSelectizeInput(session, "mirnas", selected=input$mirna)
      updateTabItems(session, "subject_type", "transcript")
      updateTabItems(session, "main_tabs", "tab_subject")
      newflag <- changeFlag()+1
      changeFlag(newflag)
      observe({
        tmp <- changeFlag()
        cs <- checksum()
        cached.hits[[cs]] <- do.scan()
        current.cs(cs)
      })
      updateTabItems(session, "main_tabs", "tab_hits")
    })
    
    output$dl_mirTargets <- downloadHandler(
      filename = function() {
        if(is.null(input$mirna)) return(NULL)
        paste0(input$mirna, ifelse(input$targetlist_utronly,"-utr-","-"), "targets.csv")
      },
      content = function(con) {
        write.csv(mirtargets_prepared(), con, col.names=TRUE)
      }
    )
    
  }
}
