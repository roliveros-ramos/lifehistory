
#' @export
check_taxon = function(x, na.return=FALSE, verbose=TRUE) {

  if(!requireNamespace("taxize", quietly = TRUE))
    stop("You need to install the 'taxize' package.")

  if(is.factor(x)) x = as.character(x)
  if(!is.character(x)) stop("x must be a character vector.")

  .checkScientificName = function(x, na.return, verbose) {
    if(is.na(x)) return(NA)
    tmp = taxize::gnr_resolve(sci = x, canonical = TRUE)$matched_name2
    tmp = names(which.max(table(tmp)))
    isNA = FALSE
    if(is.null(tmp)) {
      if(isTRUE(na.return)) {
        if(isTRUE(verbose)) message(sprintf("Name '%s' not found, returning NA.", x))
        tmp = NA_character_
        isNA = TRUE
      } else {
        if(isTRUE(verbose)) message(sprintf("Name '%s' not found, returning original name.", x))
        tmp = x
      }
    }
    id = identical(x, tmp)
    if(!id & !isNA & isTRUE(verbose)) message(sprintf("Taxon '%s' was corrected to '%s'.", x, tmp))
    return(tmp)
  }

  sx = na.omit(unique(x))

  out = unlist(sapply(sx, .checkScientificName,
                      na.return=na.return, verbose=verbose))

  out = out[x]
  names(out) = NULL

  if(is.null(out)) {
    # case?
  }

  return(out)
}

#' @export
get_taxon = function(x, rank, db="itis") {

  if(!requireNamespace("taxize", quietly = TRUE))
    stop("You need to install the 'taxize' package.")

  if(length(rank)!=1) stop("Only one taxon is allowed.")

  if(is.factor(x)) x = as.character(x)
  if(!is.character(x)) stop("x must be a character vector.")

  db = match.arg(db, c("itis", "ncbi", "both"))
  isna = is.na(x)
  xna = unique(x[!isna])
  tmp = taxize::tax_name(sci=xna, db=db, get=rank, messages=FALSE, ask=FALSE)[, rank]
  xind = which(is.na(tmp) & .is_binomial(xna))
  if(length(xind)>0) {
    last_x = .get_first(xna[xind])
    tmp[xind] = taxize::tax_name(query=last_x, db=db, get=rank, messages=FALSE, ask=FALSE)[, rank]
  }
  names(tmp) = xna
  out = tmp[x]
  names(out) = NULL
  return(out)

}

#' @export
get_classification = function(x, db="itis") {
  ddb = db
  if(db %in% c("fishbase", "fishlife")) ddb = "itis"
  ox = x
  x = check_taxon(x, verbose=FALSE, na.return = TRUE)
  if(is.na(x)) stop(sprintf("Taxon '%s' not found.", ox))
  taxa = c('kingdom', 'subkingdom', 'infrakingdom', 'phylum', 'subphylum',
           'infraphylum', 'superclass', 'class', 'superorder', 'order', 'suborder',
           'family', 'subfamily', 'genus', 'species')
  out = setNames(rep(NA, length(taxa)), nm=taxa)
  x = taxize::classification(x, db=ddb, messages=FALSE, ask=FALSE)[[1]]
  out[x$rank] = x$name
  out = switch(db,
               fishlife = .get_classification_fishlife(out),
               fishbase = .get_classification_fishbase(out),
               out)
  out = as.list(out)
  return(out)
}

# Internal ----------------------------------------------------------------

.is_binomial = function(x) {
  out = sapply(strsplit(x, split=" "), length) == 2
  return(out)
}

.get_first = function(x) {
  out = sapply(strsplit(x, split=" "), "[", i=1)
  return(out)
}


.get_classification_fishbase = function(x) {
  # dat = fishbase
  # names(dat) = tolower(names(dat))
  # dat = dat[, c("genus", "species", "subfamily", "family", "order", "class")]
  dat = .get_fishbase_list()
  out = x[c("species", "genus", "family", "order")]
  # if(!is.na(out[["species"]])) out[["species"]] = strsplit(out[["species"]], split=" ")[[1]][2]
  out = out[!sapply(out, FUN=is.na)]
  out = as.data.frame(as.list(out))
  out = merge(out, dat, all.x=TRUE)
  .myunique = function(x) {
    out = unique(x)
    if(length(out)==1) return(out)
    return(NA)
  }
  out = apply(out, 2, .myunique)
  # if(!is.na(out["species"])) out["species"] = paste(out["genus"], out["species"], sep=" ")
  out = list(rank=names(out), name=setNames(out, nm=NULL))
  x[] = NA
  x[out$rank] = out$name
  return(x)
}

.get_fishbase_list = function(server = getOption("FISHBASE_API", "fishbase"), limit=200000L) {
  taxa = c("genus", "species", "subfamily", "family", "order", "class", "superclass", "phylum", "kingdom")
  fishbase = suppressMessages(load_taxa(server=server, update = TRUE, limit=limit))
  fishbase = use_ascii(fishbase)
  fishbase = as.data.frame(lapply(fishbase, unlist))
  names(fishbase) = tolower(names(fishbase))
  if(is.null(fishbase$superclass)) fishbase$superclass = NA
  if(is.null(fishbase$phylum)) fishbase$phylum = "Chordata"
  if(is.null(fishbase$kingdom)) fishbase$kingdom = "Animalia"
  fishbase = fishbase[, taxa]
  return(fishbase)
}


.get_classification_fishlife = function(x, taxa=c("species", "genus", "family", "order")) {
  dat = fishlife
  out = x[taxa]
  if(!is.na(out[["species"]])) out[["species"]] = strsplit(out[["species"]], split=" ")[[1]][2]
  out = out[!sapply(out, FUN=is.na)]
  out = as.data.frame(as.list(out))
  out = merge(out, dat, all.x=TRUE)
  .myunique = function(x) {
    out = unique(x)
    if(length(out)==1) return(out)
    return(NA)
  }
  out = apply(out, 2, .myunique)
  if(is.na(out["class"])) {
    xo = as.list(x)
    xo$species = NA
    xout = .get_classification_fishlife(xo)
    xout$species = as.character(out["species"])
    out = unlist(xout[c(taxa, "class")])
  }
  # if(!is.na(out["species"])) out["species"] = paste(out["genus"], out["species"], sep=" ")
  out = list(rank=names(out), name=setNames(out, nm=NULL))
  x[] = NA
  x[out$rank] = out$name
  return(x)
}

# Internal functions ------------------------------------------------------

replace_non_ascii =function(string){

  .find_non_ascii = function(string){
    grep("I_WAS_NOT_ASCII",
         iconv(string, "latin1", "ASCII", sub="I_WAS_NOT_ASCII"))
  }

  i = .find_non_ascii(string)
  non_ascii = "áéíóúÁÉÍÓÚñÑüÜ’åôö"
  ascii = "aeiouAEIOUnNuU'aoo"
  translated = sapply(string[i], function(x)
    chartr(non_ascii, ascii, x))
  string[i] = unname(translated)
  return(string)
}

use_ascii = function(df){
  for(i in 1:length(df)){
    df[[i]] = replace_non_ascii(df[[i]])
  }
  return(df)
}

use_utf8 = function(df){
  chars = which(sapply(df, class) %in% "character")
  for(i in chars){
    df[[i]] = enc2utf8(df[[i]])
  }
  return(df)
}

