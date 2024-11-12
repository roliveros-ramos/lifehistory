

#' @export
lf_create = function(x, marks, date, bins=NULL) {

  if(is.null(bins)) bins = .calculate_bins(marks)
  date = as.Date(date)
  if(nrow(x)!=length(marks)) stop("Length-frequency date dimension does not match number of marks.")
  if(ncol(x)!=length(date)) stop("Length-frequency date dimension does not match number of dates.")

  x = x[, order(date)] # sorted dates

  to_keep = apply(x, 2, FUN=function(x) !all(x==0) | all(is.na(x)))

  out = list(x=x[, to_keep], marks=marks, date=date[to_keep], bins=bins)
  class(out) = "len_freq"
  return(out)

}

