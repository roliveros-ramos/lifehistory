

#' @export
lf_create = function(x, marks, date, bins=NULL) {

  if(is.null(bins)) bins = .calculate_bins(marks)
  date = as.Date(date)
  out = list(x=x, marks=marks, date=date, bins=bins)
  class(out) = "len_freq"
  return(out)

}

