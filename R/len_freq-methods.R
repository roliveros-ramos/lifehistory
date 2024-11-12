
#' @export
window.len_freq = function(x, start, end, frequency=12, index.return=FALSE, ...) {
  start = .vec2date(start, frequency = frequency)
  end = .vec2date(end, frequency = frequency)
  ind = x$date >= start & x$date < end
  x$x = x$x[, ind]
  x$date = x$date[ind]
  return(x)
}

#' @export
Ops.len_freq = function(e1, e2) {
  if(inherits(e1, "len_freq") & !inherits(e2, "len_freq")) {
    # if(identical(dim(e1$x)[1:2], dim(e2))) e2 = as.numeric(e2)
    e1$x = get(.Generic, mode="function")(e1$x, e2)
    return(e1)
  }
  if(inherits(e2, "len_freq") & !inherits(e1, "len_freq")) {
    # if(identical(dim(e2$x)[1:2], dim(e1))) e1 = as.numeric(e1)
    e2$x = get(.Generic, mode="function")(e1, e2$x)
    return(e2)
  }
  if(inherits(e1, "len_freq") & inherits(e2, "len_freq")) {
      stop(gettextf("Operation '%s' not allowed for 'len_freq' objects.",
                    .Generic), domain = NA)
    }
}

#' @export
cbind.len_freq = function(..., deparse.level=1) {
  obj = list(...)
  xs = lapply(obj, FUN="[[", i="x")
  dims = unique(sapply(xs, function(x) dim(x)[1]))
  if(length(dims)!=1) stop("Dimensions do not match!")
  xs = do.call(cbind, xs)
  marks = lapply(obj, FUN="[[", i="marks")
  check = all(sapply(marks, identical, y=marks[[1]]))
  if(!check) stop("Length marks does not match!")
  dates = do.call(c, lapply(obj, FUN="[[", i="date"))
  out = lf_create(x=xs, marks=marks[[1]], date=dates)
  return(out)
}
