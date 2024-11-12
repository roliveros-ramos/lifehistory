
.vec2date = function(x, frequency=12) {
  if(length(x)==2) {
    x[2] = (x[2]-1)/frequency
  }
  if(length(x)==1) {
    out = numeric(2)
    out[1] = as.integer(x)
    out[2] = x - out[1]
    x = out
  }
  start_date = ymd(sprintf("%d-01-01", x[1]))
  yday(start_date) = floor((365+leap_year(start_date))*x[2]) + 1
  return(start_date)
}
