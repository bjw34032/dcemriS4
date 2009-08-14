convert.datatype.anlz <- function(datatype) {
  switch(as.character(datatype),
         "0" = "UNKNOWN",
         "1" = "BINARY",
         "2" = "UNSIGNED_CHAR",
         "4" = "SIGNED_SHORT",
         "8" = "SIGNED_INT",
         "16" = "FLOAT",
         "32" = "COMPLEX",
         "64" = "DOUBLE",
         "128" = "RGB",
         "255" = "ALL")
}

convert.orient.anlz <- function(orientation) {
  switch(as.character(orientation),
         "0" = "transverse unflipped",
         "1" = "coronal unflipped",
         "2" = "sagittal unflipped",
         "3" = "transverse flipped",
         "4" = "coronal flipped",
         "5" = "sagittal flipped",
         "unknown")
}


