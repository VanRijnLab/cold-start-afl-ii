require(DBI)
require(RSQLite)

db_connect <- function(path) {
  db <- dbConnect(RSQLite::SQLite(), path)
  return(db)
}

db_get_query <- function(path, query, set_dt = TRUE) {
  db <- db_connect(path)
  out <- dbGetQuery(db, query)
  dbDisconnect(db)
  if (set_dt == TRUE) {
    setDT(out)
  }
  return(out)
}