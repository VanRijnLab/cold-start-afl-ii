library(data.table)
library(fst)
library(DBI)
library(purrr)
library(stringr)

source(file.path("scripts", "99_database_functions.R"))

db_path <- file.path("data", "noordhoff.sqlite")



# Process data ------------------------------------------------------------

courses <- c("Grandes Lignes", "Stepping Stones")

map(courses, function (course) {
  
  message(course)
  
  responses_query <- paste0("SELECT * FROM responses_noduplicates WHERE method IN ('", course, "')")
  res <- db_get_query(db_path, responses_query)
  
  setorder(res, date, user_id, book_info_id, subsession, start_time)
  
  res[, answer_sanitised := str_replace_all(str_to_lower(str_squish(answer)), "[[:punct:]]", "")]
  
  # Fact IDs are recycled; they are only tied to a unique fact if we group by date, book_info_id, and subsession.
  res[, fact_id_uniq := paste(date, book_info_id, subsession, fact_id, sep = "_")]
  
  correct_answers <- res[correct == 1,
                         .(answer_sanitised = unique(answer_sanitised)),
                         by = .(fact_id_uniq, date, book_info_id, subsession, fact_id)]
  
  # There can still be multiple correct answers to a unique fact, but these are just variations in gender/spacing/optional bits/etc:
  correct_answers[, .(answer_variant = answer_sanitised[length(answer_sanitised) > 1]), by = .(fact_id_uniq)]
  
  # For those cases, assign the most common answer variant as the standardised answer to each fact ID
  correct_answers_std <- (
    res
    [correct == 1, .N, by = .(fact_id_uniq, date, book_info_id, subsession, fact_id, answer_sanitised)]
    [order(-N)]
    [, N := NULL]
    [, .SD[1], by = .(fact_id_uniq)]
  )
  
  # An answer is sometimes also associated with more than one fact_id, but we cannot be sure that the question was the same
  correct_answers[, .(fact_id = fact_id_uniq[length(fact_id_uniq) > 1]), by = .(answer_sanitised)]
  
  # There are many cases where just the date is different, but the book chapter info, the subsession, fact ID, and correct response are the same.
  # In these cases, try to combine into a single fact_id:
  correct_answers_std <- correct_answers_std[, .(fact_id_uniq_merged = paste0(fact_id_uniq[1], "_m"), 
                                                 fact_id_uniq = fact_id_uniq), 
                                             by = .(book_info_id, subsession, fact_id, answer_sanitised)]
  
  res <- merge(res, correct_answers_std[, .(fact_id_uniq, fact_id_uniq_merged)],
               by = c("fact_id_uniq"),
               all.x = TRUE)
  
  
  # Format columns ----------------------------------------------------------
  
  res <- res[,
             .(user_id,
               fact_id = fact_id_uniq_merged,
               start_time,
               rt,
               correct,
               choices,
               method)]
  
  
  
  write_fst(res, file.path("data", paste0("formatted_", str_replace_all(course, " ", "_"), ".fst")))
  
  rm(res, correct_answers, correct_answers_std)
  gc()
  
})