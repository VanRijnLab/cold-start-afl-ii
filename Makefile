.PHONY: all clean

all: output/04_evaluate_ROF_predictions.nb.html output/05_evaluate_behavioural_predictions.nb.html output/06_data_statistics.nb.html

clean:
	rm -f output/*


# Data formatting
data/formatted_Grandes_Lignes.fst data/formatted_Stepping_Stones.fst: /Users/maarten/Documents/projects/cold-start-afl-ii-full/data/noordhoff.sqlite
	Rscript scripts/01_format_data.R


# Rate of forgetting
data/rate_of_forgetting_Grandes_Lignes.fst data/rate_of_forgetting_Stepping_Stones.fst: data/formatted_Grandes_Lignes.fst data/formatted_Stepping_Stones.fst
	Rscript scripts/02_calculate_final_rate_of_forgetting.R


# Predictions
data/predictions/* : data/rate_of_forgetting_Grandes_Lignes.fst data/rate_of_forgetting_Stepping_Stones.fst
	Rscript scripts/03_predict_rate_of_forgetting.R


# Evaluate rate of forgetting predictions
output/04_evaluate_ROF_predictions.nb.html: scripts/04_evaluate_ROF_predictions.Rmd
	Rscript -e "rmarkdown::render('scripts/04_evaluate_ROF_predictions.Rmd', output_format = 'all')"
	rm -r output/04_evaluate_ROF_predictions_files
	mv scripts/04_evaluate_ROF_predictions.nb.html scripts/04_evaluate_ROF_predictions.html scripts/04_evaluate_ROF_predictions.md scripts/04_evaluate_ROF_predictions_files output


# Evaluate behavioural predictions
output/05_evaluate_behavioural_predictions.nb.html: scripts/05_evaluate_behavioural_predictions.Rmd
	Rscript -e "rmarkdown::render('scripts/05_evaluate_behavioural_predictions.Rmd', output_format = 'all')"
	rm -r output/05_evaluate_behavioural_predictions_files
	mv scripts/05_evaluate_behavioural_predictions.nb.html scripts/05_evaluate_behavioural_predictions.html scripts/05_evaluate_behavioural_predictions.md scripts/05_evaluate_behavioural_predictions_files output


# Extra data statistics
output/06_data_statistics.nb.html: scripts/06_data_statistics.Rmd
	Rscript -e "rmarkdown::render('scripts/06_data_statistics.Rmd', output_format = 'all')"
	mv scripts/06_data_statistics.nb.html scripts/06_data_statistics.html scripts/06_data_statistics.md scripts/06_data_statistics_files output