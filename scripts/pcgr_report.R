# Load the necessary library
library(pcgrr)

# Path to the configuration file
yaml_fname <- "pcgr/PCGR.pcgr.grch38.conf.yaml"

# Function to define log layout for log4r
my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - pcgr-report-generation - ", level, " - ", ..., "\n", collapse = "")
}

# Create a logger with INFO level and console appender
log4r_logger <- log4r::logger(threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))

# Set the logger for the package
options("PCGRR_LOG4R_LOGGER" = log4r_logger)

# Generate the report content
pcg_report <- pcgrr::generate_report(yaml_fname = yaml_fname)

# Write the report to a TSV file
pcgrr::write_report_tsv(report = pcg_report, output_type = 'snv_indel')

# Write the report to an Excel file
pcgrr::write_report_excel(report = pcg_report)
