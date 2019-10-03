# Libs ----
library(restez)

# Functions ----
check_restez <- function() {
  on.exit(restez_disconnect())
  restez_connect()
  restez_ready()
}

# Vars ----
gb_section <- '16' # Rodents
restez_path <- 'data'
restez_path_set(filepath = restez_path)

# Set up ----
if (!check_restez) {
  db_download(preselection = gb_section)
  restez_connect()
  db_create(min_length = 150)
  restez_disconnect()
}

# Status check ----
restez_connect()
restez_status()
restez_disconnect()
