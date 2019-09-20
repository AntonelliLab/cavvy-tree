# Libs ----
library(restez)

# Vars ----
gb_section <- '16' # Rodents
restez_path <- file.path('data', 'restez')

# Set up ----
if (!dir.exists(restez_path)) {
  dir.create(restez_path)
  restez_path_set(filepath = restez_path)
  db_download(preselection = gb_section)
  restez_connect()
  db_create(min_length = 150)
  restez_disconnect()
}

# Status check ----
restez_path_set(filepath = restez_path)
restez_connect()
restez_status()
restez_disconnect()
