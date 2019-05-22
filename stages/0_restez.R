# Set up restez
library(restez)
restez_path_set(filepath = file.path('data'))
db_download(preselection = 13)  # "Other mammalian"
db_create()
