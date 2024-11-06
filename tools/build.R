
# build package -----------------------------------------------------------

usethis::use_mit_license(copyright_holder = "Akira Terui")
usethis::use_roxygen_md()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)


# check syntax ------------------------------------------------------------

lintr::lint_package()

