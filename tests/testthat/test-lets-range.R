
test_that("lets.range checks GEOS compatibility", {
  
  check_geos <- letsR:::.check_lets_range_geos
  
  expect_error(
    check_geos("3.10.9"),
    "requires GEOS >= 3.11.0",
    fixed = TRUE
  )
  
  expect_error(
    check_geos("3.9.0"),
    "requires GEOS >= 3.11.0",
    fixed = TRUE
  )
  
  expect_error(
    check_geos(NA_character_),
    "Cannot determine GEOS version"
  )
  
  expect_silent(check_geos("3.11.0"))
  expect_silent(check_geos("3.12.0"))
  expect_silent(check_geos("3.13.0"))
  
  expect_silent(
    check_geos("3.12.1-CAPI-1.18.1")
  )
})


test_that("lets.range estimates valid geographic ranges", {
  
  geos <- sf::sf_extSoftVersion()[["GEOS"]]
  geos_num <- sub(
    "^([0-9]+\\.[0-9]+\\.[0-9]+).*",
    "\\1",
    geos
  )
  
  skip_if(
    is.na(geos) ||
      utils::compareVersion(geos_num, "3.11.0") < 0L,
    "Requires GEOS >= 3.11.0"
  )
  
  occs <- data.frame(
    scientificName = rep("Species_alpha", 3),
    decimalLongitude = c(0, 0.1, 0.2),
    decimalLatitude = c(0, 0.1, 0.2)
  )
  
  mask <- sf::st_as_sf(
    data.frame(
      id = 1,
      wkt = paste0(
        "POLYGON ((-2 -2, 2 -2, 2 2, ",
        "-2 2, -2 -2))"
      )
    ),
    wkt = "wkt",
    crs = 4326
  )
  
  result <- letsR::lets.range(
    occs,
    crop_shp = mask,
    eps_m = 100000,
    min_pts = 1,
    buf = 10000,
    plot_map = FALSE
  )
  
  expect_s3_class(result, "sf")
  expect_gt(nrow(result), 0)
  expect_equal(sf::st_crs(result)$epsg, 4326)
  
  expect_true(
    all(c("binomial", "cluster") %in% names(result))
  )
  
  expect_equal(
    unique(result$binomial),
    "Species_alpha"
  )
  
  expect_true(
    all(sf::st_is_valid(result))
  )
  
  expect_true(
    all(
      sf::st_geometry_type(result) == "MULTIPOLYGON"
    )
  )
})


test_that("lets.range validates occurrence data", {
  
  invalid_occs <- data.frame(
    longitude = 0,
    latitude = 0
  )
  
  mask <- sf::st_sfc(
    sf::st_polygon(
      list(
        matrix(
          c(-2, -2, 2, -2, 2, 2,
            -2, 2, -2, -2),
          ncol = 2,
          byrow = TRUE
        )
      )
    ),
    crs = 4326
  )
  
  expect_error(
    letsR::lets.range(
      invalid_occs,
      crop_shp = mask,
      plot_map = FALSE
    ),
    "occs must contain"
  )
})
