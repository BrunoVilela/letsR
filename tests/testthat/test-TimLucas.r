context("Test letsR functions")

test_that("Basic usage works", {

	xy <- cbind(1:10, 1:10)
	species <- rep(c('Milvus milvus', 'Buteo buteo'), each = 5)

	PresAbMat <- lets.presab.points(xy, species)
	expect_equal(class(PresAbMat), "PresenceAbsence")
	
	expect_equal(dim(PresAbMat[[1]]), c(10, 4))
	expect_equal(sum(PresAbMat[[1]][,3]), 5)
	expect_equal(sum(PresAbMat[[1]][,4]), 5)

	expect_true(inherits(PresAbMat[[2]], "SpatRaster"))
	
	# Species list as factor not character
	species <- factor(rep(c('Milvus milvus', 'Buteo buteo'), each = 5))

	PresAbMat <- lets.presab.points(xy, species)
	expect_equal(class(PresAbMat), "PresenceAbsence")
	
	
	# # Check other projections
	# crsdif <- terra::crs("+proj=lcc +lat_1=48 +lat_2=33 +lon_0=-100 +ellps=WGS84")
	# PresAbMat <- lets.presab.points(xy, species, crs = crsdif)
	# expect_equal(class(PresAbMat), "PresenceAbsence")

	
})

# Doesn't work
test_that("Single species works", {

	xy <- cbind(1:10, 1:10)
	species <- rep('Buteo buteo', 10)

	PresAbMat <- lets.presab.points(xy, species)
	expect_equal(class(PresAbMat), "PresenceAbsence")

})


# Single records doesn't work. Line 87/98 of lets_presab_points.R
#  Probably not a big deal, noone uses single records.
test_that("Single records works", {

	xy <- cbind(1:3, 1:3)
	species <- c('Buteo buteo', 'Milvus milvus', 'Meles meles')

	PresAbMat <- lets.presab.points(xy, species, 
                                  crs=terra::crs("+proj=longlat +datum=WGS84"))
	expect_equal(class(PresAbMat), "PresenceAbsence")

})


test_that("Other parameters work", {

	xy <- cbind(1:10, 1:10)
	species <- rep(c('Milvus milvus', 'Buteo buteo'), each = 5)

	PresAbMat1 <- lets.presab.points(xy, species, resol = 5 )
	expect_equal(class(PresAbMat1), "PresenceAbsence")

	
	PresAbMat2 <- lets.presab.points(xy, species, remove.cells = FALSE )
	expect_equal(class(PresAbMat2), "PresenceAbsence")
	
	
	PresAbMat3 <- lets.presab.points(xy, species, remove.sp = FALSE )
	expect_equal(class(PresAbMat3), "PresenceAbsence")


	PresAbMat4 <- lets.presab.points(xy, species, show.matrix = TRUE )
	expect_true(is.matrix(PresAbMat4))


})
