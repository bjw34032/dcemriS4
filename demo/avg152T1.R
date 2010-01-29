avg152T1 <- readANALYZE(system.file("nifti/avg152T1.img.gz", package="dcemriS4"))
par(bg="black")
image(avg152T1)
orthographic(avg152T1)

