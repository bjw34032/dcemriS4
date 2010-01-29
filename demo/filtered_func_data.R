ffd <- readNIfTI(system.file("nifti/filtered_func_data.nii.gz", package="dcemriS4"))
par(bg="black")
image(ffd)
orthographic(ffd)
