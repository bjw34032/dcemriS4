ffd <- readNIfTI(system.file("nifti/filtered_func_data.nii.gz",
                             package="dcemriS4"))
image(ffd)
orthographic(ffd)
