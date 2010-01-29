mniRL <- readNIfTI(system.file("nifti/avg152T1_RL_nifti.nii.gz", package="dcemriS4"))
par(bg="black")
image(mniRL)
orthographic(mniRL)
