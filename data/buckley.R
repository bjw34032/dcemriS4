buckley <- utils::read.table("buckley.txt", header=TRUE)
breast <- list(data = buckley[,4:16])
breast$fp <- c(.17, .37, .57, .77, .97, rep(.57, 8))
breast$ps <- c(rep(.33, 5), .01, .17, .49, .65, rep(.33, 4))
breast$vp <- c(rep(.06, 9), .0001, .03, .09, .12)
breast$ve <- rep(.45, 13)
breast$ktrans <- (1 - exp(-breast$ps/breast$fp)) * breast$fp
breast$kep <- breast$ktrans / breast$ve
meningioma <- list(data = buckley[,17:29])
meningioma$fp <- c(.4, .8, 1.2, 1.6, 2.0, rep(1.2, 8))
meningioma$ps <- c(rep(.34, 5), 0, .17, .51, .68, rep(.34, 4))
meningioma$vp <- c(rep(.08, 9), .0001, .04, .12, .16)
meningioma$ve <- rep(.4, 13)
meningioma$ktrans <- (1 - exp(-meningioma$ps/meningioma$fp)) * meningioma$fp
meningioma$kep <- meningioma$ktrans / meningioma$ve
